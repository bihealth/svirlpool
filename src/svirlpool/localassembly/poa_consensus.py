"""
In-memory POA consensus assembly.

An alternative to the ``lamassemble`` consensus step that never touches the
filesystem: abPOA (via ``pyabpoa``) builds a partial-order graph over the
trimmed reads of one cluster and emits a draft consensus, which is then
polished with a racon-equivalent windowed-POA pass driven by minimap2's
Python bindings (``mappy``).  Both libraries work on Python strings, so a
cluster is assembled without a single temporary file or subprocess.

The polishing pass mirrors what racon does: align every read of the cluster
back to the draft, cut the draft into windows, re-run POA per window with the
draft window as the graph backbone, and concatenate the window consensuses.
Running POA twice is not redundant -- the first pass threads raw, noisy reads
into a graph in input order, which mis-threads in tandem repeats, while the
second pass places every read by a proper minimap2 alignment first.
"""

import logging
from dataclasses import dataclass, field

import mappy as mp
import pyabpoa as pa

log = logging.getLogger(__name__)

# Defaults chosen to match racon's behaviour where an equivalent exists.
DEFAULT_WINDOW_SIZE = 500
DEFAULT_MIN_WINDOW_READS = 3
# Minimum fraction of a window a read must cover before it may vote on it.
DEFAULT_MIN_WINDOW_OVERLAP = 0.5
DEFAULT_MAX_POA_READS = 12
DEFAULT_POLISH_ROUNDS = 1
# Backbone extension: how far the layout may grow past the seed read.
DEFAULT_EXTENSION_ROUNDS = 6
DEFAULT_MIN_OVERHANG = 50
DEFAULT_MIN_EXTENSION_OVERLAP = 500
DEFAULT_MIN_EXTENSION_IDENTITY = 0.70
DEFAULT_END_TOLERANCE = 500

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# htslib CIGAR opcodes, as reported by mappy
_CIGAR_MATCH = {0, 7, 8}  # M, =, X -> consume reference and query
_CIGAR_INS = 1  # I -> consume query only
_CIGAR_DEL = {2, 3}  # D, N -> consume reference only
_CIGAR_SOFTCLIP = 4  # S -> consume query only


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class PoaAssemblyStats:
    """Bookkeeping for one cluster assembly, for logging and benchmarking."""

    n_reads: int = 0
    n_reads_in_poa: int = 0
    n_reads_flipped: int = 0
    n_reads_unoriented: int = 0
    draft_length: int = 0
    final_length: int = 0
    polish_rounds: int = 0
    windows_total: int = 0
    windows_polished: int = 0
    windows_rejected: int = 0
    backbone_from_poa: bool = False
    reads_mapped_per_round: list[int] = field(default_factory=list)


def _new_aligner(aln_mode: str = "g") -> pa.msa_aligner:
    """abPOA aligner with ONT-appropriate affine gap penalties.

    The scoring is abPOA's own long-read default (match 2 / mismatch 4,
    two-piece affine gaps 4,2 and 24,1), which is also minimap2's ``map-ont``
    scheme and the right choice over tandem repeats: the second, cheap
    extension piece lets a whole repeat unit be skipped without paying a
    linear price per base.
    """
    return pa.msa_aligner(
        aln_mode=aln_mode,
        match=2,
        mismatch=4,
        gap_open1=4,
        gap_ext1=2,
        gap_open2=24,
        gap_ext2=1,
    )


def _poa(
    seqs: list[str], aligner: pa.msa_aligner | None = None, aln_mode: str = "g"
) -> str | None:
    """Run abPOA over ``seqs`` and return the consensus, or None if it failed."""
    if len(seqs) == 0:
        return None
    if len(seqs) == 1:
        return seqs[0]
    try:
        res = (aligner or _new_aligner(aln_mode)).msa(
            seqs, out_cons=True, out_msa=False, max_n_cons=1
        )
    except Exception as e:  # pyabpoa raises bare exceptions on malformed input
        log.warning(f"abPOA failed on {len(seqs)} sequences: {e}")
        return None
    if not res.cons_seq:
        return None
    cons = res.cons_seq[0]
    return cons if len(cons) > 0 else None


def orient_reads(
    sequences: dict[str, str], stats: PoaAssemblyStats | None = None
) -> dict[str, str]:
    """Flip reads onto a common strand so POA can thread them into one graph.

    POA has no notion of reverse complements: a read on the opposite strand
    would be threaded as an unrelated sequence and wreck the graph.  The read
    closest to the median length is used as the strand reference -- a median
    read is far less likely to be a chimera than the longest one.

    Reads that do not map to the reference read are passed through unchanged;
    dropping them here would silently discard a divergent haplotype, and the
    polishing pass ignores whatever it cannot align anyway.
    """
    if len(sequences) < 2:
        return dict(sequences)

    by_length = sorted(sequences.items(), key=lambda kv: len(kv[1]))
    ref_name, ref_seq = by_length[len(by_length) // 2]

    try:
        aligner = mp.Aligner(seq=ref_seq, preset="map-ont", best_n=1)
    except Exception as e:
        log.warning(f"mappy failed to index the strand reference read: {e}")
        return dict(sequences)
    if not aligner:
        return dict(sequences)

    oriented: dict[str, str] = {}
    for name, seq in sequences.items():
        if name == ref_name:
            oriented[name] = seq
            continue
        hit = next((h for h in aligner.map(seq) if h.is_primary), None)
        if hit is None:
            oriented[name] = seq
            if stats is not None:
                stats.n_reads_unoriented += 1
            continue
        if hit.strand == -1:
            oriented[name] = _revcomp(seq)
            if stats is not None:
                stats.n_reads_flipped += 1
        else:
            oriented[name] = seq
    return oriented


def _select_poa_reads(sequences: dict[str, str], max_reads: int) -> list[str]:
    """Pick the reads that seed the draft graph, median-length read first.

    abPOA builds its graph from the first sequence and aligns the rest into
    it, so the first sequence is the backbone.  A median-length read is the
    safest backbone: the longest read of a cluster is disproportionately
    likely to be a chimera or to carry an over-expanded repeat array.
    """
    ordered = sorted(sequences.values(), key=len)
    if len(ordered) == 0:
        return []
    median_idx = len(ordered) // 2
    backbone = ordered[median_idx]
    rest = ordered[:median_idx] + ordered[median_idx + 1 :]
    if max_reads > 0 and len(rest) > max_reads - 1:
        # Keep the reads closest in length to the backbone: they align into
        # the graph most cleanly and dominate the consensus anyway.
        rest = sorted(rest, key=lambda s: abs(len(s) - len(backbone)))[: max_reads - 1]
    return [backbone] + rest


def _window_boundaries(length: int, window_size: int) -> list[int]:
    """Window edges over [0, length], merging away a runt final window."""
    if length <= window_size:
        return [0, length]
    bounds = list(range(0, length, window_size))
    if length - bounds[-1] < window_size // 2:
        bounds.pop()
    bounds.append(length)
    return bounds


def _ref_to_query_anchors(
    hit: mp.Alignment, boundaries: list[int], query_length: int
) -> dict[int, int]:
    """Map backbone positions to query offsets by walking one alignment's CIGAR.

    Query offsets are returned in *alignment* orientation, i.e. against the
    reverse complement of the read for a reverse-strand hit.  mappy reports
    ``q_st``/``q_en`` on the original read (as PAF does) while its CIGAR walks
    the query in alignment orientation, so the start offset is mirrored for
    reverse-strand hits.

    Only boundaries the alignment actually spans are returned.  A boundary
    falling inside a deletion is anchored to the query position where the
    deletion starts, which keeps neighbouring windows contiguous.
    """
    anchors: dict[int, int] = {}
    ref_pos = hit.r_st
    query_pos = hit.q_st if hit.strand == 1 else query_length - hit.q_en

    idx = 0
    while idx < len(boundaries) and boundaries[idx] < ref_pos:
        idx += 1

    for length, op in hit.cigar:
        if idx >= len(boundaries):
            break
        if op in _CIGAR_MATCH:
            while idx < len(boundaries) and boundaries[idx] < ref_pos + length:
                anchors[boundaries[idx]] = query_pos + (boundaries[idx] - ref_pos)
                idx += 1
            ref_pos += length
            query_pos += length
        elif op == _CIGAR_INS or op == _CIGAR_SOFTCLIP:
            query_pos += length
        elif op in _CIGAR_DEL:
            while idx < len(boundaries) and boundaries[idx] < ref_pos + length:
                anchors[boundaries[idx]] = query_pos
                idx += 1
            ref_pos += length
        # hard clips consume neither reference nor query

    if idx < len(boundaries) and boundaries[idx] == ref_pos:
        anchors[boundaries[idx]] = query_pos
    return anchors


def _primary_hit(aligner: mp.Aligner, seq: str):
    return next((h for h in aligner.map(seq) if h.is_primary), None)


def layout_score(consensus: str, sequences: dict[str, str]) -> float:
    """How many read bases the consensus explains, over all reads.

    This is the objective every layout decision is judged against.  A longer
    consensus only helps if the extra sequence lets more read bases align;
    a chimeric join makes reads align in pieces and the score falls.
    """
    if len(consensus) == 0:
        return 0.0
    try:
        aligner = mp.Aligner(seq=consensus, preset="map-ont", best_n=1)
    except Exception:
        return 0.0
    if not aligner:
        return 0.0
    total = 0
    for seq in sequences.values():
        hit = _primary_hit(aligner, seq)
        if hit is not None:
            total += hit.q_en - hit.q_st
    return float(total)


def select_backbone(sequences: dict[str, str], n_candidates: int = 5) -> str:
    """Pick the read that explains the most read bases, not simply the longest.

    The longest read of a cluster is often the one carrying an artefact -- an
    untrimmed adapter, a chimeric join, an over-expanded repeat array -- and
    seeding the layout with it propagates that artefact into the consensus.
    Scoring the few longest candidates against the whole cluster costs one
    minimap2 index each and picks the read the cluster actually agrees with.
    """
    ordered = sorted(sequences.values(), key=len, reverse=True)
    candidates = ordered[: max(1, n_candidates)]
    if len(candidates) == 1:
        return candidates[0]
    return max(candidates, key=lambda c: layout_score(c, sequences))


def extend_backbone(
    backbone: str,
    sequences: dict[str, str],
    max_rounds: int = DEFAULT_EXTENSION_ROUNDS,
    min_overhang: int = DEFAULT_MIN_OVERHANG,
    min_overlap: int = DEFAULT_MIN_EXTENSION_OVERLAP,
    min_overlap_identity: float = DEFAULT_MIN_EXTENSION_IDENTITY,
    end_tolerance: int = DEFAULT_END_TOLERANCE,
) -> str:
    """Grow the backbone outwards along reads that hang off its ends.

    A cluster's reads tile a locus; no single read need span it.  Polishing
    alone can never reach past the backbone read, so the locus gets truncated
    to the longest read -- which is how a backbone-only assembly loses the
    part of an insertion that lamassemble's overlap layout recovers.

    Each round splices on the read with the largest overhang at either end,
    the read replacing the backbone across the overlap.  Over a tandem repeat
    an overlap can be genuine yet out of phase by whole repeat units, which
    welds on a shifted copy and leaves a consensus that its own reads align to
    only in pieces.  Two guards make that self-correcting: a candidate join
    needs a long, reasonably similar overlap, and every accepted round has to
    raise :func:`layout_score`.  The first round that does not is discarded
    and extension stops.
    """
    best_score = layout_score(backbone, sequences)
    for _ in range(max_rounds):
        try:
            aligner = mp.Aligner(seq=backbone, preset="map-ont", best_n=1)
        except Exception as e:
            log.warning(f"mappy failed to index the backbone during extension: {e}")
            return backbone
        if not aligner:
            return backbone

        best_left: tuple[int, str] = (0, "")
        best_right: tuple[int, str] = (0, "")
        for seq in sequences.values():
            hit = _primary_hit(aligner, seq)
            if hit is None or hit.r_en - hit.r_st < min_overlap:
                continue
            if (1.0 - hit.NM / max(1, hit.blen)) < min_overlap_identity:
                continue
            query = seq if hit.strand == 1 else _revcomp(seq)
            q_st = hit.q_st if hit.strand == 1 else len(seq) - hit.q_en
            q_en = hit.q_en if hit.strand == 1 else len(seq) - hit.q_st

            if hit.r_st <= end_tolerance and q_st > max(min_overhang, best_left[0]):
                best_left = (q_st, query[:q_en] + backbone[hit.r_en :])
            right_overhang = len(seq) - q_en
            if hit.r_en >= len(backbone) - end_tolerance and right_overhang > max(
                min_overhang, best_right[0]
            ):
                best_right = (right_overhang, backbone[: hit.r_st] + query[q_st:])

        if best_left[0] >= best_right[0] and best_left[0] > 0:
            candidate = best_left[1]
        elif best_right[0] > 0:
            candidate = best_right[1]
        else:
            break

        score = layout_score(candidate, sequences)
        if score <= best_score:
            # The join explained no extra read bases, so it was noise or a
            # repeat-phase slip.  Keep the shorter, coherent backbone.
            break
        backbone, best_score = candidate, score
    return backbone


def polish_windowed_poa(
    backbone: str,
    sequences: dict[str, str],
    window_size: int = DEFAULT_WINDOW_SIZE,
    min_window_reads: int = DEFAULT_MIN_WINDOW_READS,
    min_window_overlap: float = DEFAULT_MIN_WINDOW_OVERLAP,
    stats: PoaAssemblyStats | None = None,
) -> str:
    """Polish ``backbone`` with the cluster reads, racon-style, entirely in memory.

    Every read is aligned to the backbone with minimap2, the backbone is cut
    into windows, and each window is re-assembled by POA from the backbone
    window plus the read fragments covering it.  Windows without enough
    fragments are emitted unchanged, so low-coverage stretches are carried
    over rather than invented.

    A read that covers only part of a window still contributes: its fragment
    is padded with the backbone's own sequence out to the window edges.
    Requiring reads to span a whole window would silence every read shorter
    than ``window_size``, which in a cut-read cluster is most of them.  The
    padding makes the backbone win wherever a read does not reach, which is
    the conservative direction.
    """
    if len(backbone) == 0 or len(sequences) == 0:
        return backbone

    try:
        aligner = mp.Aligner(seq=backbone, preset="map-ont", best_n=1)
    except Exception as e:
        log.warning(f"mappy failed to index the draft consensus: {e}")
        return backbone
    if not aligner:
        return backbone

    boundaries = _window_boundaries(len(backbone), window_size)
    n_windows = len(boundaries) - 1
    fragments: list[list[str]] = [[] for _ in range(n_windows)]

    n_mapped = 0
    for seq in sequences.values():
        hit = next((h for h in aligner.map(seq) if h.is_primary), None)
        if hit is None:
            continue
        n_mapped += 1
        query = seq if hit.strand == 1 else _revcomp(seq)
        # The read's own alignment ends are anchor positions too, so that a
        # read covering part of a window can be cut at exactly that point.
        positions = sorted(set(boundaries) | {hit.r_st, hit.r_en})
        anchors = _ref_to_query_anchors(hit, positions, len(seq))
        for w in range(n_windows):
            ws, we = boundaries[w], boundaries[w + 1]
            lo, hi = max(ws, hit.r_st), min(we, hit.r_en)
            if hi - lo < min_window_overlap * (we - ws):
                continue
            start, end = anchors.get(lo), anchors.get(hi)
            if start is None or end is None or end <= start:
                continue
            fragments[w].append(backbone[ws:lo] + query[start:end] + backbone[hi:we])

    if stats is not None:
        stats.reads_mapped_per_round.append(n_mapped)

    poa_aligner = _new_aligner()
    pieces: list[str] = []
    for w in range(n_windows):
        window_seq = backbone[boundaries[w] : boundaries[w + 1]]
        frags = fragments[w]
        if len(frags) < min_window_reads:
            pieces.append(window_seq)
            continue
        # Fragments wildly off the window length are mis-anchored, not variant;
        # letting them into the graph is how a polished window blows up.
        lo, hi = len(window_seq) / 3.0, len(window_seq) * 3.0
        frags = [f for f in frags if lo <= len(f) <= hi]
        if len(frags) < min_window_reads:
            pieces.append(window_seq)
            continue

        polished = _poa([window_seq] + frags, aligner=poa_aligner)
        if polished is None or not (lo <= len(polished) <= hi):
            pieces.append(window_seq)
            if stats is not None:
                stats.windows_rejected += 1
            continue
        pieces.append(polished)
        if stats is not None:
            stats.windows_polished += 1

    if stats is not None:
        stats.windows_total += n_windows
    return "".join(pieces)


def _polish_rounds(
    backbone: str,
    oriented: dict[str, str],
    rounds: int,
    window_size: int,
    min_window_reads: int,
    stats: PoaAssemblyStats | None,
) -> str:
    consensus = backbone
    for _ in range(max(0, rounds)):
        polished = polish_windowed_poa(
            backbone=consensus,
            sequences=oriented,
            window_size=window_size,
            min_window_reads=min_window_reads,
            stats=stats,
        )
        if len(polished) == 0:
            break
        consensus = polished
        if stats is not None:
            stats.polish_rounds += 1
    return consensus


def assemble_consensus_poa(
    sequences: dict[str, str],
    name: str = "",
    polish_rounds: int = DEFAULT_POLISH_ROUNDS,
    window_size: int = DEFAULT_WINDOW_SIZE,
    min_window_reads: int = DEFAULT_MIN_WINDOW_READS,
    max_poa_reads: int = DEFAULT_MAX_POA_READS,
    extension_rounds: int = DEFAULT_EXTENSION_ROUNDS,
    aln_mode: str = "g",
    stats: PoaAssemblyStats | None = None,
) -> str | None:
    """Assemble one cluster's reads into a consensus without touching disk.

    Two backbones are built and the one that explains more read bases wins:

    * the abPOA consensus over the whole cluster, which is the heaviest path
      through the graph -- right when the reads stack, but in a cluster of cut
      reads that tile a locus it follows the depth peak and collapses the
      locus to it;
    * the read that hosts the cluster best, extended along reads hanging off
      its ends, which follows the tiling but cannot represent a region no
      single read reaches into.

    Neither is reliably better, and which one wins is a property of the
    cluster rather than of the dataset, so both are built and
    :func:`layout_score` decides.  The winner is then polished window by
    window, where every fragment really does describe the same stretch of
    sequence and POA is on solid ground.

    Args:
        sequences: Read name -> sequence for the reads of a single cluster.
        name: Cluster name, used only for log messages.
        polish_rounds: Windowed-POA passes over the chosen backbone.
        window_size: Backbone window size for polishing.
        min_window_reads: Fragments required before a window is re-assembled.
        max_poa_reads: Cap on reads entering the whole-cluster draft graph;
            0 means all.  Polishing always uses every read regardless.
        extension_rounds: Rounds of backbone extension along overhanging reads.
        aln_mode: abPOA alignment mode for the whole-cluster draft.
        stats: Optional object collecting per-cluster bookkeeping.

    Returns:
        The consensus sequence, or None if no consensus could be built.
    """
    clean = {n: s.upper() for n, s in sequences.items() if len(s) > 0}
    if len(clean) == 0:
        log.warning(f"No usable reads for POA consensus '{name}'.")
        return None
    if stats is not None:
        stats.n_reads = len(clean)
    if len(clean) == 1:
        only = next(iter(clean.values()))
        if stats is not None:
            stats.draft_length = stats.final_length = len(only)
        return only

    oriented = orient_reads(clean, stats=stats)

    candidates: list[str] = []

    poa_reads = _select_poa_reads(oriented, max_poa_reads)
    if stats is not None:
        stats.n_reads_in_poa = len(poa_reads)
    draft = _poa(poa_reads, aln_mode=aln_mode)
    if draft is not None:
        candidates.append(draft)

    laid_out = extend_backbone(
        select_backbone(oriented), oriented, max_rounds=extension_rounds
    )
    if len(laid_out) > 0:
        candidates.append(laid_out)

    if len(candidates) == 0:
        log.warning(f"No backbone could be built for POA consensus '{name}'.")
        return None

    backbone = max(candidates, key=lambda c: layout_score(c, oriented))
    if stats is not None:
        stats.draft_length = len(backbone)
        stats.backbone_from_poa = backbone is draft

    consensus = _polish_rounds(
        backbone, oriented, polish_rounds, window_size, min_window_reads, stats
    )
    if stats is not None:
        stats.final_length = len(consensus)
    return consensus if len(consensus) > 0 else None
