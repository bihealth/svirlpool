"""Cross-sample merge of SVcomposites by aligning the consensus cores (prototype).

Switched on with ``svirlpool sv-calling --core-alignment-merge align|or``; the
default ``off`` leaves the merge untouched.

Idea: every SVpattern sits on a consensus sequence (``consensusID``,
``read_start``/``read_end`` on the padded consensus). Two consensuses whose
reference placements overlap (different samples, or two haplotypes of one
sample) are aligned with minimap2, core against core. An SV ``x`` on consensus
A is *present* on consensus B when the A-B alignment across ``x`` and its flanks
is near-identical (few indel bases); its interval then projects through the
alignment onto B. ``x`` and an SV ``y`` of B are the same allele when one of
them is present on the other's consensus and its projection overlaps the other.
This judges "same allele" by sequence, so

- tandem-repeat length variants are compared by the length difference of the
  haplotypes, not by the sizes of the reference-alignment indels;
- one haplotype decomposed differently (INS 149 + 199 + 597 vs INS 944) links
  every piece to the single event;
- two different events near each other (or sharing a repeatID) are not linked
  unless the haplotypes say so.

Pairs that the alignment cannot judge (an SV without aligned flanks, cores that
do not align) fall back to the existing pairwise rule. See
``experiments/core_ava_merge/DESIGN.md``.
"""

from __future__ import annotations

import logging
import os
import pickle
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Callable
from pathlib import Path

import attrs

from ..localassembly import SVpatterns
from ..util.datastructures import UnionFind
from .SVcomposite import SVcomposite

log = logging.getLogger(__name__)

#: "off" (default, the existing merge), "align" (the alignment decides every
#: pair it can judge, the existing rule the others) or "or" (a pair merges when
#: either the alignment or the existing rule says so). Set once from the CLI
#: before any worker process is forked.
MODE: str = "off"
MODES = ("off", "align", "or")
#: when set, every judged composite pair is written to <dir>/pairs.<chr>.jsonl
PAIR_LOG_DIR: str | None = None

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@attrs.define(frozen=True)
class CoreAlnParams:
    #: bases of flank on each side of an SV that must agree between the two
    #: haplotypes (clipped to the aligned part of the two cores)
    flank: int = 150
    #: at least this many aligned bases on each side of the SV, else this side
    #: of the pair cannot be judged
    min_flank: int = 30
    #: tolerance for indel bases between the two haplotypes inside the window:
    #: max(abs_tol, frac_tol * SV size)
    abs_tol: int = 10
    frac_tol: float = 0.1
    #: an inserted sequence may differ by at most this fraction of mismatches
    max_mismatch_frac: float = 0.2
    #: slack around the projected interval when looking for the partner SV
    link_slack: int = 50
    minimap2: str = ""
    preset: str = "asm20"
    #: minimizers occurring more often than max(10, occ_per_core * cores in the
    #: group) are not used as anchors. A sequence unique to the locus occurs once
    #: per core; without a cap of this order the cores of satellite loci produce
    #: millions of anchors (minutes per group).
    occ_per_core: int = 3
    #: seconds per minimap2 call; a group that takes longer is not aligned (its
    #: pairs fall back to the pairwise rule)
    timeout: float = 10.0
    #: groups with more core bases than this are not aligned (fallback)
    max_group_bp: int = 300_000
    #: band around the diagonal implied by the reference placements of two cores
    #: (plus the indel bases of their SVs) in which alignment blocks are kept
    diagonal_band: int = 1000


PARAMS: CoreAlnParams = CoreAlnParams()


def minimap2_binary(params: CoreAlnParams) -> str:
    if params.minimap2:
        return params.minimap2
    found = shutil.which("minimap2")
    if found:
        return found
    local = Path(sys.executable).parent / "minimap2"
    if local.exists():
        return str(local)
    raise FileNotFoundError("minimap2 not found on PATH or next to the interpreter")


# --------------------------------------------------------------------------- #
# core sequences from the per-sample DBs
# --------------------------------------------------------------------------- #


def parse_core_description(description: str) -> tuple[int, int]:
    """core=S-E from a consensus description (positions on the padded consensus)."""
    for token in description.split():
        if token.startswith("core="):
            s, e = token[5:].split("-")
            return int(s), int(e)
    raise ValueError(f"no core= in consensus description {description!r}")


def _load_cores_job(args) -> tuple[str, dict[str, tuple[int, str]], float]:
    samplename, db_path, consensusIDs = args
    t0 = time.perf_counter()
    return samplename, load_core_sequences(db_path, consensusIDs), time.perf_counter() - t0


def load_cores_for_composites(
    composites: list[SVcomposite],
    db_paths: list[Path | str],
    samplenames: list[str],
    threads: int = 1,
) -> dict[tuple[str, str], tuple[int, str]]:
    """Cores of every consensus that carries an INS or DEL composite, all samples."""
    from concurrent.futures import ProcessPoolExecutor

    wanted: dict[str, set[str]] = defaultdict(set)
    for svc in composites:
        if issubclass(svc.sv_type, (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion)):
            for p in svc.svPatterns:
                wanted[p.samplename].add(p.consensusID)
    jobs = [(s, db, wanted.get(s, set())) for s, db in zip(samplenames, db_paths, strict=True)]
    out: dict[tuple[str, str], tuple[int, str]] = {}
    with ProcessPoolExecutor(max_workers=max(1, min(threads, len(jobs)))) as ex:
        for samplename, cores, seconds in ex.map(_load_cores_job, jobs):
            log.info(
                "CORE_ALN_MERGE|CORES|sample=%s|n=%d|bp=%d|seconds=%.1f",
                samplename,
                len(cores),
                sum(len(s) for _, s in cores.values()),
                seconds,
            )
            for cid, value in cores.items():
                out[(samplename, cid)] = value
    return out


def load_core_sequences(
    db_path: Path | str, consensusIDs: set[str] | None = None
) -> dict[str, tuple[int, str]]:
    """{consensusID: (core_start on the padded consensus, core sequence)}.

    Streams the consensus_sequences table (the padded sequences are ~75 kb each,
    ~5 Gb per 30x sample) and keeps only the cores of the wanted consensuses.
    """
    out: dict[str, tuple[int, str]] = {}
    uri = f"file:{db_path}?mode=ro"
    conn = sqlite3.connect(uri, uri=True, timeout=60)
    try:
        cur = conn.execute(
            "SELECT consensusID, sequence, description FROM consensus_sequences"
        )
        for consensusID, blob, description in cur:
            if consensusIDs is not None and consensusID not in consensusIDs:
                continue
            s, e = parse_core_description(description or "")
            padded: str = pickle.loads(blob)
            out[consensusID] = (s, padded[s:e].upper())
    finally:
        conn.close()
    return out


# --------------------------------------------------------------------------- #
# where a composite sits on its consensus core
# --------------------------------------------------------------------------- #


@attrs.define
class CoreSV:
    """A composite's interval on its consensus core, reference-forward oriented.

    [start, end) is where the event sits (the inserted bases of an insertion, a
    point for a deletion); [amb_start, amb_end) extends it over the bases across
    which the same event could be placed (an indel inside a tandem repeat can
    shift by whole periods without changing the haplotype).
    """

    key: tuple[str, str]  # (samplename, consensusID)
    start: int
    end: int
    size: int
    amb_start: int = -1
    amb_end: int = -1
    #: the composite's reference interval (chr dropped; one chromosome per group)
    ref_start: int = -1
    ref_end: int = -1

    def __attrs_post_init__(self):
        if self.amb_start < 0:
            self.amb_start = self.start
        if self.amb_end < 0:
            self.amb_end = self.end


#: cap of the placement-ambiguity extension on each side (bp)
MAX_SHIFT: int = 1000


def insertion_shift_range(seq: str, s: int, e: int, cap: int = MAX_SHIFT) -> tuple[int, int]:
    """[s - L, e + R]: how far the inserted bases seq[s:e] can slide in seq."""
    n = len(seq)
    r = 0
    while r < cap and e + r < n and seq[s + r] == seq[e + r]:
        r += 1
    left = 0
    while left < cap and s - 1 - left >= 0 and seq[s - 1 - left] == seq[e - 1 - left]:
        left += 1
    return s - left, e + r


def deletion_shift_range(
    seq: str, s: int, deleted: str, cap: int = MAX_SHIFT
) -> tuple[int, int]:
    """[s - L, s + R]: how far a deletion of `deleted` at point s can slide in seq."""
    m = len(deleted)
    if m == 0:
        return s, s
    n = len(seq)
    r = 0
    while r < cap and s + r < n and seq[s + r] == deleted[r % m]:
        r += 1
    left = 0
    while left < cap and s - 1 - left >= 0 and seq[s - 1 - left] == deleted[m - 1 - (left % m)]:
        left += 1
    return s - left, s + r


def composite_consensus_key(svc: SVcomposite) -> tuple[str, str] | None:
    """(samplename, consensusID) if all patterns come from one consensus, else None."""
    keys = {(p.samplename, p.consensusID) for p in svc.svPatterns}
    return next(iter(keys)) if len(keys) == 1 else None


def composite_strand(svc: SVcomposite) -> bool | None:
    """aln_is_reverse shared by all primitives of the composite, else None."""
    strands = {
        prim.aln_is_reverse for p in svc.svPatterns for prim in p.SVprimitives
    }
    return next(iter(strands)) if len(strands) == 1 else None


def core_interval_of_composite(
    svc: SVcomposite, core_start: int, oriented_core: str, reverse: bool
) -> CoreSV | None:
    """The composite's interval on its core, oriented to the reference strand.

    read_start/read_end of the primitives are forward positions on the padded
    consensus; ``reverse`` says whether the core was reverse complemented to
    match the reference strand (the consensus aligns in reverse).
    """
    key = composite_consensus_key(svc)
    if key is None or composite_strand(svc) is not reverse:
        return None
    core_len = len(oriented_core)
    s = min(p.read_start for p in svc.svPatterns) - core_start
    e = max(p.read_end for p in svc.svPatterns) - core_start
    if s > e:
        s, e = e, s
    if s < 0 or e > core_len:
        return None
    if reverse:
        s, e = core_len - e, core_len - s
    amb_s, amb_e = s, e
    if len(svc.svPatterns) == 1:
        pattern = svc.svPatterns[0]
        if isinstance(pattern, SVpatterns.SVpatternInsertion) and e > s:
            amb_s, amb_e = insertion_shift_range(oriented_core, s, e)
        elif isinstance(pattern, SVpatterns.SVpatternDeletion) and e - s <= 1:
            deleted = pattern.get_sequence()
            if deleted:
                amb_s, amb_e = deletion_shift_range(oriented_core, s, deleted.upper())
                amb_e = max(amb_e, e)
    return CoreSV(
        key=key,
        start=s,
        end=e,
        size=abs(svc.get_size()),
        amb_start=amb_s,
        amb_end=amb_e,
        ref_start=svc.ref_start[1],
        ref_end=svc.ref_end[1],
    )


# --------------------------------------------------------------------------- #
# minimap2 alignments of cores
# --------------------------------------------------------------------------- #


@attrs.define
class Block:
    """One minimap2 alignment block, + strand, with its CIGAR (=/X/I/D)."""

    qs: int
    qe: int
    ts: int
    te: int
    score: int
    cigar: list[tuple[int, str]]


def parse_cigar(cg: str) -> list[tuple[int, str]]:
    out: list[tuple[int, str]] = []
    n = 0
    for ch in cg:
        if ch.isdigit():
            n = n * 10 + ord(ch) - 48
        else:
            out.append((n, ch))
            n = 0
    return out


def _parse_paf(stdout: str, blocks: dict[tuple[str, str], list[Block]]) -> None:
    for line in stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[4] != "+":
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        if cg is None:
            continue
        score = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), int(f[9]))
        blocks[(f[0], f[5])].append(
            Block(
                qs=int(f[2]),
                qe=int(f[3]),
                ts=int(f[7]),
                te=int(f[8]),
                score=score,
                cigar=parse_cigar(cg),
            )
        )


def run_minimap2(
    seqs: dict[str, str], params: CoreAlnParams, tmp_dir: str | None = None
) -> dict[tuple[str, str], list[Block]] | None:
    """Pairwise alignments of the given sequences; {(query, target): [+ blocks]}.

    One call per target with the later sequences as queries and
    --secondary=no: each pair yields its primary alignment and the
    supplementary ones (the flanks of a large indel come back as separate
    blocks). All-vs-all in one call (-X) reports every secondary chain, which in
    tandem repeats are thousands per pair (minutes per group). None when the
    group exceeds params.timeout seconds in total.
    """
    names = list(seqs)
    if len(names) < 2:
        return {}
    occ = max(10, params.occ_per_core * len(names))
    blocks: dict[tuple[str, str], list[Block]] = defaultdict(list)
    t_end = time.perf_counter() + params.timeout
    with tempfile.TemporaryDirectory(dir=tmp_dir) as td:
        tf = os.path.join(td, "t.fa")
        qf = os.path.join(td, "q.fa")
        for i, tn in enumerate(names[:-1]):
            with open(tf, "w") as f:
                f.write(f">{tn}\n{seqs[tn]}\n")
            with open(qf, "w") as f:
                for qn in names[i + 1 :]:
                    f.write(f">{qn}\n{seqs[qn]}\n")
            cmd = [
                minimap2_binary(params),
                "-x",
                params.preset,
                "-c",
                "--eqx",
                "--secondary=no",
                "-U",
                f"{occ},{occ}",
                "-t",
                "1",
                tf,
                qf,
            ]
            remaining = t_end - time.perf_counter()
            if remaining <= 0:
                return None
            try:
                res = subprocess.run(
                    cmd, capture_output=True, text=True, check=True, timeout=remaining
                )
            except subprocess.TimeoutExpired:
                return None
            _parse_paf(res.stdout, blocks)
    return dict(blocks)


def best_colinear_chain(blocks: list[Block]) -> list[Block]:
    """Highest-scoring chain of blocks increasing on both sequences, no overlaps."""
    if not blocks:
        return []
    bl = sorted(blocks, key=lambda b: (b.qs, b.ts))
    best = [b.score for b in bl]
    prev = [-1] * len(bl)
    for i, bi in enumerate(bl):
        for j in range(i):
            bj = bl[j]
            if bj.qe <= bi.qs and bj.te <= bi.ts and best[j] + bi.score > best[i]:
                best[i] = best[j] + bi.score
                prev[i] = j
    i = max(range(len(bl)), key=lambda k: best[k])
    chain = []
    while i >= 0:
        chain.append(bl[i])
        i = prev[i]
    return chain[::-1]


@attrs.define
class PairAlignment:
    """A co-linear alignment of core q against core t.

    runs: (q, t, length) of each gap-free =/X stretch; mm_q/mm_t: mismatch
    positions; gaps: (q, t, qlen, tlen) of each I (qlen > 0), D (tlen > 0) and
    each stretch between two blocks (both may be > 0: unaligned on both).
    """

    q0: int
    q1: int
    t0: int
    t1: int
    runs: list[tuple[int, int, int]]
    mm_q: list[int]
    mm_t: list[int]
    gaps: list[tuple[int, int, int, int]]

    @classmethod
    def from_chain(cls, chain: list[Block]) -> PairAlignment:
        runs: list[tuple[int, int, int]] = []
        mm_q: list[int] = []
        mm_t: list[int] = []
        gaps: list[tuple[int, int, int, int]] = []
        prev: Block | None = None
        for b in chain:
            if prev is not None and (b.qs > prev.qe or b.ts > prev.te):
                gaps.append((prev.qe, prev.te, b.qs - prev.qe, b.ts - prev.te))
            qp, tp = b.qs, b.ts
            for n, op in b.cigar:
                if op in "=XM":
                    if op == "X":
                        mm_q.extend(range(qp, qp + n))
                        mm_t.extend(range(tp, tp + n))
                    if runs and runs[-1][0] + runs[-1][2] == qp and runs[-1][1] + runs[-1][2] == tp:
                        r = runs[-1]
                        runs[-1] = (r[0], r[1], r[2] + n)
                    else:
                        runs.append((qp, tp, n))
                    qp += n
                    tp += n
                elif op == "I":
                    gaps.append((qp, tp, n, 0))
                    qp += n
                elif op == "D":
                    gaps.append((qp, tp, 0, n))
                    tp += n
            prev = b
        return cls(
            q0=chain[0].qs,
            q1=chain[-1].qe,
            t0=chain[0].ts,
            t1=chain[-1].te,
            runs=runs,
            mm_q=mm_q,
            mm_t=mm_t,
            gaps=gaps,
        )

    def inverted(self) -> PairAlignment:
        return PairAlignment(
            q0=self.t0,
            q1=self.t1,
            t0=self.q0,
            t1=self.q1,
            runs=[(t, q, n) for q, t, n in self.runs],
            mm_q=self.mm_t,
            mm_t=self.mm_q,
            gaps=[(t, q, tl, ql) for q, t, ql, tl in self.gaps],
        )

    def map_pos(self, pos: int, right: bool = False) -> int:
        """Position on t of position pos on q (a gap maps to its left/right anchor)."""
        if not self.runs:
            return self.t0
        i = bisect_right([r[0] for r in self.runs], pos) - 1
        if i < 0:
            return self.runs[0][1]
        q0, t0, n = self.runs[i]
        if pos < q0 + n:
            return t0 + (pos - q0)
        if right and i + 1 < len(self.runs):
            return self.runs[i + 1][1]
        return t0 + n

    def window_differences(self, w0: int, w1: int) -> tuple[int, int]:
        """(indel bases, mismatches) of the alignment inside [w0, w1) on q."""
        indel = 0
        for gq, _gt, ql, tl in self.gaps:
            if ql > 0:
                inside = gq < w1 and gq + ql > w0
            else:
                inside = w0 <= gq <= w1
            if inside:
                indel += max(ql, tl)
        lo = bisect_right(self.mm_q, w0 - 1)
        hi = bisect_right(self.mm_q, w1 - 1)
        return indel, hi - lo


def project(aln: PairAlignment, start: int, end: int) -> tuple[int, int] | None:
    """[start, end) on q projected onto t, or None if it is not inside the chain."""
    if start < aln.q0 or end > aln.q1:
        return None
    ps = aln.map_pos(start)
    pe = aln.map_pos(end, right=True)
    return (ps, pe) if ps <= pe else (pe, ps)


def _overlaps(a0: int, a1: int, b0: int, b1: int, slack: int) -> bool:
    """Intervals [a0, a1] and [b0, b1] (points allowed) within slack of each other."""
    return a0 - slack <= b1 and b0 <= a1 + slack


@attrs.define
class PairDecision:
    #: the alignment covers both SVs with flanks
    evaluable: bool
    #: same allele: corresponding places, haplotypes agree across both SVs
    link: bool
    #: corresponding places, haplotypes differ there (different alleles)
    distinct: bool
    #: indel bases / SV size (lower is better)
    score: float
    detail: dict
    #: no corresponding places, but the haplotypes agree across both SVs: two
    #: placements of one event, or two events both haplotypes carry
    both_same: bool = False


def decide_pair(
    x: CoreSV,
    y: CoreSV,
    aln_xy: PairAlignment | None,
    len_x: int,
    len_y: int,
    params: CoreAlnParams,
    span_x: tuple[int, int] | None = None,
    span_y: tuple[int, int] | None = None,
) -> PairDecision:
    """Alignment-based decision for composite x (on core A = x.key) and y (on B).

    1. Correspondence: x's placement range projected onto B overlaps y's (within
       link_slack), or y's projected onto A overlaps x's.
    2. Haplotype agreement: the A-B alignment across the hull of x's range and
       y's projection on A, plus `flank` on each side, has at most
       max(abs_tol, frac_tol * larger size) indel bases (and the inserted bases
       of x at most max_mismatch_frac mismatches). The hull makes a difference
       anywhere between or around the two SVs count: in a long, imperfect VNTR
       a window around x alone can look identical while the other haplotype has
       extra units next to it.
    Evaluable only if both ranges project and the hull has min_flank aligned
    bases on both sides. One exception: when the reference spans of both cores
    (span_x, span_y) cover both SVs, the two cores hold the locus, and a range
    that does not project means the alignment breaks there (a large difference
    or an unalignable repeat): judged, not linked.
    """
    lo_ref = min(x.ref_start, y.ref_start)
    hi_ref = max(x.ref_end, y.ref_end)
    covered = (
        span_x is not None
        and span_y is not None
        and span_x[0] <= lo_ref
        and hi_ref <= span_x[1]
        and span_y[0] <= lo_ref
        and hi_ref <= span_y[1]
    )
    if aln_xy is None:
        return PairDecision(False, False, False, 1.0, {"reason": "no_alignment"})
    inv = aln_xy.inverted()
    px = project(aln_xy, x.amb_start, x.amb_end)  # x on B
    py = project(inv, y.amb_start, y.amb_end)  # y on A
    detail = {
        "x_iv": (x.start, x.end, x.amb_start, x.amb_end),
        "y_iv": (y.start, y.end, y.amb_start, y.amb_end),
        "x_proj": px,
        "y_proj": py,
    }
    if px is None or py is None:
        detail["reason"] = "not_projected" + ("_covered" if covered else "")
        return PairDecision(covered, False, False, 1.0, detail)
    corr = _overlaps(px[0], px[1], y.amb_start, y.amb_end, params.link_slack) or _overlaps(
        py[0], py[1], x.amb_start, x.amb_end, params.link_slack
    )
    lo = min(x.amb_start, py[0])
    hi = max(x.amb_end, py[1])
    if lo - aln_xy.q0 < params.min_flank or aln_xy.q1 - hi < params.min_flank:
        detail["reason"] = "no_flank"
        return PairDecision(False, False, False, 1.0, detail)
    w0 = max(aln_xy.q0, lo - params.flank)
    w1 = min(aln_xy.q1, hi + params.flank, len_x)
    indel, _ = aln_xy.window_differences(w0, w1)
    _, mm = aln_xy.window_differences(x.start, max(x.end, x.start + 1))
    size = max(x.size, y.size)
    tol = max(params.abs_tol, params.frac_tol * size)
    same = indel <= tol
    if same and x.end - x.start > 0 and mm > params.max_mismatch_frac * (x.end - x.start) + 2:
        same = False
    detail.update({"corr": corr, "indel": indel, "mm": mm, "window": (w0, w1), "same": same})
    return PairDecision(
        evaluable=True,
        link=corr and same,
        distinct=corr and not same,
        score=indel / max(1, size),
        both_same=(not corr) and same,
        detail=detail,
    )


# --------------------------------------------------------------------------- #
# merging one candidate group
# --------------------------------------------------------------------------- #


@attrs.define
class GroupStats:
    groups: int = 0
    minimap2_calls: int = 0
    minimap2_seconds: float = 0.0
    minimap2_timeouts: int = 0
    groups_too_large: int = 0
    cores_aligned: int = 0
    consensus_pairs_aligned: int = 0
    composite_pairs: int = 0
    evaluable_pairs: int = 0
    links_alignment: int = 0
    links_old_rule: int = 0
    links_both: int = 0
    links_alignment_only: int = 0
    links_old_only: int = 0
    links_fallback: int = 0
    links_both_same_old: int = 0
    blocks_off_diagonal: int = 0
    distinct_pairs: int = 0
    cannot_link_blocked: int = 0
    unplaced_composites: int = 0

    def add(self, other: GroupStats) -> None:
        for f in attrs.fields(GroupStats):
            setattr(self, f.name, getattr(self, f.name) + getattr(other, f.name))


def _core_name(key: tuple[str, str]) -> str:
    return f"{key[0]}|{key[1]}"


def _consensus_strands(group: list[SVcomposite]) -> dict[tuple[str, str], bool]:
    """Reference strand of each consensus: the majority aln_is_reverse of its composites."""
    votes: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    for svc in group:
        key = composite_consensus_key(svc)
        strand = composite_strand(svc)
        if key is not None and strand is not None:
            votes[key][int(strand)] += 1
    return {k: v[1] > v[0] for k, v in votes.items()}


def place_composites(
    group: list[SVcomposite],
    cores: dict[tuple[str, str], tuple[int, str]],
) -> tuple[list[CoreSV | None], dict[tuple[str, str], str]]:
    """Each composite's CoreSV (or None) and the oriented core of each consensus."""
    strands = _consensus_strands(group)
    oriented: dict[tuple[str, str], str] = {}
    placed: list[CoreSV | None] = []
    for svc in group:
        key = composite_consensus_key(svc)
        if key is None or key not in cores or key not in strands:
            placed.append(None)
            continue
        core_start, seq = cores[key]
        if key not in oriented:
            oriented[key] = revcomp(seq) if strands[key] else seq
        placed.append(
            core_interval_of_composite(svc, core_start, oriented[key], strands[key])
        )
    return placed, oriented


def core_reference_spans(
    composites: list[SVcomposite], placed: list[CoreSV | None]
) -> tuple[dict[tuple[str, str], tuple[int, int]], dict[tuple[str, str], int]]:
    """Reference span of each core (its consensus_aln_interval, most common
    value) and the summed SV sizes of each core's composites."""
    spans: dict[tuple[str, str], Counter] = defaultdict(Counter)
    drift: dict[tuple[str, str], int] = defaultdict(int)
    for svc, p in zip(composites, placed, strict=True):
        if p is None:
            continue
        for pattern in svc.svPatterns:
            for prim in pattern.SVprimitives:
                interval = prim.consensus_aln_interval
                if interval is not None:
                    a, b = interval[1], interval[2]
                    spans[p.key][(min(a, b), max(a, b))] += 1
        drift[p.key] += p.size
    return {k: c.most_common(1)[0][0] for k, c in spans.items()}, dict(drift)


def align_group_cores(
    oriented: dict[tuple[str, str], str],
    keys: set[tuple[str, str]],
    params: CoreAlnParams,
    stats: GroupStats,
    ref_start: dict[tuple[str, str], int] | None = None,
    drift: dict[tuple[str, str], int] | None = None,
) -> dict[tuple[tuple[str, str], tuple[str, str]], PairAlignment]:
    """Co-linear alignments between all pairs of the given cores (both orders).

    With ref_start (reference position of each oriented core's first base) a
    block is kept only near the diagonal the reference placements imply:
    |(t - q) - (ref_start[q] - ref_start[t])| <= diagonal_band + drift[q] +
    drift[t] (drift: the indel bases the core's own SVs shift it by). This drops
    blocks to another copy of a duplicated or repeated sequence in the cores.
    """
    seqs = {_core_name(k): oriented[k] for k in keys}
    if len(seqs) < 2:
        return {}
    if sum(len(s) for s in seqs.values()) > params.max_group_bp:
        stats.groups_too_large += 1
        return {}
    t0 = time.perf_counter()
    raw = run_minimap2(seqs, params)
    stats.minimap2_seconds += time.perf_counter() - t0
    stats.minimap2_calls += 1
    if raw is None:
        stats.minimap2_timeouts += 1
        return {}
    stats.cores_aligned += len(seqs)
    by_name = {_core_name(k): k for k in keys}
    out: dict[tuple[tuple[str, str], tuple[str, str]], PairAlignment] = {}
    for (qn, tn), blocks in raw.items():
        if qn == tn:
            continue
        if ref_start is not None:
            qk0, tk0 = by_name[qn], by_name[tn]
            if qk0 in ref_start and tk0 in ref_start:
                expected = ref_start[qk0] - ref_start[tk0]
                band = params.diagonal_band + (drift or {}).get(qk0, 0) + (drift or {}).get(tk0, 0)
                kept = [b for b in blocks if abs((b.ts - b.qs) - expected) <= band]
                stats.blocks_off_diagonal += len(blocks) - len(kept)
                blocks = kept
        chain = best_colinear_chain(blocks)
        if not chain:
            continue
        aln = PairAlignment.from_chain(chain)
        qk, tk = by_name[qn], by_name[tn]
        if (qk, tk) in out:
            continue
        out[(qk, tk)] = aln
        out[(tk, qk)] = aln.inverted()
        stats.consensus_pairs_aligned += 1
    return out


def merge_same_type_with_core_alignments(
    composites: list[SVcomposite],
    placed: list[CoreSV | None],
    oriented: dict[tuple[str, str], str],
    alignments: dict[tuple[tuple[str, str], tuple[str, str]], PairAlignment],
    old_rule: Callable[[SVcomposite, SVcomposite], bool],
    mode: str,
    params: CoreAlnParams,
    stats: GroupStats,
    pair_log: list[dict] | None = None,
    spans: dict[tuple[str, str], tuple[int, int]] | None = None,
) -> list[SVcomposite]:
    """Union of composites of one SV type: alignment links, old rule, constraints.

    Links are applied best first (alignment links by their divergence, then the
    old rule's); a link is skipped when it would put two composites of the same
    sample from different consensuses into one record although the alignment
    found them to be different alleles.
    """
    n = len(composites)
    if n == 0:
        return []
    spans = spans or {}
    links: list[tuple[float, int, int]] = []
    cannot: set[tuple[int, int]] = set()
    decisions: dict[tuple[int, int], PairDecision | None] = {}
    # represented[(i, key)]: composite i lands on a composite of consensus key
    # that is present there (an x_on_y with the haplotypes agreeing)
    represented: set[tuple[int, tuple[str, str]]] = set()
    for i in range(n):
        for j in range(i + 1, n):
            pa, pb = placed[i], placed[j]
            if pa is not None and pb is not None and pa.key == pb.key:
                continue  # one haplotype: two composites are two events
            dec: PairDecision | None = None
            if pa is not None and pb is not None:
                aln = alignments.get((pa.key, pb.key))
                dec = decide_pair(
                    pa,
                    pb,
                    aln,
                    len(oriented[pa.key]),
                    len(oriented[pb.key]),
                    params,
                    spans.get(pa.key),
                    spans.get(pb.key),
                )
                if dec.link:
                    represented.add((i, pb.key))
                    represented.add((j, pa.key))
            decisions[(i, j)] = dec
    for (i, j), dec in decisions.items():
        a, b = composites[i], composites[j]
        pa, pb = placed[i], placed[j]
        stats.composite_pairs += 1
        evaluable = dec is not None and dec.evaluable
        old = old_rule(a, b)
        # both SVs present on the other haplotype, neither represented there by
        # another composite: the same event placed differently, if the pairwise
        # rule agrees
        unrepresented_twin = (
            evaluable
            and dec.both_same
            and (i, pb.key) not in represented
            and (j, pa.key) not in represented
        )
        if evaluable:
            stats.evaluable_pairs += 1
            if dec.link:
                stats.links_alignment += 1
            if dec.distinct:
                stats.distinct_pairs += 1
            if dec.link and old:
                stats.links_both += 1
            elif dec.link:
                stats.links_alignment_only += 1
            elif old:
                stats.links_old_only += 1
        if old:
            stats.links_old_rule += 1
        if mode == "align":
            if evaluable:
                link = dec.link or (unrepresented_twin and old)
                if not dec.link and unrepresented_twin and old:
                    stats.links_both_same_old += 1
            else:
                link = old
                if old:
                    stats.links_fallback += 1
        else:  # "or"
            link = old or (evaluable and dec.link)
        if link:
            score = dec.score if (dec is not None and dec.link) else 1.0
            links.append((score, i, j))
        if evaluable and dec.distinct and pa.key[0] == pb.key[0]:
            cannot.add((i, j))
        if pair_log is not None:
            pair_log.append({
                "a": _pair_log_id(a),
                "b": _pair_log_id(b),
                "size_a": a.get_size(),
                "size_b": b.get_size(),
                "old": old,
                "evaluable": evaluable,
                "link_aln": bool(dec.link) if dec is not None else None,
                "distinct": bool(dec.distinct) if dec is not None else None,
                "both_same": bool(dec.both_same) if dec is not None else None,
                "twin": bool(unrepresented_twin),
                "final": link,
                **({k: v for k, v in dec.detail.items()} if dec is not None else {}),
            })
    for idx, p in enumerate(placed):
        if p is None:
            stats.unplaced_composites += 1
    uf = UnionFind(range(n))
    members: dict[int, set[int]] = {i: {i} for i in range(n)}
    for _score, i, j in sorted(links):
        ri, rj = uf.find(i), uf.find(j)
        if ri == rj:
            continue
        mi, mj = members[ri], members[rj]
        if any((min(u, v), max(u, v)) in cannot for u in mi for v in mj):
            stats.cannot_link_blocked += 1
            continue
        uf.union(ri, rj)
        r = uf.find(ri)
        merged = mi | mj
        members.pop(ri, None)
        members.pop(rj, None)
        members[r] = merged
    out: list[SVcomposite] = []
    for cc in uf.get_connected_components(allow_singletons=True):
        patterns = [p for idx in cc for p in composites[idx].svPatterns]
        out.append(SVcomposite.from_SVpatterns(svPatterns=patterns))
    return out


def _pair_log_id(svc: SVcomposite) -> str:
    keys = sorted({f"{p.samplename}:{p.consensusID}" for p in svc.svPatterns})
    chrom, pos = svc.ref_start
    return f"{svc.sv_type.get_sv_type()}@{chrom}:{pos}[{','.join(keys)}]"


def merge_group_with_core_alignments(
    group: list[SVcomposite],
    cores: dict[tuple[str, str], tuple[int, str]],
    old_rule_ins: Callable[[SVcomposite, SVcomposite], bool],
    old_rule_del: Callable[[SVcomposite, SVcomposite], bool],
    mode: str,
    params: CoreAlnParams,
    stats: GroupStats,
    pair_log: list[dict] | None = None,
) -> tuple[list[SVcomposite], list[SVcomposite]]:
    """(merged INS and DEL composites, the group's other composites untouched)."""
    stats.groups += 1
    ins = [
        svc for svc in group if issubclass(svc.sv_type, SVpatterns.SVpatternInsertion)
    ]
    dels = [
        svc for svc in group if issubclass(svc.sv_type, SVpatterns.SVpatternDeletion)
    ]
    rest = [
        svc
        for svc in group
        if not issubclass(svc.sv_type, SVpatterns.SVpatternInsertion)
        and not issubclass(svc.sv_type, SVpatterns.SVpatternDeletion)
    ]
    indels = ins + dels
    placed, oriented = place_composites(indels, cores)
    keys_ins = {p.key for p in placed[: len(ins)] if p is not None}
    keys_del = {p.key for p in placed[len(ins) :] if p is not None}
    keys = keys_ins | keys_del
    alignments = {}
    spans, drift = core_reference_spans(indels, placed)
    if len(keys_ins) > 1 or len(keys_del) > 1:
        alignments = align_group_cores(
            oriented, keys, params, stats, {k: v[0] for k, v in spans.items()}, drift
        )
    merged = merge_same_type_with_core_alignments(
        ins,
        placed[: len(ins)],
        oriented,
        alignments,
        old_rule_ins,
        mode,
        params,
        stats,
        pair_log,
        spans,
    ) + merge_same_type_with_core_alignments(
        dels,
        placed[len(ins) :],
        oriented,
        alignments,
        old_rule_del,
        mode,
        params,
        stats,
        pair_log,
        spans,
    )
    return merged, rest
