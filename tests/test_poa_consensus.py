"""Tests for the in-memory POA consensus backend."""

import random

import mappy as mp
import pytest

from svirlpool.localassembly import poa_consensus as pc

ALPHABET = "ACGT"


def random_sequence(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(ALPHABET) for _ in range(length))


def corrupt(seq: str, err: float, rng: random.Random) -> str:
    """Apply deletion-dominated noise, in the shape ONT reads have."""
    out: list[str] = []
    for base in seq:
        r = rng.random()
        if r < err * 0.55:
            continue
        if r < err * 0.80:
            out.append(rng.choice(ALPHABET))
            continue
        if r < err:
            out.append(base)
            out.append(rng.choice(ALPHABET))
            continue
        out.append(base)
    return "".join(out)


def divergence(reference: str, query: str) -> float:
    """Gap-compressed divergence of query against reference, 1.0 if unalignable."""
    aligner = mp.Aligner(seq=reference, preset="map-ont")
    hit = next((h for h in aligner.map(query) if h.is_primary), None)
    if hit is None:
        return 1.0
    return hit.NM / max(1, hit.blen)


# =============================================================================
# helpers
# =============================================================================


def test_revcomp_roundtrip():
    seq = random_sequence(500, seed=1)
    assert pc._revcomp(pc._revcomp(seq)) == seq
    assert pc._revcomp("ACGTN") == "NACGT"


@pytest.mark.parametrize(
    ("length", "window", "expected"),
    [
        (100, 500, [0, 100]),
        (500, 500, [0, 500]),
        (1400, 500, [0, 500, 1000, 1400]),
        # A final window shorter than half a window is merged into its
        # predecessor rather than left as a runt.
        (1100, 500, [0, 500, 1100]),
        (1200, 500, [0, 500, 1200]),
    ],
)
def test_window_boundaries(length, window, expected):
    assert pc._window_boundaries(length, window) == expected


def test_window_boundaries_are_monotonic_and_cover_the_backbone():
    for length in (600, 999, 1000, 1001, 5000, 20615):
        bounds = pc._window_boundaries(length, 500)
        assert bounds[0] == 0
        assert bounds[-1] == length
        assert all(b < c for b, c in zip(bounds[:-1], bounds[1:], strict=True))


# =============================================================================
# CIGAR walking
# =============================================================================


def _anchors_for(reference: str, read: str, boundaries: list[int]):
    aligner = mp.Aligner(seq=reference, preset="map-ont")
    hit = next((h for h in aligner.map(read) if h.is_primary), None)
    assert hit is not None, "test read failed to align"
    return hit, pc._ref_to_query_anchors(hit, boundaries, len(read))


def test_ref_to_query_anchors_forward_strand():
    reference = random_sequence(3000, seed=2)
    read = corrupt(reference, 0.05, random.Random(3))
    boundaries = pc._window_boundaries(len(reference), 500)
    hit, anchors = _anchors_for(reference, read, boundaries)

    assert hit.strand == 1
    # Anchors must be inside the read and ordered the same way as the
    # backbone positions they came from.
    spanned = [b for b in boundaries if b in anchors]
    assert len(spanned) >= 2
    values = [anchors[b] for b in spanned]
    assert values == sorted(values)
    assert all(0 <= v <= len(read) for v in values)
    # Each anchor should sit near its backbone position: the read is the
    # backbone plus ~5% noise, so offsets stay proportional.
    for b in spanned:
        assert abs(anchors[b] - b) < 0.25 * len(reference) + 100


def test_ref_to_query_anchors_reverse_strand_uses_alignment_orientation():
    """A reverse-strand read's anchors must index its reverse complement.

    mappy reports q_st/q_en on the original read while its CIGAR walks the
    reverse complement; getting this backwards silently extracts garbage
    fragments from every flipped read.
    """
    reference = random_sequence(3000, seed=4)
    forward = corrupt(reference, 0.05, random.Random(5))
    read = pc._revcomp(forward)
    boundaries = pc._window_boundaries(len(reference), 500)
    hit, anchors = _anchors_for(reference, read, boundaries)

    assert hit.strand == -1
    query = pc._revcomp(read)
    spanned = sorted(b for b in boundaries if b in anchors)
    assert len(spanned) >= 2
    # The fragment cut between two anchors must match the backbone stretch
    # between the same two boundaries.
    lo, hi = spanned[0], spanned[-1]
    fragment = query[anchors[lo] : anchors[hi]]
    assert len(fragment) > 0
    assert divergence(reference[lo:hi], fragment) < 0.15


# =============================================================================
# orientation
# =============================================================================


def test_orient_reads_flips_reverse_complemented_reads():
    truth = random_sequence(2000, seed=6)
    rng = random.Random(7)
    reads, flipped = {}, set()
    for i in range(10):
        seq = corrupt(truth, 0.05, rng)
        if i % 2:
            seq = pc._revcomp(seq)
            flipped.add(f"r{i}")
        reads[f"r{i}"] = seq

    stats = pc.PoaAssemblyStats()
    oriented = pc.orient_reads(reads, stats=stats)

    assert set(oriented) == set(reads)
    assert stats.n_reads_flipped == len(flipped)
    # Every read now agrees with the truth in the same orientation.
    for seq in oriented.values():
        assert divergence(truth, seq) < 0.15


def test_orient_reads_passes_through_tiny_inputs():
    assert pc.orient_reads({}) == {}
    assert pc.orient_reads({"a": "ACGT"}) == {"a": "ACGT"}


def test_orient_reads_keeps_unalignable_reads():
    """An unrelated read is kept, not dropped: it may be a divergent haplotype."""
    truth = random_sequence(1500, seed=8)
    rng = random.Random(9)
    reads = {f"r{i}": corrupt(truth, 0.05, rng) for i in range(6)}
    reads["stranger"] = random_sequence(1500, seed=10)

    stats = pc.PoaAssemblyStats()
    oriented = pc.orient_reads(reads, stats=stats)

    assert "stranger" in oriented
    assert stats.n_reads_unoriented == 1


# =============================================================================
# polishing
# =============================================================================


def test_polish_removes_backbone_errors():
    truth = random_sequence(4000, seed=11)
    rng = random.Random(12)
    reads = {f"r{i}": corrupt(truth, 0.08, rng) for i in range(20)}
    backbone = reads["r0"]

    polished = pc.polish_windowed_poa(backbone, reads)

    assert divergence(truth, polished) < divergence(truth, backbone) / 2


def test_polish_uses_reads_shorter_than_a_window():
    """Reads shorter than window_size must still be able to vote.

    Requiring a read to span a whole window silences every read shorter than
    ``window_size``, which in a cut-read cluster is most of them -- the
    backbone then comes back completely unpolished.
    """
    truth = random_sequence(2000, seed=13)
    rng = random.Random(14)
    backbone = corrupt(truth, 0.08, rng)
    # 300 bp reads against a 500 bp window: none of them spans a window.
    reads = {}
    for i in range(30):
        start = (i * 60) % (len(truth) - 300)
        reads[f"r{i}"] = corrupt(truth[start : start + 300], 0.05, rng)

    polished = pc.polish_windowed_poa(backbone, reads, window_size=500)

    assert polished != backbone
    assert divergence(truth, polished) < divergence(truth, backbone)


def test_polish_returns_backbone_without_reads():
    backbone = random_sequence(800, seed=15)
    assert pc.polish_windowed_poa(backbone, {}) == backbone
    assert pc.polish_windowed_poa("", {"r": "ACGT"}) == ""


# =============================================================================
# layout
# =============================================================================


def test_layout_score_prefers_the_consensus_that_hosts_more_reads():
    truth = random_sequence(4000, seed=16)
    rng = random.Random(17)
    reads = {f"r{i}": corrupt(truth, 0.05, rng) for i in range(8)}

    full = pc.layout_score(truth, reads)
    half = pc.layout_score(truth[:2000], reads)

    assert full > half
    assert pc.layout_score("", reads) == 0.0


def test_extend_backbone_grows_along_overhanging_reads():
    truth = random_sequence(6000, seed=18)
    rng = random.Random(19)
    # A seed read covering the middle, plus reads overhanging both ends.
    reads = {
        "seed": corrupt(truth[1500:4500], 0.05, rng),
        "left": corrupt(truth[0:3000], 0.05, rng),
        "right": corrupt(truth[3000:6000], 0.05, rng),
    }
    backbone = reads["seed"]

    extended = pc.extend_backbone(backbone, reads)

    assert len(extended) > len(backbone)
    assert pc.layout_score(extended, reads) > pc.layout_score(backbone, reads)


def test_extend_backbone_refuses_a_join_that_explains_nothing():
    """An unrelated read must not be welded onto the backbone."""
    truth = random_sequence(3000, seed=20)
    rng = random.Random(21)
    reads = {f"r{i}": corrupt(truth, 0.05, rng) for i in range(5)}
    reads["stranger"] = random_sequence(3000, seed=22)

    backbone = reads["r0"]
    extended = pc.extend_backbone(backbone, reads)

    # Whatever it does, it must not have spliced the unrelated sequence on.
    assert len(extended) < len(backbone) + 1000


# =============================================================================
# end to end
# =============================================================================


def test_assemble_consensus_poa_on_stacked_reads():
    truth = random_sequence(4000, seed=23)
    rng = random.Random(24)
    reads = {f"r{i}": corrupt(truth, 0.08, rng) for i in range(20)}

    stats = pc.PoaAssemblyStats()
    consensus = pc.assemble_consensus_poa(reads, name="stacked", stats=stats)

    assert consensus is not None
    assert divergence(truth, consensus) < 0.01
    assert abs(len(consensus) - len(truth)) < 0.05 * len(truth)
    assert stats.n_reads == 20
    assert stats.polish_rounds == 1


def test_assemble_consensus_poa_does_not_collapse_a_tiled_cluster():
    """A depth peak must not swallow the rest of the locus.

    This is the muc1 failure mode: a handful of long reads span an insertion
    while many short reads pile up on one flank.  abPOA's consensus is the
    heaviest path through the graph, so on its own it follows the pile and
    returns a consensus the length of the short reads, losing the insertion
    entirely.
    """
    truth = random_sequence(9000, seed=25)
    rng = random.Random(26)
    reads = {f"long{i}": corrupt(truth, 0.06, rng) for i in range(3)}
    for i in range(20):
        reads[f"short{i}"] = corrupt(truth[8400:9000], 0.06, rng)

    consensus = pc.assemble_consensus_poa(reads, name="tiled")

    assert consensus is not None
    assert len(consensus) > 0.8 * len(truth), (
        f"locus collapsed to {len(consensus)} bp of {len(truth)}"
    )
    assert divergence(truth, consensus) < 0.05


def test_assemble_consensus_poa_edge_cases():
    assert pc.assemble_consensus_poa({}) is None
    assert pc.assemble_consensus_poa({"a": ""}) is None
    only = random_sequence(300, seed=27)
    assert pc.assemble_consensus_poa({"a": only}) == only


def test_assemble_consensus_poa_uppercases_input():
    truth = random_sequence(1200, seed=28)
    rng = random.Random(29)
    reads = {f"r{i}": corrupt(truth, 0.05, rng).lower() for i in range(8)}

    consensus = pc.assemble_consensus_poa(reads)

    assert consensus is not None
    assert consensus == consensus.upper()


def test_assemble_consensus_poa_handles_mixed_strands():
    truth = random_sequence(3000, seed=30)
    rng = random.Random(31)
    reads = {}
    for i in range(12):
        seq = corrupt(truth, 0.06, rng)
        reads[f"r{i}"] = pc._revcomp(seq) if i % 2 else seq

    consensus = pc.assemble_consensus_poa(reads)

    assert consensus is not None
    assert divergence(truth, consensus) < 0.02
