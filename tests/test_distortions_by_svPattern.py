"""Tests for ``SVpatterns.distortions_by_svPattern`` and ``util.exponential_weight``.

Coordinate-space invariant under test
------------------------------------
``ConsensusDistortion.position`` is ``signal.ref_start`` of a *cut read aligned to the
core consensus FASTA* (``consensus.final_consensus``).  Cut reads are aligned before
``create_padding_for_consensus`` runs, so those positions -- and the
``intervals_cutread_alignments`` produced from the very same alignments -- live in
**core-consensus** coordinates, i.e. ``[0, len(consensus.consensus_sequence)]``.

``SVprimitive.read_start`` / ``read_end`` come from the *padded* consensus-to-reference
alignment, so they live in **padded-consensus** coordinates.  The conversion between
the two is ``core = padded - consensus_padding.padding_size_left`` -- exactly what
``SVprimitives.add_genotypeMeasurements_to_SVprimitives`` does (it is handed
``core_interval_start=consensus.consensus_padding.padding_size_left``).

``SVprimitive.ref_start`` / ``ref_end`` are **reference-genome** coordinates and must
never be compared against a distortion position.
"""

import math

import numpy as np
import pytest

from svirlpool.localassembly.consensus_class import Consensus, ConsensusPadding
from svirlpool.localassembly.SVpatterns import (
    SVpatternDeletion,
    SVpatternInsertion,
    distortions_by_svPattern,
)
from svirlpool.localassembly.SVprimitives import SVprimitive
from svirlpool.svcalling.genotyping import GenotypeMeasurement
from svirlpool.util.datatypes import MergedSVSignal, ReadAlignmentSignals, SVsignal
from svirlpool.util.util import exponential_weight

# Production defaults, see consensus_align.py --distance-scale / --falloff
DISTANCE_SCALE = 5000.0
FALLOFF = 1.0


# ---------------------------------------------------------------------------
# Builders -- these construct the same objects the production code constructs
# ---------------------------------------------------------------------------


def _make_consensus(
    core_length: int,
    padding_size_left: int,
    padding_size_right: int,
    distortions: list[tuple[str, int, int, int]],
    consensusID: str = "7.0",
    original_regions: list[tuple[str, int, int]] | None = None,
) -> Consensus:
    """Build a ``Consensus`` carrying *distortions* as cut-read alignment signals.

    ``distortions`` is a list of ``(readname, core_position, size, sv_type)`` tuples,
    where ``core_position`` is a **core-consensus** offset -- the space
    ``final_consensus`` records ``signal.ref_start`` in.
    """
    if original_regions is None:
        original_regions = [("chr1", 1_000, 1_500)]

    signals_by_read: dict[str, list[SVsignal]] = {}
    for readname, core_position, size, sv_type in distortions:
        signals_by_read.setdefault(readname, []).append(
            SVsignal(
                ref_start=core_position,
                ref_end=core_position + (size if sv_type == 0 else 0),
                read_start=core_position,
                read_end=core_position + size,
                size=size,
                sv_type=sv_type,
            )
        )

    cut_read_alignment_signals = [
        ReadAlignmentSignals(
            samplename="testsample",
            read_name=readname,
            reference_name=consensusID,  # cut reads align to the consensus, not the genome
            alignment_forward=True,
            SV_signals=signals,
        )
        for readname, signals in signals_by_read.items()
    ]

    # cut reads span the whole core consensus -- core-relative, as in final_consensus
    intervals_cutread_alignments = [
        (0, core_length, readname, True) for readname in signals_by_read
    ]

    padded_length = padding_size_left + core_length + padding_size_right
    padding = ConsensusPadding(
        sequence="n" * padding_size_left + "A" * core_length + "n" * padding_size_right,
        readname_left="pad_read_left",
        readname_right="pad_read_right",
        padding_size_left=padding_size_left,
        padding_size_right=padding_size_right,
        consensus_interval_on_sequence_with_padding=(
            padding_size_left,
            padded_length - padding_size_right,
        ),
    )

    return Consensus(
        ID=consensusID,
        crIDs=[int(consensusID.split(".")[0])],
        original_regions=original_regions,
        consensus_sequence="A" * core_length,
        consensus_padding=padding,
        intervals_cutread_alignments=intervals_cutread_alignments,
        cut_read_alignment_signals=cut_read_alignment_signals,
        clustering_meta_data={},
    )


def _make_svprimitive(
    *,
    chrom: str,
    ref_start: int,
    ref_end: int,
    core_start: int,
    core_end: int,
    padding_size_left: int,
    sv_type: int,
    supporting_reads: list[str],
    supporting_reads_end: list[str] | None = None,
    consensusID: str = "7.0",
    svID: int = 0,
) -> SVprimitive:
    """Build an ``SVprimitive`` whose ``read_*`` fields are padded-consensus coords.

    ``core_start`` / ``core_end`` are given in core-consensus space for readability;
    they are shifted by ``padding_size_left`` on the way in, exactly as the padded
    consensus-to-reference alignment would have produced them.
    """
    merged = MergedSVSignal(
        chr=chrom,
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=padding_size_left + core_start,
        read_end=padding_size_left + core_end,
        size=abs(core_end - core_start),
        sv_type=sv_type,
        repeatIDs=[],
        original_alt_sequences=[],
        original_ref_sequences=[],
    )
    svp = SVprimitive.from_merged_sv_signal(
        merged_sv_signal=merged,
        samplename="testsample",
        consensusID=consensusID,
        alignmentID=0,
        svID=svID,
        aln_is_reverse=False,
        consensus_aln_interval=(chrom, ref_start, ref_end),
    )
    svp.genotypeMeasurement = GenotypeMeasurement(
        start_on_consensus=ref_start,
        supporting_reads_start=list(supporting_reads),
        end_on_consensus=ref_end if supporting_reads_end is not None else None,
        supporting_reads_end=(
            list(supporting_reads_end) if supporting_reads_end is not None else None
        ),
    )
    return svp


def _insertion_pattern(
    *,
    chrom: str = "chr1",
    ref_start: int,
    core_start: int,
    core_end: int,
    padding_size_left: int,
    supporting_reads: list[str],
    consensusID: str = "7.0",
) -> SVpatternInsertion:
    return SVpatternInsertion(
        SVprimitives=[
            _make_svprimitive(
                chrom=chrom,
                ref_start=ref_start,
                ref_end=ref_start,  # insertions are a point on the reference
                core_start=core_start,
                core_end=core_end,
                padding_size_left=padding_size_left,
                sv_type=1,
                supporting_reads=supporting_reads,
                consensusID=consensusID,
            )
        ]
    )


def _weighted_mean(pairs: list[tuple[float, float]]) -> float:
    total_w = sum(w for _, w in pairs)
    return sum(s * w for s, w in pairs) / total_w


# ---------------------------------------------------------------------------
# exponential_weight -- pins the numeric claims of the defect catalogue (F1)
# ---------------------------------------------------------------------------


def test_exponential_weight_is_one_at_distance_zero():
    assert exponential_weight(distance=0, scale=DISTANCE_SCALE, falloff=FALLOFF) == 1.0


def test_exponential_weight_production_defaults_locality_table():
    """The intended locality: 0.905 at 500 bp, 0.82 at 1 kb."""
    w500 = exponential_weight(distance=500, scale=DISTANCE_SCALE, falloff=FALLOFF)
    w1000 = exponential_weight(distance=1000, scale=DISTANCE_SCALE, falloff=FALLOFF)
    assert w500 == pytest.approx(0.9048374, abs=1e-6)
    assert w1000 == pytest.approx(0.8187308, abs=1e-6)
    assert round(float(w500), 3) == 0.905
    assert round(float(w1000), 2) == 0.82


def test_exponential_weight_is_monotonically_decreasing():
    distances = [0, 100, 500, 1_000, 5_000, 20_000, 100_000]
    weights = [
        float(exponential_weight(distance=d, scale=DISTANCE_SCALE, falloff=FALLOFF))
        for d in distances
    ]
    assert weights == sorted(weights, reverse=True)


def test_exponential_weight_underflows_to_exact_zero_beyond_3_73_mb():
    """F1's numeric claim: exp(-d/5000) is 0.0 in double precision past ~3.73 Mb."""
    w_3_5mb = exponential_weight(
        distance=3_500_000, scale=DISTANCE_SCALE, falloff=FALLOFF
    )
    w_3_70mb = exponential_weight(
        distance=3_700_000, scale=DISTANCE_SCALE, falloff=FALLOFF
    )
    w_3_73mb = exponential_weight(
        distance=3_730_000, scale=DISTANCE_SCALE, falloff=FALLOFF
    )

    assert w_3_5mb == pytest.approx(9.86e-305, rel=1e-2)
    assert w_3_5mb > 0.0
    # subnormal but still non-zero
    assert 0.0 < w_3_70mb < np.finfo(np.float64).tiny
    assert w_3_73mb == 0.0

    # Locate the exact underflow boundary by bisection and pin it to ~3.72 Mb.
    lo, hi = 3_000_000, 4_000_000
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if exponential_weight(distance=mid, scale=DISTANCE_SCALE, falloff=FALLOFF) > 0:
            lo = mid
        else:
            hi = mid
    assert 3_720_000 < hi < 3_730_000
    assert exponential_weight(distance=lo, scale=DISTANCE_SCALE, falloff=FALLOFF) > 0.0
    assert exponential_weight(distance=hi, scale=DISTANCE_SCALE, falloff=FALLOFF) == 0.0
    # ... and the analytic estimate agrees: exp underflows past d/scale ~ 745.13
    assert hi == pytest.approx(DISTANCE_SCALE * -math.log(5e-324), rel=1e-3)


def test_exponential_weight_underflows_for_every_realistic_genomic_coordinate():
    """A genomic coordinate on any human chromosome arm past ~3.73 Mb kills the weight."""
    for genomic_coordinate in (10_000_000, 157_313_244, 248_000_000):
        assert (
            exponential_weight(
                distance=genomic_coordinate, scale=DISTANCE_SCALE, falloff=FALLOFF
            )
            == 0.0
        )


# ---------------------------------------------------------------------------
# Characterisation of the pre-fix behaviour (F1) -- documents the regression
# ---------------------------------------------------------------------------


def test_CHARACTERISATION_prefix_reference_coordinate_distance_underflows():
    """CHARACTERISATION TEST -- documents the *broken* pre-fix arithmetic of F1.

    Before the fix ``distortions_by_svPattern`` computed

        distance = min(|distortion.position - svp.ref_start|,
                       |distortion.position - svp.ref_end|)

    with ``distortion.position`` in core-consensus space (order 1e3-1e4) and
    ``ref_start``/``ref_end`` in reference-genome space (order 1e7-1e8).  The distance
    was therefore essentially the genomic coordinate, every weight underflowed to
    exactly 0.0, ``total_weight > 0`` failed and every read fell back to 0.0.

    This test reproduces that expression against the same fixture the fixed function
    is tested with, so the historical behaviour stays pinned even after the fix.
    """
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=41_279,
        padding_size_right=40_044,
        distortions=[
            ("read_a", 5_000, 30, 0),
            ("read_a", 5_600, 24, 1),
            ("read_b", 5_100, 28, 0),
        ],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_400,
        padding_size_left=41_279,
        supporting_reads=["read_a", "read_b"],
    )
    sv_start = pattern.SVprimitives[0].ref_start
    sv_end = pattern.SVprimitives[-1].ref_end

    prefix_weights = []
    for distortion in consensus.get_consensus_distortions():
        distance = min(
            abs(distortion.position - sv_start), abs(distortion.position - sv_end)
        )
        # the distance is, to within a few kb, the genomic coordinate itself
        assert distance > 149_000_000
        prefix_weights.append(
            float(
                exponential_weight(
                    distance=distance, scale=DISTANCE_SCALE, falloff=FALLOFF
                )
            )
        )

    assert prefix_weights == [0.0, 0.0, 0.0]
    # ... hence total_weight == 0 and the function fell back to 0.0 for every read.


def test_CHARACTERISATION_prefix_below_the_underflow_threshold_was_meaningless():
    """CHARACTERISATION TEST -- documents the *broken* pre-fix arithmetic of F1 for a
    locus **below** the ~3.72 Mb underflow threshold.

    There the pre-fix weights did not underflow, so the arm looked alive.  It was still
    meaningless: because every distortion position is consensus-local while sv_start /
    sv_end are genomic, the distance reduced to ``sv_start - position``, a monotone
    gradient along the consensus that ignores where the variant actually sits.  A
    distortion 500 bp from the variant and one 20 kb from it on the *other* side got
    weights ordered purely by consensus offset, and both were far below the intended
    0.905 / 0.018.

    This test recomputes that expression and stays green after the fix.
    """
    pad_left = 12_345
    core_start, core_end = 10_000, 10_100
    genomic_start = 730_000  # chr2:0.73 Mb, well below the underflow threshold
    consensus = _make_consensus(
        core_length=40_000,
        padding_size_left=pad_left,
        padding_size_right=6_789,
        distortions=[
            ("read_a", core_start - 500, 10, 0),  # 500 bp *before* the variant
            ("read_a", core_end + 20_000, 200, 0),  # 20 kb *after* the variant
        ],
    )
    pattern = _insertion_pattern(
        ref_start=genomic_start,
        core_start=core_start,
        core_end=core_end,
        padding_size_left=pad_left,
        supporting_reads=["read_a"],
    )
    sv_start = pattern.SVprimitives[0].ref_start
    sv_end = pattern.SVprimitives[-1].ref_end

    prefix = []
    for distortion in consensus.get_consensus_distortions():
        distance = min(
            abs(distortion.position - sv_start), abs(distortion.position - sv_end)
        )
        prefix.append((
            distortion.position,
            distance,
            float(
                exponential_weight(
                    distance=distance, scale=DISTANCE_SCALE, falloff=FALLOFF
                )
            ),
        ))

    near_pos, near_distance, near_weight = prefix[0]
    far_pos, far_distance, far_weight = prefix[1]

    # the weights do not underflow here ...
    assert near_weight > 0.0 and far_weight > 0.0
    # ... but the "distance" is just the genomic coordinate minus the consensus offset
    assert near_distance == genomic_start - near_pos
    assert far_distance == genomic_start - far_pos
    # ... so the distortion 20 kb away from the variant is ranked *closer* than the one
    # 500 bp away, purely because it sits earlier in the consensus
    assert far_weight > near_weight
    # ... and neither is anywhere near the intended locality of 0.905 / 0.018
    assert near_weight < 1e-60
    assert far_weight < 1e-60

    # The fix restores the intended ordering and magnitudes at the very same locus.
    fixed = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    w_near = float(
        exponential_weight(distance=500, scale=DISTANCE_SCALE, falloff=FALLOFF)
    )
    w_far = float(
        exponential_weight(distance=20_000, scale=DISTANCE_SCALE, falloff=FALLOFF)
    )
    assert fixed["read_a"] == pytest.approx(
        _weighted_mean([(10.0, w_near), (200.0, w_far)])
    )


# ---------------------------------------------------------------------------
# The fix: distortions must be non-zero at a large genomic coordinate
# ---------------------------------------------------------------------------


def test_distortions_are_non_zero_at_a_large_genomic_coordinate():
    """F1 verification: a locus at chr1:150,000,000 must still yield real estimates."""
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=41_279,
        padding_size_right=40_044,
        distortions=[
            ("read_a", 5_000, 30, 0),
            ("read_a", 5_600, 24, 1),
            ("read_b", 5_100, 28, 0),
            ("read_c", 18_000, 90, 0),
        ],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_400,
        padding_size_left=41_279,
        supporting_reads=["read_a", "read_b", "read_c"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )

    assert set(result) == {"read_a", "read_b", "read_c"}
    assert all(value != 0.0 for value in result.values()), result
    # read_a carries a 30 bp and a 24 bp distortion -> weighted mean between them
    assert 24.0 < result["read_a"] < 30.0
    assert result["read_b"] == pytest.approx(28.0)
    assert result["read_c"] == pytest.approx(90.0)


def test_distortion_values_depend_on_consensus_offset_not_genomic_coordinate():
    """Moving a distortion *within the consensus* changes the answer; moving the
    locus along the genome does not."""
    kwargs = {
        "core_length": 20_000,
        "padding_size_left": 1_000,
        "padding_size_right": 1_000,
    }

    near = _make_consensus(
        distortions=[("read_a", 5_100, 10, 0), ("read_a", 15_000, 100, 0)], **kwargs
    )
    far = _make_consensus(
        distortions=[("read_a", 9_000, 10, 0), ("read_a", 15_000, 100, 0)], **kwargs
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_050,
        padding_size_left=1_000,
        supporting_reads=["read_a"],
    )

    near_result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=near,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    far_result = distortions_by_svPattern(
        svPattern=pattern, consensus=far, distance_scale=DISTANCE_SCALE, falloff=FALLOFF
    )

    assert near_result["read_a"] != far_result["read_a"]
    # the 10 bp distortion sits closer to the variant in `near`, so it pulls the
    # weighted mean further down there
    assert near_result["read_a"] < far_result["read_a"]


def test_translation_invariance_along_the_genome():
    """The identical consensus-local configuration at chr1:1,000 and chr1:150,000,000
    must produce byte-identical weighted means."""
    distortions = [
        ("read_a", 5_000, 30, 0),
        ("read_a", 5_600, 24, 1),
        ("read_b", 5_100, 28, 0),
        ("read_b", 12_000, 61, 1),
    ]
    consensus_low = _make_consensus(
        core_length=20_000,
        padding_size_left=41_279,
        padding_size_right=40_044,
        distortions=distortions,
        original_regions=[("chr1", 1_000, 1_500)],
    )
    consensus_high = _make_consensus(
        core_length=20_000,
        padding_size_left=41_279,
        padding_size_right=40_044,
        distortions=distortions,
        original_regions=[("chr1", 150_000_000, 150_000_500)],
    )
    pattern_low = _insertion_pattern(
        ref_start=1_000,
        core_start=5_000,
        core_end=5_400,
        padding_size_left=41_279,
        supporting_reads=["read_a", "read_b"],
    )
    pattern_high = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_400,
        padding_size_left=41_279,
        supporting_reads=["read_a", "read_b"],
    )

    low = distortions_by_svPattern(
        svPattern=pattern_low,
        consensus=consensus_low,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    high = distortions_by_svPattern(
        svPattern=pattern_high,
        consensus=consensus_high,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )

    assert low == high
    assert all(v != 0.0 for v in low.values())


def test_padding_size_left_is_the_offset_that_makes_the_result_invariant():
    """The only offset that makes the answer independent of how much padding a
    consensus happened to acquire is ``consensus_padding.padding_size_left``.

    Two consensuses with the same core-local layout but wildly different padding must
    agree.  A fix that used 0, or the padded length, or the right padding would not.
    """
    distortions = [("read_a", 5_000, 30, 0), ("read_a", 9_000, 90, 0)]
    results = []
    for pad_left, pad_right in ((0, 0), (1_000, 25_000), (100_000, 37_033)):
        consensus = _make_consensus(
            core_length=20_000,
            padding_size_left=pad_left,
            padding_size_right=pad_right,
            distortions=distortions,
        )
        pattern = _insertion_pattern(
            ref_start=150_000_000,
            core_start=5_000,
            core_end=5_400,
            padding_size_left=pad_left,
            supporting_reads=["read_a"],
        )
        results.append(
            distortions_by_svPattern(
                svPattern=pattern,
                consensus=consensus,
                distance_scale=DISTANCE_SCALE,
                falloff=FALLOFF,
            )["read_a"]
        )

    assert results[0] == results[1] == results[2]
    assert results[0] != 0.0


def test_a_near_distortion_dominates_a_far_one_with_production_defaults():
    """The locality the noise model was designed around: with scale=5000, falloff=1.0
    a distortion 500 bp from the variant carries weight 0.905 and one 20 kb away
    carries weight 0.018, so the near one dominates the weighted mean."""
    consensus = _make_consensus(
        core_length=40_000,
        padding_size_left=12_345,
        padding_size_right=6_789,
        distortions=[
            ("read_a", 5_500, 10, 0),  # 500 bp past the pattern end
            ("read_a", 25_000, 200, 0),  # 20 kb past the pattern end
        ],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=4_900,
        core_end=5_000,
        padding_size_left=12_345,
        supporting_reads=["read_a"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )

    w_near = float(
        exponential_weight(distance=500, scale=DISTANCE_SCALE, falloff=FALLOFF)
    )
    w_far = float(
        exponential_weight(distance=20_000, scale=DISTANCE_SCALE, falloff=FALLOFF)
    )
    assert round(w_near, 3) == 0.905
    assert w_far == pytest.approx(0.0183156, abs=1e-6)

    expected = _weighted_mean([(10.0, w_near), (200.0, w_far)])
    assert result["read_a"] == pytest.approx(expected)
    # the near 10 bp distortion dominates the far 200 bp one
    assert result["read_a"] < 15.0
    # ... and it is emphatically not the unweighted mean of 105
    assert abs(result["read_a"] - 105.0) > 80.0


def test_distortion_exactly_at_the_pattern_start_gets_distance_zero():
    """A distortion at the pattern's consensus-local start must get weight exactly 1."""
    core_start, core_end = 5_000, 5_400
    pad_left = 41_279
    far_position = 15_000
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=pad_left,
        padding_size_right=40_044,
        distortions=[
            ("read_a", core_start, 12, 0),  # exactly at the pattern start
            ("read_a", far_position, 300, 0),
        ],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=core_start,
        core_end=core_end,
        padding_size_left=pad_left,
        supporting_reads=["read_a"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )

    w_far = float(
        exponential_weight(
            distance=far_position - core_end, scale=DISTANCE_SCALE, falloff=FALLOFF
        )
    )
    expected = _weighted_mean([(12.0, 1.0), (300.0, w_far)])
    assert result["read_a"] == pytest.approx(expected)


def test_distortion_exactly_at_the_pattern_end_gets_distance_zero():
    core_start, core_end = 5_000, 5_400
    pad_left = 41_279
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=pad_left,
        padding_size_right=40_044,
        distortions=[("read_a", core_end, 12, 0), ("read_a", 15_000, 300, 0)],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=core_start,
        core_end=core_end,
        padding_size_left=pad_left,
        supporting_reads=["read_a"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    w_far = float(
        exponential_weight(
            distance=15_000 - core_end, scale=DISTANCE_SCALE, falloff=FALLOFF
        )
    )
    expected = _weighted_mean([(12.0, 1.0), (300.0, w_far)])
    assert result["read_a"] == pytest.approx(expected)


def test_deletion_pattern_uses_both_primitive_boundaries():
    """For a multi-primitive pattern the distance is measured to the *nearer* of
    ``SVprimitives[0].read_start`` and ``SVprimitives[-1].read_end``."""
    pad_left = 3_000
    consensus = _make_consensus(
        core_length=30_000,
        padding_size_left=pad_left,
        padding_size_right=1_000,
        distortions=[("read_a", 20_400, 40, 1)],  # 400 bp past the pattern *end*
    )
    first = _make_svprimitive(
        chrom="chr1",
        ref_start=150_000_000,
        ref_end=150_000_500,
        core_start=5_000,
        core_end=5_010,
        padding_size_left=pad_left,
        sv_type=2,
        supporting_reads=["read_a"],
        supporting_reads_end=["read_a"],
        svID=0,
    )
    last = _make_svprimitive(
        chrom="chr1",
        ref_start=150_010_000,
        ref_end=150_010_500,
        core_start=19_990,
        core_end=20_000,
        padding_size_left=pad_left,
        sv_type=2,
        supporting_reads=["read_a"],
        supporting_reads_end=["read_a"],
        svID=1,
    )
    pattern = SVpatternDeletion(SVprimitives=[first, last])

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    # single distortion -> weighted mean is the size itself, but it must be reached
    # without any underflow, i.e. the weight has to be > 0
    assert result["read_a"] == pytest.approx(40.0)


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------


def test_no_distortions_returns_zero_for_every_supporting_read():
    consensus = _make_consensus(
        core_length=500,
        padding_size_left=28_639,
        padding_size_right=30_639,
        distortions=[],
    )
    pattern = _insertion_pattern(
        ref_start=157_299_125,
        core_start=250,
        core_end=280,
        padding_size_left=28_639,
        supporting_reads=["read_a", "read_b", "read_c"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    assert result == {"read_a": 0.0, "read_b": 0.0, "read_c": 0.0}


def test_distortions_on_non_supporting_reads_are_ignored():
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=1_000,
        padding_size_right=1_000,
        distortions=[
            ("read_a", 5_100, 30, 0),
            ("read_other", 5_100, 999, 0),  # not a supporting read
        ],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_050,
        padding_size_left=1_000,
        supporting_reads=["read_a", "read_b"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    assert set(result) == {"read_a", "read_b"}
    assert result["read_a"] == pytest.approx(30.0)
    # read_b supports the pattern but carries no distortion of its own
    assert result["read_b"] == 0.0


def test_only_non_supporting_reads_carry_distortions():
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=1_000,
        padding_size_right=1_000,
        distortions=[("read_other", 5_100, 999, 0)],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_050,
        padding_size_left=1_000,
        supporting_reads=["read_a"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    assert result == {"read_a": 0.0}


def test_single_supporting_read_with_a_single_distortion():
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=100_000,
        padding_size_right=37_033,
        distortions=[("read_a", 5_100, 42, 1)],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_050,
        padding_size_left=100_000,
        supporting_reads=["read_a"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    assert result == {"read_a": pytest.approx(42.0)}


def test_breakend_signals_are_not_counted_as_distortions():
    """``get_consensus_distortions`` only surfaces sv_type 0 (INS) and 1 (DEL)."""
    consensus = _make_consensus(
        core_length=20_000,
        padding_size_left=1_000,
        padding_size_right=1_000,
        distortions=[("read_a", 5_100, 30, 0), ("read_a", 5_200, 5_000, 3)],
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=5_000,
        core_end=5_050,
        padding_size_left=1_000,
        supporting_reads=["read_a"],
    )

    result = distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=DISTANCE_SCALE,
        falloff=FALLOFF,
    )
    assert result["read_a"] == pytest.approx(30.0)


# ---------------------------------------------------------------------------
# The coordinate-space guard
# ---------------------------------------------------------------------------


def test_distance_never_exceeds_the_consensus_length(caplog):
    consensus = _make_consensus(
        core_length=28_163,
        padding_size_left=39_224,
        padding_size_right=26_643,
        distortions=[
            ("read_a", 588, 30, 0),
            ("read_a", 27_656, 24, 1),
            ("read_b", 14_000, 28, 0),
        ],
    )
    pattern = _insertion_pattern(
        ref_start=157_311_138,
        core_start=847,
        core_end=1_733,
        padding_size_left=39_224,
        supporting_reads=["read_a", "read_b"],
    )

    with caplog.at_level("WARNING"):
        result = distortions_by_svPattern(
            svPattern=pattern,
            consensus=consensus,
            distance_scale=DISTANCE_SCALE,
            falloff=FALLOFF,
        )

    assert all(value != 0.0 for value in result.values())
    assert "coordinate space" not in caplog.text.lower()


def test_out_of_range_distortion_warns_but_does_not_raise(caplog):
    """A distortion further from the pattern than the consensus is long can only mean
    the two operands have drifted apart in coordinate space again -- warn loudly, but
    do not abort a whole-genome run."""
    consensus = _make_consensus(
        core_length=1_000,
        padding_size_left=2_000,
        padding_size_right=2_000,
        distortions=[("read_a", 50_000, 30, 0)],  # impossible: far past the core
    )
    pattern = _insertion_pattern(
        ref_start=150_000_000,
        core_start=100,
        core_end=150,
        padding_size_left=2_000,
        supporting_reads=["read_a"],
    )

    with caplog.at_level("WARNING"):
        result = distortions_by_svPattern(
            svPattern=pattern,
            consensus=consensus,
            distance_scale=DISTANCE_SCALE,
            falloff=FALLOFF,
        )

    assert set(result) == {"read_a"}
    assert "coordinate space" in caplog.text.lower()
    assert consensus.ID in caplog.text
