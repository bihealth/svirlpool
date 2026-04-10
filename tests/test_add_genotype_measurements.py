"""
Tests for add_genotypeMeasurements_to_SVprimitives.

Key invariant being tested:
  intervals_cutread_alignments stores core-relative coordinates [0, core_size).
  svp.read_start / svp.read_end are padded-consensus coordinates.
  The function must subtract core_interval_start before comparing.
"""

import pysam
import pytest

from svirlpool.localassembly.SVprimitives import (
    SVprimitive,
    add_genotypeMeasurements_to_SVprimitives,
)
from svirlpool.util.datatypes import MergedSVSignal


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_pysam_aln(seq_len: int, ref_start: int = 0) -> pysam.AlignedSegment:
    """Return a simple all-match pysam alignment of *seq_len* bases."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "ref", "LN": seq_len + ref_start + 1}],
    })
    seg = pysam.AlignedSegment(header)
    seg.query_name = "consensus"
    seg.query_sequence = "A" * seq_len
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = ref_start
    seg.mapping_quality = 60
    seg.cigar = [(0, seq_len)]  # seq_len M
    seg.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
    return seg


def _make_svp(
    read_start: int,
    read_end: int,
    sv_type: int,
    consensusID: str = "1.0",
) -> SVprimitive:
    """Return a minimal SVprimitive with the given padded-consensus coordinates."""
    merged = MergedSVSignal(
        ref_start=0,
        ref_end=0,
        read_start=read_start,
        read_end=read_end,
        size=abs(read_end - read_start),
        sv_type=sv_type,
        chr="chr1",
        repeatIDs=[],
        original_alt_sequences=[],
        original_ref_sequences=[],
    )
    return SVprimitive.from_merged_sv_signal(
        merged_sv_signal=merged,
        samplename="test_sample",
        consensusID=consensusID,
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=("chr1", 0, 0),
    )


# ---------------------------------------------------------------------------
# Tests: coordinate correctness (core-relative vs padded)
# ---------------------------------------------------------------------------


def test_insertion_reads_found_with_core_relative_coordinates():
    """SVprimitive at padded position 510 (core_interval_start=47009) → core pos 510.
    A read interval covering core pos 510 should match. Mirrors real log data for crID=526.
    Uses sv_type=1 (start-only path) to keep the assertion simple."""
    core_interval_start = 47009
    read_start_padded = core_interval_start + 510  # = 47519

    svp = _make_svp(read_start=read_start_padded, read_end=read_start_padded, sv_type=1)
    pysam_aln = _make_pysam_aln(seq_len=read_start_padded + 100, ref_start=0)

    # All three intervals start at or before core pos 510
    intervals = [
        (500, 2679, "read1", True),  # starts at 500 < 510
        (100, 2960, "read2", True),  # starts at 100 < 510
        (210, 1009, "read3", True),  # starts at 210 < 510
    ]

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=intervals,
        core_interval_start=core_interval_start,
    )

    assert svp.genotypeMeasurement is not None
    assert set(svp.genotypeMeasurement.supporting_reads_start) == {
        "read1",
        "read2",
        "read3",
    }
    assert (
        svp.genotypeMeasurement.supporting_reads_end is None
    )  # sv_type=1 is start-only


def test_insertion_no_reads_when_position_outside_all_intervals():
    """Core position 5000 is beyond all intervals [0, 4568] → 0 supporting reads."""
    core_interval_start = 1000
    read_start_padded = core_interval_start + 5000  # core pos 5000

    svp = _make_svp(read_start=read_start_padded, read_end=read_start_padded, sv_type=0)
    pysam_aln = _make_pysam_aln(seq_len=read_start_padded + 100, ref_start=0)

    intervals = [
        (0, 4568, "read1", True),
        (100, 2000, "read2", True),
    ]

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=intervals,
        core_interval_start=core_interval_start,
    )

    assert svp.genotypeMeasurement is not None
    assert svp.genotypeMeasurement.supporting_reads_start == []


def test_sv_type_1_uses_start_only():
    """sv_type=1 (small inline DEL) uses only the start breakpoint — end is None."""
    core_interval_start = 0
    svp = _make_svp(read_start=100, read_end=200, sv_type=1)
    pysam_aln = _make_pysam_aln(seq_len=300, ref_start=0)

    intervals = [(50, 250, "read1", True)]

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=intervals,
        core_interval_start=core_interval_start,
    )

    assert svp.genotypeMeasurement is not None
    assert svp.genotypeMeasurement.supporting_reads_start == ["read1"]
    assert svp.genotypeMeasurement.supporting_reads_end is None
    assert svp.genotypeMeasurement.end_on_consensus is None


def test_deletion_both_breakpoints_counted():
    """Large deletions are represented as BND pairs (sv_type=3 or 4), not sv_type=1.
    Both start and end breakpoints get independent supporting read counts."""
    core_interval_start = 200
    # sv spans core positions 100-300 (padded: 300-500)
    read_start_padded = core_interval_start + 100  # = 300
    read_end_padded = core_interval_start + 300  # = 500

    svp = _make_svp(
        read_start=read_start_padded, read_end=read_end_padded, sv_type=3
    )  # BND (large DEL breakpoint)
    pysam_aln = _make_pysam_aln(seq_len=600, ref_start=0)

    intervals = [
        (50, 200, "read_start_only", True),  # covers start (100) but not end (300)
        (250, 350, "read_end_only", True),  # covers end (300) but not start (100)
        (0, 400, "read_both", True),  # covers both
    ]

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=intervals,
        core_interval_start=core_interval_start,
    )

    assert svp.genotypeMeasurement is not None
    assert set(svp.genotypeMeasurement.supporting_reads_start) == {
        "read_start_only",
        "read_both",
    }
    assert set(svp.genotypeMeasurement.supporting_reads_end) == {
        "read_end_only",
        "read_both",
    }


def test_empty_intervals_gives_no_supporting_reads():
    """No intervals → no supporting reads for any SV type."""
    core_interval_start = 100
    read_start_padded = core_interval_start + 50

    svp = _make_svp(read_start=read_start_padded, read_end=read_start_padded, sv_type=0)
    pysam_aln = _make_pysam_aln(seq_len=300, ref_start=0)

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=[],
        core_interval_start=core_interval_start,
    )

    assert svp.genotypeMeasurement is not None
    assert svp.genotypeMeasurement.supporting_reads_start == []


def test_reversed_interval_tuple_is_handled():
    """Intervals stored as (end, start) (start > end) must still match correctly."""
    core_interval_start = 0
    read_start_padded = 500

    svp = _make_svp(read_start=read_start_padded, read_end=read_start_padded, sv_type=0)
    pysam_aln = _make_pysam_aln(seq_len=1000, ref_start=0)

    # Deliberately store interval reversed
    intervals = [(800, 100, "read1", False)]  # should be normalised to (100, 800)

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=intervals,
        core_interval_start=core_interval_start,
    )

    assert svp.genotypeMeasurement is not None
    assert "read1" in svp.genotypeMeasurement.supporting_reads_start


def test_multiple_svps_processed_independently():
    """Multiple SVprimitives on the same consensus get independent read lists."""
    core_interval_start = 0

    svp_a = _make_svp(read_start=100, read_end=100, sv_type=0, consensusID="5.0")
    svp_b = _make_svp(read_start=900, read_end=900, sv_type=0, consensusID="5.0")
    pysam_aln = _make_pysam_aln(seq_len=1100, ref_start=0)

    intervals = [
        (50, 200, "reads_a", True),  # only covers svp_a at 100
        (800, 1000, "reads_b", True),  # only covers svp_b at 900
        (0, 1100, "reads_both", True),
    ]

    add_genotypeMeasurements_to_SVprimitives(
        svps=[svp_a, svp_b],
        pysam_aln=pysam_aln,
        intervals_cutread_alignments=intervals,
        core_interval_start=core_interval_start,
    )

    assert set(svp_a.genotypeMeasurement.supporting_reads_start) == {
        "reads_a",
        "reads_both",
    }
    assert set(svp_b.genotypeMeasurement.supporting_reads_start) == {
        "reads_b",
        "reads_both",
    }
