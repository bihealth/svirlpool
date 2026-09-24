"""--repeat-collapse-mode: which indels of one consensus are merged horizontally.

All patterns below share repeatID 7 and consensusID 1.0, so every pair is a
candidate for the union; the mode alone decides which pairs are taken.
"""

import pytest

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.multisample_sv_calling import (
    REPEAT_COLLAPSE_MODES,
    svPatterns_to_horizontally_merged_svComposites,
)

SV_TYPES = {SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion}


def _primitive(*, sv_type, ref_start, ref_end, read_start, read_end, size):
    return SVprimitives.SVprimitive(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        size=size,
        sv_type=sv_type,
        chr="chr1",
        repeatIDs=[7],
        original_alt_sequences=["A" * size] if sv_type == 0 else [],
        original_ref_sequences=["A" * size] if sv_type == 1 else [],
        samplename="S1",
        consensusID="1.0",
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=("chr1", 0, 10_000),
        genotypeMeasurement=genotyping.GenotypeMeasurement(
            start_on_consensus=read_start,
            supporting_reads_start=["r1", "r2", "r3"],
        ),
    )


def _ins(read_start, size):
    """An insertion occupying [read_start, read_start + size) on the consensus."""
    pattern = SVpatterns.SVpatternInsertion(
        SVprimitives=[
            _primitive(
                sv_type=0,
                ref_start=1000 + read_start,
                ref_end=1001 + read_start,
                read_start=read_start,
                read_end=read_start + size,
                size=size,
            )
        ]
    )
    pattern.set_sequence("A" * size)
    return pattern


def _del(read_start, size):
    """A deletion at a single consensus position."""
    pattern = SVpatterns.SVpatternDeletion(
        SVprimitives=[
            _primitive(
                sv_type=1,
                ref_start=1000 + read_start,
                ref_end=1000 + read_start + size,
                read_start=read_start,
                read_end=read_start,
                size=size,
            )
        ]
    )
    pattern.set_sequence("A" * size)
    return pattern


def _component_sizes(patterns, mode, max_gap=50):
    composites = svPatterns_to_horizontally_merged_svComposites(
        patterns, sv_types=SV_TYPES, collapse_repeats=mode, collapse_max_gap=max_gap
    )
    return sorted(len(c.svPatterns) for c in composites)


def test_modes_are_exported():
    assert REPEAT_COLLAPSE_MODES == ("all", "same-type", "adjacent", "proximal", "none")


def test_all_unions_every_indel_sharing_a_repeat():
    patterns = [_ins(0, 60), _del(400, 80), _ins(900, 70)]
    assert _component_sizes(patterns, "all") == [3]
    # the legacy boolean spelling is the same mode
    assert _component_sizes(patterns, True) == [3]


def test_none_unions_nothing():
    patterns = [_ins(0, 60), _ins(70, 60)]
    assert _component_sizes(patterns, "none") == [1, 1]
    assert _component_sizes(patterns, False) == [1, 1]


def test_same_type_never_joins_an_insertion_with_a_deletion():
    patterns = [_ins(0, 60), _del(100, 80), _ins(900, 70)]
    assert _component_sizes(patterns, "same-type") == [1, 2]


def test_adjacent_joins_only_nearby_fragments_of_the_same_type():
    # 0-60 and 90-150 are 30 bp apart: fragments of one event.
    # 900-970 is 750 bp further on: a distinct event in the same repeat.
    patterns = [_ins(0, 60), _ins(90, 60), _ins(900, 70)]
    assert _component_sizes(patterns, "adjacent", max_gap=50) == [1, 2]
    assert _component_sizes(patterns, "adjacent", max_gap=1000) == [3]


def test_adjacent_does_not_join_across_types():
    patterns = [_ins(0, 60), _del(70, 80)]
    assert _component_sizes(patterns, "adjacent", max_gap=50) == [1, 1]


def test_unknown_mode_is_rejected():
    with pytest.raises(ValueError, match="collapse_repeats"):
        _component_sizes([_ins(0, 60), _ins(70, 60)], "sometimes")


def test_proximal_joins_nearby_fragments_without_a_shared_repeat():
    first, second, far = _ins(0, 33), _ins(49, 13), _ins(900, 70)
    for p in (first, second, far):
        for prim in p.SVprimitives:
            prim.repeatIDs = []
    patterns = [first, second, far]
    # outside any repeat, adjacent has nothing to join on
    assert _component_sizes(patterns, "adjacent") == [1, 1, 1]
    assert _component_sizes(patterns, "proximal") == [1, 2]
