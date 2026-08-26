"""Unit tests for the vertical-merge size-similarity gate.

The gate is two-armed and the arms are OR-ed, so the *weaker* arm decides:

  1) a fractional bound on the two sizes;
  2) a population-driven Cohen's D on the background size-distortion signals.

The same gate is reached from three entry points -- insertions, deletions and
inversions -- and the tests below deliberately exercise all three with the same
inputs, because the gate was once triplicated by copy-paste and drifted.
"""

import pickle
import random

import numpy as np
import pytest

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.svcalling import genotyping, svcomposite_merging
from svirlpool.svcalling.SVcomposite import SVcomposite
from svirlpool.svcalling.svcomposite_merging import (
    can_merge_svComposites_deletions,
    can_merge_svComposites_insertions,
    can_merge_svComposites_inversions,
    sizetolerance_from_SVcomposite,
)
from svirlpool.svcalling.svcomposite_utils import cohens_d

# ---------------------------------------------------------------------------
# Helper factories
# ---------------------------------------------------------------------------


def _make_genotype(reads: list[str] | None = None) -> genotyping.GenotypeMeasurement:
    if reads is None:
        reads = ["read1", "read2", "read3"]
    return genotyping.GenotypeMeasurement(
        start_on_consensus=0,
        supporting_reads_start=reads,
    )


def _make_svprimitive_ins(
    *,
    chr: str = "chr1",
    ref_start: int = 1000,
    ref_end: int = 1001,
    read_start: int = 0,
    read_end: int = 500,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    reads: list[str] | None = None,
) -> SVprimitives.SVprimitive:
    return SVprimitives.SVprimitive(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        size=abs(read_end - read_start),
        sv_type=0,  # INS
        chr=chr,
        repeatIDs=[],
        original_alt_sequences=["A" * abs(read_end - read_start)],
        original_ref_sequences=[],
        samplename=samplename,
        consensusID=consensusID,
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=(chr, ref_start - 500, ref_end + 500),
        genotypeMeasurement=_make_genotype(reads),
    )


def _make_svprimitive_del(
    *,
    chr: str = "chr1",
    ref_start: int = 1000,
    ref_end: int = 1500,
    read_start: int = 0,
    read_end: int = 10,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    reads: list[str] | None = None,
) -> SVprimitives.SVprimitive:
    return SVprimitives.SVprimitive(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        size=abs(ref_end - ref_start),
        sv_type=1,  # DEL
        chr=chr,
        repeatIDs=[],
        original_alt_sequences=[],
        original_ref_sequences=["A" * abs(ref_end - ref_start)],
        samplename=samplename,
        consensusID=consensusID,
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=(chr, ref_start - 500, ref_end + 500),
        genotypeMeasurement=_make_genotype(reads),
    )


def _make_insertion_composite(
    *,
    size: int = 500,
    size_distortions: dict[str, float] | None = None,
    chr: str = "chr1",
    ref_start: int = 1000,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    sequence: str | None = None,
    reads: list[str] | None = None,
) -> SVcomposite:
    """Create an insertion SVcomposite with controllable size and size_distortions."""
    svp = _make_svprimitive_ins(
        chr=chr,
        ref_start=ref_start,
        ref_end=ref_start + 1,
        read_start=0,
        read_end=size,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    pattern = SVpatterns.SVpatternInsertion(
        SVprimitives=[svp],
        size_distortions=size_distortions,
    )
    if sequence is None:
        sequence = "A" * size
    pattern.set_sequence(sequence)
    return SVcomposite.from_SVpattern(pattern)


def _make_deletion_composite(
    *,
    size: int = 500,
    size_distortions: dict[str, float] | None = None,
    chr: str = "chr1",
    ref_start: int = 1000,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    sequence: str | None = None,
    reads: list[str] | None = None,
) -> SVcomposite:
    """Create a deletion SVcomposite with controllable size and size_distortions."""
    svp = _make_svprimitive_del(
        chr=chr,
        ref_start=ref_start,
        ref_end=ref_start + size,
        read_start=0,
        read_end=10,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    pattern = SVpatterns.SVpatternDeletion(
        SVprimitives=[svp],
        size_distortions=size_distortions,
    )
    if sequence is None:
        sequence = "A" * size
    pattern.set_sequence(sequence)
    return SVcomposite.from_SVpattern(pattern)


def _random_dna(length: int, seed: int) -> str:
    """A high-complexity DNA string, reproducible from `seed`.

    The merge criterion is parameterised by sequence complexity, so a test that
    leaves the sequence at the default homopolymer is not testing the threshold
    it names -- it is testing the maximum-complexity-allowance regime. Use this
    wherever the intent is "well-resolved sequence"; use an explicit "A" * n
    wherever the intent is a repeat.

    Note that complexity is only computed from the sequence for lengths up to
    `sequence_complexity_max_length` (300 bp, see SVpatterns.set_sequence).
    Above that a *dummy all-ones* track is stored instead, which reads as
    maximum complexity regardless of the actual sequence, so keep sizes at or
    below 300 in tests that exercise complexity.
    """
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _make_svprimitive_inv(
    *,
    chr: str = "chr1",
    ref_start: int = 1000,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    reads: list[str] | None = None,
) -> SVprimitives.SVprimitive:
    return SVprimitives.SVprimitive(
        ref_start=ref_start,
        ref_end=ref_start + 1,
        read_start=0,
        read_end=1,
        size=1,
        sv_type=3,  # breakend-like primitive, four of them make an inversion
        chr=chr,
        repeatIDs=[],
        original_alt_sequences=[],
        original_ref_sequences=[],
        samplename=samplename,
        consensusID=consensusID,
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=(chr, ref_start - 500, ref_start + 500),
        genotypeMeasurement=_make_genotype(reads),
    )


def _make_inversion_composite(
    *,
    size: int = 500,
    size_distortions: dict[str, float] | None = None,
    chr: str = "chr1",
    ref_start: int = 1000,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    sequence: str | None = None,
    reads: list[str] | None = None,
) -> SVcomposite:
    """Create an inversion SVcomposite. Inversions need exactly four primitives.

    `SVpatternInversion.get_size()` measures the *inner* interval, primitive[1]
    to primitive[2], so the two inner primitives carry the size.
    """
    outer_left = _make_svprimitive_inv(
        chr=chr,
        ref_start=ref_start,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    inner_left = _make_svprimitive_inv(
        chr=chr,
        ref_start=ref_start,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    inner_right = _make_svprimitive_inv(
        chr=chr,
        ref_start=ref_start + size,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    outer_right = _make_svprimitive_inv(
        chr=chr,
        ref_start=ref_start + size,
        samplename=samplename,
        consensusID=consensusID,
        reads=reads,
    )
    pattern = SVpatterns.SVpatternInversion(
        SVprimitives=[outer_left, inner_left, inner_right, outer_right],
        size_distortions=size_distortions,
    )
    if sequence is None:
        sequence = "A" * size
    pattern.set_sequence(sequence)
    return SVcomposite.from_SVpattern(pattern)


# The three entry points into the shared size gate, keyed by SV type. Anything
# that claims "the size gate does X" should be asserted through all three.
_KINDS = ("insertion", "deletion", "inversion")

_MAKERS = {
    "insertion": _make_insertion_composite,
    "deletion": _make_deletion_composite,
    "inversion": _make_inversion_composite,
}

_MERGERS = {
    "insertion": can_merge_svComposites_insertions,
    "deletion": can_merge_svComposites_deletions,
    "inversion": can_merge_svComposites_inversions,
}


def _make_composite(kind: str, **kwargs) -> SVcomposite:
    return _MAKERS[kind](**kwargs)


def _size_gate(
    kind: str,
    size_a: int,
    size_b: int,
    *,
    tolerance: float,
    scale_by_complexity_factor: float = 0.0,
    d: float = 2.0,
    size_distortions: dict[str, float] | None = None,
    sequence_a: str | None = None,
    sequence_b: str | None = None,
) -> bool:
    """Merge decision with every gate other than the size gate neutralised.

    `near` is set far beyond both events and `min_kmer_overlap` to 0.0, so the
    proximity and k-mer arms always pass and the return value *is* the size
    gate's verdict. With `size_distortions=None` the population arm is off as
    well (empty populations are treated as "not similar"), so the verdict is
    exactly the fractional arm.
    """
    a = _make_composite(
        kind,
        size=size_a,
        size_distortions=size_distortions,
        samplename="sample1",
        consensusID="1.0",
        sequence=sequence_a,
    )
    b = _make_composite(
        kind,
        size=size_b,
        size_distortions=size_distortions,
        samplename="sample2",
        consensusID="2.0",
        sequence=sequence_b,
    )
    return _MERGERS[kind](
        a=a,
        b=b,
        apriori_size_difference_fraction_tolerance=tolerance,
        d=d,
        near=10 * max(size_a, size_b, 100),
        min_kmer_overlap=0.0,
        scale_by_complexity_factor=scale_by_complexity_factor,
    )


# Size pairs spanning the plausible range, from near-identical to wildly
# dissimilar. (50, 5000) is the catalogue's example of a pair that no size
# criterion should ever call similar.
_SIZE_PAIRS = [
    (50, 5000),
    (12, 4096),
    (100, 200),
    (100, 250),
    (250, 1000),
    (900, 1000),
    (500, 525),
    (300, 301),
]


# ===========================================================================
# INSERTION TESTS
# ===========================================================================


class TestCanMergeInsertions:
    """Tests for can_merge_svComposites_insertions."""

    def test_identical_insertions_merge(self):
        """Two identical-size insertions at the same locus with similar populations should merge."""
        a = _make_insertion_composite(
            size=500,
            size_distortions={"r1": 5, "r2": -3, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=500,
            size_distortions={"r1": 4, "r2": -2, "r3": 1},
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )

    def test_similar_size_within_fraction_tolerance(self):
        """Two insertions within 10% size difference should merge via fractional test."""
        # 500 vs 540 → diff=40, 10% of 540=54 → within tolerance
        a = _make_insertion_composite(
            size=500,
            size_distortions={"r1": 10, "r2": -10},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=540,
            size_distortions={"r1": 10, "r2": -10},
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )

    def test_different_sizes_beyond_fraction_but_populations_overlap(self):
        """Sizes differ by >10%, but populations have overlapping distributions → merge via Cohen's D."""
        # sizes: 500 vs 600 → diff/max = 100/600 ≈ 16.7% → fraction test fails
        # but populations with large spread should overlap and have small Cohen's D
        a = _make_insertion_composite(
            size=500,
            size_distortions={
                f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])
            },
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=600,
            size_distortions={
                f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])
            },
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )

    def test_very_different_sizes_reject(self):
        """Very different sizes are rejected when the sequence is well resolved.

        Sizes 100 vs 200 (a 2x difference) fail the fraction test, which requires
        the relative size difference to be within `tol` (here 10%). With a
        high-complexity sequence the complexity allowance is close to zero, so the
        shifted Cohen's D test fails as well and the merge is correctly rejected.

        The sequence matters: complexity is what decides how much size disagreement
        is tolerated (see test_very_different_sizes_merge_in_low_complexity for the
        same sizes in a homopolymer).
        """
        a = _make_insertion_composite(
            size=100,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(100, seed=1),
        )
        b = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(200, seed=2),
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=300,
            min_kmer_overlap=0.0,
            scale_by_complexity_factor=1.0,
        )

    def test_very_different_sizes_merge_in_low_complexity(self):
        """The same sizes DO merge inside a homopolymer, and that is intended.

        Complexity is a proxy for how much placement and size ambiguity the aligner
        introduces at a locus. In a poly-A run it has near-total freedom in both,
        so `sizetolerance_from_SVcomposite` grants a tolerance approaching the full
        event size and a 100 bp and a 200 bp insertion at the same position are
        treated as one VNTR allele family.

        This is the mirror image of test_very_different_sizes_reject: identical
        sizes, identical thresholds, opposite outcome, decided only by sequence
        complexity. Setting `scale_by_complexity_factor=0.0` withdraws the
        allowance and restores rejection.
        """
        a = _make_insertion_composite(
            size=100,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
            sequence="A" * 100,
        )
        b = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
            sequence="A" * 200,
        )
        kwargs = {
            "a": a,
            "b": b,
            "apriori_size_difference_fraction_tolerance": 0.1,
            "d": 2.0,
            "near": 300,
            "min_kmer_overlap": 0.0,
        }
        assert can_merge_svComposites_insertions(
            **kwargs, scale_by_complexity_factor=1.0
        )
        assert not can_merge_svComposites_insertions(
            **kwargs, scale_by_complexity_factor=0.0
        )

    def test_identical_sizes_merge_regardless_of_complexity(self):
        """Regression guard: a complexity difference is not a size difference.

        Until this was fixed the gate compared complexity-ADJUSTED sizes,
        `lerp(size, size * complexity, scale)`. Because the two composites'
        complexity tracks are estimated from different consensus sequences in
        different samples, that turned a complexity difference into a size
        difference out of nothing: two events of IDENTICAL size were rejected once
        their complexity estimates differed by more than about 0.090 — scale
        invariantly, and hardest inside the repeats where the estimates diverge
        most, which is precisely where merging matters.

        Same size, maximally different sequence complexity, must merge.
        """
        a = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(200, seed=3),
        )
        b = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
            sequence="A" * 200,
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=300,
            min_kmer_overlap=0.0,
            scale_by_complexity_factor=1.0,
        )

    def test_empty_populations_fraction_pass(self):
        """When populations are empty, only the fractional test is used. Similar sizes merge."""
        a = _make_insertion_composite(
            size=500,
            size_distortions=None,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=520,
            size_distortions=None,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )

    def test_empty_populations_fraction_fail(self):
        """With no populations, a large size difference fails the fraction test and is rejected.

        With no size_distortions, population_similar=False (empty populations). The
        corrected fraction test also fails because a 100 vs 200 size difference (100%)
        exceeds the 10% tolerance, so similar_size=False and the merge is correctly
        rejected.
        """
        a = _make_insertion_composite(
            size=100,
            size_distortions=None,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=200,
            size_distortions=None,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=300,
            min_kmer_overlap=0.0,
            scale_by_complexity_factor=1.0,
        )

    def test_not_near_rejects(self):
        """Insertions on different chromosomes should not merge."""
        a = _make_insertion_composite(
            size=500,
            size_distortions={"r1": 5},
            chr="chr1",
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=500,
            size_distortions={"r1": 5},
            chr="chr2",
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )

    def test_kmer_similarity_rejects(self):
        """Same-size insertions with different sequence content should be rejected by k-mer check."""
        a = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 5, "r2": -3},
            samplename="sample1",
            consensusID="1.0",
            sequence="ATCGATCG" * 25,  # 200bp AT-rich
        )
        b = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 5, "r2": -3},
            samplename="sample2",
            consensusID="2.0",
            sequence="GCGCGCGC" * 25,  # 200bp GC-rich
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )

    def test_strict_tolerance_rejects_borderline(self):
        """A borderline 4.8% size difference merges under a lenient tolerance but is rejected under a strict one."""
        # 500 vs 525 → diff/max = 25/525 ≈ 4.8%
        # The corrected check requires diff <= tol * max_size, so whether this pair
        # merges now genuinely depends on the tolerance value, as the test name implies.
        a = _make_insertion_composite(
            size=500,
            size_distortions={"r1": 1, "r2": -1},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=525,
            size_distortions={"r1": 1, "r2": -1},
            samplename="sample2",
            consensusID="2.0",
        )
        # With 10% tolerance: 4.8% difference is within tolerance -> merges
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )
        # With 1% tolerance: 4.8% difference exceeds tolerance -> rejected
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.01,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )


# ===========================================================================
# DELETION TESTS
# ===========================================================================


class TestCanMergeDeletions:
    """Tests for can_merge_svComposites_deletions."""

    def test_identical_deletions_merge(self):
        """Two identical-size deletions at the same locus should merge."""
        a = _make_deletion_composite(
            size=500,
            size_distortions={"r1": 5, "r2": -3, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=500,
            size_distortions={"r1": 4, "r2": -2, "r3": 1},
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_similar_size_within_fraction_tolerance(self):
        """Deletions within 10% size difference should merge via fractional test."""
        # 1000 vs 1080 → diff=80, 10% of 1080=108 → within tolerance
        a = _make_deletion_composite(
            size=1000,
            size_distortions={"r1": 10, "r2": -10},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=1080,
            size_distortions={"r1": 10, "r2": -10},
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_different_sizes_beyond_fraction_but_populations_overlap(self):
        """Sizes differ by >10%, but wide population distributions → merge via Cohen's D."""
        # 500 vs 600 → 16.7% → fraction fails
        # wide spread populations should produce small Cohen's D
        a = _make_deletion_composite(
            size=500,
            size_distortions={
                f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])
            },
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=600,
            size_distortions={
                f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])
            },
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_very_different_sizes_reject(self):
        """A 150% size difference is rejected when the reference is well resolved.

        Sizes 100 vs 250 fail the fraction test (the relative size difference far
        exceeds the 10% tolerance). With a high-complexity reference span the
        complexity allowance is small, so the shifted Cohen's D test fails too and
        the merge is correctly rejected.

        Deletions take their complexity from the REFERENCE span rather than an
        assembled sequence (`get_reference_complexity_tracks`), which is the right
        substrate: it is the reference the aligner is placing the event against.
        """
        a = _make_deletion_composite(
            size=100,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(100, seed=11),
        )
        b = _make_deletion_composite(
            size=250,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(250, seed=12),
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=350,
            min_kmer_overlap=0.0,
        )

    def test_very_different_sizes_merge_in_low_complexity(self):
        """The same sizes merge inside a low-complexity reference span.

        The mirror of the test above, and the deletion counterpart of
        TestCanMergeInsertions.test_very_different_sizes_merge_in_low_complexity.
        In a homopolymer the aligner's choice of deletion boundaries — and
        therefore of size — is close to arbitrary, so a 100 bp and a 250 bp
        deletion at the same position are treated as one allele family.

        This is also the regime where the k-mer test cannot help: a smaller
        deletion nested in a larger one at the same locus is k-mer identical by
        construction, so size is the only criterion left, and it is the one that
        complexity relaxes here.
        """
        a = _make_deletion_composite(
            size=100,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
            sequence="A" * 100,
        )
        b = _make_deletion_composite(
            size=250,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
            sequence="A" * 250,
        )
        kwargs = {
            "a": a,
            "b": b,
            "apriori_size_difference_fraction_tolerance": 0.1,
            "d": 2.0,
            "near": 350,
            "min_kmer_overlap": 0.0,
        }
        assert can_merge_svComposites_deletions(
            **kwargs, scale_by_complexity_factor=1.0
        )
        assert not can_merge_svComposites_deletions(
            **kwargs, scale_by_complexity_factor=0.0
        )

    def test_empty_populations_fraction_pass(self):
        """Empty populations fall back to fractional test only. Similar sizes merge."""
        a = _make_deletion_composite(
            size=500,
            size_distortions=None,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=530,
            size_distortions=None,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_empty_populations_fraction_fail(self):
        """With no populations, a large size difference fails the fraction test and is rejected.

        With no size_distortions, population_similar=False (empty populations). The
        corrected fraction test also fails because a 100 vs 200 size difference (100%)
        exceeds the 10% tolerance, so similar_size=False and the merge is correctly
        rejected.
        """
        a = _make_deletion_composite(
            size=100,
            size_distortions=None,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=200,
            size_distortions=None,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=300,
            min_kmer_overlap=0.0,
        )

    def test_not_near_rejects(self):
        """Deletions on different chromosomes should not merge."""
        a = _make_deletion_composite(
            size=500,
            size_distortions={"r1": 5},
            chr="chr1",
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=500,
            size_distortions={"r1": 5},
            chr="chr2",
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_kmer_similarity_rejects_deletions(self):
        """Same-size deletions with different reference sequence content should be rejected."""
        a = _make_deletion_composite(
            size=200,
            size_distortions={"r1": 5, "r2": -3},
            samplename="sample1",
            consensusID="1.0",
            sequence="ATCGATCG" * 25,  # 200bp AT-rich
        )
        b = _make_deletion_composite(
            size=200,
            size_distortions={"r1": 5, "r2": -3},
            samplename="sample2",
            consensusID="2.0",
            sequence="GCGCGCGC" * 25,  # 200bp GC-rich
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_strict_cohens_d_threshold(self):
        """With strict Cohen's D and strict fraction tolerance together, a borderline size difference is rejected."""
        # 500 vs 560 → diff/max = 60/560 ≈ 10.7% → fails both the 10% and the 5% fraction tolerance
        # wide population spread so Cohen's D is moderate (~1.1)
        a = _make_deletion_composite(
            size=500,
            size_distortions={
                f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])
            },
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=560,
            size_distortions={
                f"r{i}": v for i, v in enumerate([-80, -50, -20, 0, 20, 50, 80])
            },
            samplename="sample2",
            consensusID="2.0",
        )
        # With lenient d=2.0: fraction test fails (10.7% > 10%), but Cohen's D ~1.1 is
        # within the lenient threshold -> merge via the population test
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )
        # With strict d=0.5 and a strict 5% fraction tolerance: both tests fail
        # (10.7% > 5%, and Cohen's D ~1.1 > 0.5) -> the merge is rejected
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.05,
            d=0.5,
            near=150,
            min_kmer_overlap=0.7,
        )


# ===========================================================================
# THE SIZE GATE ITSELF
#
# The classes below test the shared size-similarity gate directly, through all
# three entry points, with every other merge criterion neutralised.
# ===========================================================================


class TestSizeGateIsNoLongerInert:
    """The former tautology, F3, restated as the property it has become.

    Before the fix the fractional arm read

        |a_adj - b_adj| <= (1.0 + tol) * max(a_adj, b_adj)

    and `get_size()` returns a magnitude, so `|x - y| <= max(x, y)` already held
    for every non-negative pair -- before the `1.0 +` even added its bonus. The
    arm was therefore true whenever `max_size > 0`, at every tolerance, which
    made the whole size gate inert and the Cohen's-D arm it is OR-ed with dead
    code.

    The predecessor of this class asserted that every pair in `_SIZE_PAIRS`
    merged at every tolerance, and passed. The same sweep is kept here with the
    verdict inverted, so the sweep that documented the defect now documents its
    absence. The old assertion survives verbatim as the control setting,
    TestSizeGateControlSetting, where tol = 1.0 makes it true again by design.
    """

    def test_dissimilar_sizes_no_longer_merge_at_a_strict_tolerance(self):
        """At tol = 0.06 only pairs genuinely within 6% may pass."""
        wrong = []
        for kind in _KINDS:
            for size_a, size_b in _SIZE_PAIRS:
                want = abs(size_a - size_b) <= 0.06 * max(size_a, size_b)
                got = _size_gate(kind, size_a, size_b, tolerance=0.06)
                if got is not want:
                    wrong.append((kind, size_a, size_b, got, want))
        assert wrong == [], f"the gate does not follow its own arithmetic: {wrong}"


class TestSizeGateControlSetting:
    """`--apriori-size-difference-fraction-tolerance 1.0` is vacuous by construction.

    This is the control arm the re-benchmark needs: at tol = 1.0 the corrected
    bound `|a - b| <= 1.0 * max(a, b)` is exactly the inequality that made the
    pre-fix gate inert, so the fixed code with tol = 1.0 reproduces the pre-fix
    behaviour and any measured difference is attributable to the tolerance
    rather than to the code change.

    It must therefore hold both before and after the fix.
    """

    def test_tolerance_one_accepts_every_size_pair(self):
        failures = []
        for kind in _KINDS:
            for size_a, size_b in _SIZE_PAIRS:
                if not _size_gate(kind, size_a, size_b, tolerance=1.0):
                    failures.append((kind, size_a, size_b))
        assert failures == [], f"tol=1.0 must be vacuous, but it rejected: {failures}"

    def test_tolerance_one_is_vacuous_at_full_complexity_weight(self):
        """The control must not depend on the complexity weight either."""
        failures = []
        for kind in _KINDS:
            for size_a, size_b in _SIZE_PAIRS:
                if not _size_gate(
                    kind,
                    size_a,
                    size_b,
                    tolerance=1.0,
                    scale_by_complexity_factor=1.0,
                ):
                    failures.append((kind, size_a, size_b))
        assert failures == [], f"tol=1.0 must be vacuous, but it rejected: {failures}"


class TestSizeGateIsSharedByAllSVTypes:
    """Regression guard for F6: one gate, three entry points, one verdict.

    Insertions, deletions and inversions each used to carry their own copy of
    the two-armed test. That is how a single defect came to exist in three
    places at once, and how the three copies could drift apart. Whatever the
    gate decides, it must decide identically for all three.
    """

    def test_all_three_paths_agree_on_every_size_pair(self):
        disagreements = []
        for size_a, size_b in _SIZE_PAIRS:
            for tolerance in (0.0, 0.01, 0.06, 0.1, 0.5, 1.0):
                verdicts = {
                    kind: _size_gate(kind, size_a, size_b, tolerance=tolerance)
                    for kind in _KINDS
                }
                if len(set(verdicts.values())) != 1:
                    disagreements.append((size_a, size_b, tolerance, verdicts))
        assert disagreements == [], f"the size gate is not shared: {disagreements}"

    def test_all_three_paths_agree_with_populations_and_complexity(self):
        """Same, with both arms live and the complexity term at full weight."""
        distortions = {f"r{i}": v for i, v in enumerate([-30, -10, 0, 10, 30])}
        disagreements = []
        for size_a, size_b in [(100, 200), (200, 205), (100, 250)]:
            for tolerance in (0.0, 0.06, 0.5):
                for sequence in ("A", "R"):
                    seq_a = (
                        "A" * size_a
                        if sequence == "A"
                        else _random_dna(size_a, seed=100 + size_a)
                    )
                    seq_b = (
                        "A" * size_b
                        if sequence == "A"
                        else _random_dna(size_b, seed=200 + size_b)
                    )
                    verdicts = {
                        kind: _size_gate(
                            kind,
                            size_a,
                            size_b,
                            tolerance=tolerance,
                            scale_by_complexity_factor=1.0,
                            size_distortions=distortions,
                            sequence_a=seq_a,
                            sequence_b=seq_b,
                        )
                        for kind in _KINDS
                    }
                    if len(set(verdicts.values())) != 1:
                        disagreements.append((
                            size_a,
                            size_b,
                            tolerance,
                            sequence,
                            verdicts,
                        ))
        assert disagreements == [], f"the size gate is not shared: {disagreements}"


class TestFractionArmIsATrueFractionalBound:
    """The corrected arm: `|size_a - size_b| <= tol * max(size_a, size_b)`."""

    def test_dissimilar_sizes_are_rejected(self):
        """Pairs no tolerance below 1.0 should accept."""
        rejected_everywhere = []
        for kind in _KINDS:
            for size_a, size_b in [(50, 5000), (12, 4096), (100, 250)]:
                for tolerance in (0.0, 0.01, 0.06, 0.1):
                    rejected_everywhere.append(
                        _size_gate(kind, size_a, size_b, tolerance=tolerance)
                    )
        assert not any(rejected_everywhere)

    def test_boundary_lands_where_the_arithmetic_says(self):
        """900 vs 1000: diff = 100, max = 1000, so the boundary is exactly 0.1.

        The bound is inclusive (`<=`), so tol = 0.1 must accept and anything
        below it must reject. Asserted for all three entry points.
        """
        size_a, size_b = 900, 1000
        expected = {
            0.0: False,
            0.05: False,
            0.09: False,
            0.0999: False,
            0.1: True,
            0.100001: True,
            0.5: True,
            1.0: True,
        }
        wrong = []
        for kind in _KINDS:
            for tolerance, want in expected.items():
                got = _size_gate(kind, size_a, size_b, tolerance=tolerance)
                if got is not want:
                    wrong.append((kind, tolerance, got, want))
        assert wrong == [], f"boundary is not at 0.1: {wrong}"

    def test_bound_is_scale_invariant(self):
        """The same relative difference decides the same way at any absolute size."""
        wrong = []
        for kind in _KINDS:
            for base in (100, 1000, 10000):
                # 20% apart: rejected at 0.1, accepted at 0.25
                small, large = int(base * 0.8), base
                if _size_gate(kind, small, large, tolerance=0.1):
                    wrong.append(("accepted at 0.1", kind, small, large))
                if not _size_gate(kind, small, large, tolerance=0.25):
                    wrong.append(("rejected at 0.25", kind, small, large))
        assert wrong == [], f"the bound is not scale invariant: {wrong}"


class TestComplexityIsATolerance:
    """F4: complexity must GRANT size tolerance, in bp, not rescale the sizes.

    `sizetolerance_from_SVcomposite` returns the allowance in base pairs:
    `(1 - mean_complexity) * |size|`. Low complexity -- a homopolymer, a perfect
    repeat -- is where the aligner has the most freedom in where it places an
    indel and how large it calls it, so it must grant the *most* tolerance.
    """

    def test_low_complexity_grants_more_tolerance_than_high(self):
        homopolymer = _make_insertion_composite(size=200, sequence="A" * 200)
        well_resolved = _make_insertion_composite(
            size=200, sequence=_random_dna(200, seed=42)
        )

        tol_low_complexity = sizetolerance_from_SVcomposite(homopolymer)
        tol_high_complexity = sizetolerance_from_SVcomposite(well_resolved)

        assert tol_low_complexity > tol_high_complexity, (
            f"homopolymer allowance {tol_low_complexity:.1f} bp should exceed "
            f"well-resolved allowance {tol_high_complexity:.1f} bp"
        )
        # An allowance in bp, bounded by the event's own size.
        assert 0.0 <= tol_high_complexity <= 200.0
        assert 0.0 <= tol_low_complexity <= 200.0

    def test_the_two_no_evidence_regimes_are_distinguished(self):
        """A present all-zero track is minimum complexity, not missing data.

        Those two situations are opposite and the fix deliberately separates
        them:

        * track ABSENT -- no evidence of aligner ambiguity. Fail closed:
          mean_complexity = 1.0, allowance 0 bp.
        * track PRESENT and all zero -- a homopolymer or a perfect repeat, i.e.
          *minimum* complexity, which is exactly where placement is most
          arbitrary. Allowance = the full event size.
        """
        absent = _make_insertion_composite(size=200)
        # Built by the factory, then stripped: no `set_sequence`, no track.
        absent.svPatterns[0].sequence_complexity = None
        assert sizetolerance_from_SVcomposite(absent) == 0.0

        all_zero = _make_insertion_composite(size=200)
        all_zero.svPatterns[0].sequence_complexity = pickle.dumps(
            np.zeros(200, dtype=np.float32)
        )
        assert sizetolerance_from_SVcomposite(all_zero) == pytest.approx(200.0)

    def test_events_above_300_bp_get_no_allowance_whatever_their_sequence(self):
        """A documented limitation, pinned so it cannot change unnoticed.

        `set_sequence` computes a complexity track only for sequences up to
        `sequence_complexity_max_length` (300 bp). Above that it stores a *dummy
        all-ones* array instead -- a track that is present and well formed, and
        that reads as maximum complexity. The allowance is therefore exactly
        0 bp for every event above 300 bp, whatever its sequence actually is: a
        5 kb poly-A insertion is treated as perfectly well-resolved sequence.

        This is not what the complexity term is meant to express, and it is not
        F4; it is a property of where the track is computed. Assert it so that
        anyone who moves the 300 bp threshold sees this test rather than an
        unexplained change in merge behaviour.
        """
        long_homopolymer = _make_insertion_composite(size=1000, sequence="A" * 1000)
        assert sizetolerance_from_SVcomposite(long_homopolymer) == 0.0

        short_homopolymer = _make_insertion_composite(size=300, sequence="A" * 300)
        assert sizetolerance_from_SVcomposite(short_homopolymer) > 290.0

    def test_low_complexity_merges_what_high_complexity_rejects(self):
        """The sign of the effect, asserted on merge decisions.

        Identical sizes, identical thresholds, identical populations; only the
        sequence differs. A 100 bp and a 200 bp insertion inside a poly-A run
        are one VNTR allele family and must merge; the same two events in
        well-resolved sequence must not.
        """
        distortions = {"r1": 1, "r2": -1, "r3": 2}
        low_complexity = _size_gate(
            "insertion",
            100,
            200,
            tolerance=0.1,
            scale_by_complexity_factor=1.0,
            size_distortions=distortions,
            sequence_a="A" * 100,
            sequence_b="A" * 200,
        )
        high_complexity = _size_gate(
            "insertion",
            100,
            200,
            tolerance=0.1,
            scale_by_complexity_factor=1.0,
            size_distortions=distortions,
            sequence_a=_random_dna(100, seed=1),
            sequence_b=_random_dna(200, seed=2),
        )
        assert low_complexity, "a homopolymer must grant tolerance"
        assert not high_complexity, "well-resolved sequence must not"

    def test_scale_by_complexity_factor_withdraws_the_allowance(self):
        """`--scale-by-complexity-factor 0.0` must turn the allowance off."""
        distortions = {"r1": 1, "r2": -1, "r3": 2}
        kwargs = {
            "tolerance": 0.1,
            "size_distortions": distortions,
            "sequence_a": "A" * 100,
            "sequence_b": "A" * 200,
        }
        assert _size_gate(
            "insertion", 100, 200, scale_by_complexity_factor=1.0, **kwargs
        )
        assert not _size_gate(
            "insertion", 100, 200, scale_by_complexity_factor=0.0, **kwargs
        )

    def test_a_complexity_difference_is_not_a_size_difference(self):
        """Regression guard: two events of the SAME size must always merge.

        `main` compared complexity-ADJUSTED sizes,
        `lerp(size, size * complexity, scale)`. The two composites' complexity
        tracks come from different consensus sequences in different samples, so
        a difference between the two estimates was indistinguishable from a
        difference in size and could only ever suppress merges -- hardest inside
        the repeats where the estimates diverge most.

        This passes on `main` only because the gate is inert there (F3); it
        becomes load-bearing the moment the fractional bound is corrected, which
        is why it is asserted at tol = 0.0, where nothing but an exact size
        match may pass.
        """
        for kind in _KINDS:
            assert _size_gate(
                kind,
                200,
                200,
                tolerance=0.0,
                scale_by_complexity_factor=1.0,
                sequence_a=_random_dna(200, seed=3),
                sequence_b="A" * 200,
            ), f"{kind}: identical sizes must merge whatever the complexity"


# ===========================================================================
# F2 -- A DEGENERATE SIZE POPULATION IS NOT AN EFFECT SIZE
#
# `cohens_d` used to fabricate a value whenever the pooled standard deviation
# was zero: `inf` if the two means differed, `0.0` if they did not. Both are
# statements about the separation of two distributions *relative to their
# spread*, made from data that carry no information about spread at all. The
# `inf` in particular reads as "maximally separated, therefore reject", and it
# is what let F1 -- every distortion value stuck at exactly 0.0, i.e. every
# population a constant vector -- pass through an entire benchmark campaign
# without a single visible symptom.
#
# The contract now is: the utility reports "not estimable" (None) and the
# caller decides what to do about it.
# ===========================================================================


class TestCohensDIsNotEstimableOnDegenerateInput:
    """Direct unit tests for the effect-size utility itself."""

    def test_two_constant_populations_with_different_means(self):
        """Both groups constant, means apart.

        Pre-fix this returned `float("inf")`: an infinitely large effect size
        inferred from zero observed variability. That is the value that made the
        population arm reject silently everywhere F1 was active.
        """
        assert cohens_d([500] * 8, [600] * 8) is None

    def test_two_constant_populations_with_equal_means(self):
        """Both groups constant, means identical.

        Pre-fix this returned `0.0`. It is no more defensible than the `inf`:
        the quotient is 0/0, and 0.0 asserts "no difference relative to the
        spread" from data with no spread. It merely fails in the permissive
        direction instead of the restrictive one, which is why nobody noticed.
        """
        assert cohens_d([500] * 8, [500] * 8) is None

    def test_single_observation_in_each_group(self):
        """One observation each: no within-group variability can be estimated.

        Pre-fix: `inf` for differing values, `0.0` for equal ones. This is not a
        contrived input -- it is what two composites with one supporting read
        each produce.
        """
        assert cohens_d([500], [600]) is None
        assert cohens_d([500], [500]) is None

    def test_a_constant_group_against_a_spread_group_is_estimable(self):
        """Degeneracy is a property of the *pair*, not of one group.

        With one group constant and the other spread, the pooled standard
        deviation is positive and the effect size is a real number.
        """
        value = cohens_d([500] * 8, [400, 450, 500, 550, 600, 650, 700, 750])
        assert value is not None
        assert np.isfinite(value)

    def test_a_genuine_effect_size_is_unchanged(self):
        """Regression guard: the non-degenerate arithmetic must not move."""
        x = [10.0, 12.0, 14.0, 16.0, 18.0]
        y = [20.0, 22.0, 24.0, 26.0, 28.0]
        expected = (np.mean(x) - np.mean(y)) / np.sqrt(
            (
                (len(x) - 1) * np.std(x, ddof=1) ** 2
                + (len(y) - 1) * np.std(y, ddof=1) ** 2
            )
            / (len(x) + len(y) - 2)
        )
        assert cohens_d(x, y) == pytest.approx(float(expected))

    def test_empty_input_still_raises(self):
        """An empty sample is a caller bug, not a degenerate population."""
        with pytest.raises(ValueError):
            cohens_d([], [1, 2, 3])
        with pytest.raises(ValueError):
            cohens_d([1, 2, 3], [])


class TestDegeneratePopulationsDoNotDecideMerges:
    """The same thing seen through the size gate.

    Which inputs actually reach the degenerate branch is narrower than it looks,
    because arm 2 short-circuits to `population_similar = True` whenever the two
    population means are already within the granted complexity tolerance. What
    is left is: *both* populations constant AND the size gap strictly greater
    than that tolerance.
    """

    @pytest.fixture(autouse=True)
    def _reset_rate_limiter(self):
        svcomposite_merging._reset_degenerate_population_warnings()
        yield
        svcomposite_merging._reset_degenerate_population_warnings()

    @staticmethod
    def _pair(size_a: int, size_b: int, distortions_a, distortions_b):
        a = _make_insertion_composite(
            size=size_a,
            size_distortions=distortions_a,
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(size_a, seed=11),
        )
        b = _make_insertion_composite(
            size=size_b,
            size_distortions=distortions_b,
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(size_b, seed=12),
        )
        return a, b

    def test_degenerate_pair_is_rejected_and_warned_about(self, caplog):
        """Both populations constant, gap far beyond any tolerance.

        This is exactly the state F1 produces genome-wide: every distortion
        value 0.0, so every population a constant vector at the composite's own
        size. Pre-fix `cohensD` came back as `inf` and `population_similar` was
        False -- the right verdict reached by an invented number, and reached in
        complete silence.
        """
        a, b = self._pair(
            100,
            200,
            {"r1": 0.0, "r2": 0.0, "r3": 0.0},
            {"r1": 0.0, "r2": 0.0, "r3": 0.0},
        )
        with caplog.at_level(
            "WARNING", logger="svirlpool.svcalling.svcomposite_merging"
        ):
            similar, fraction_similar, population_similar, cohensD = (
                svcomposite_merging._similar_size(a, b, 0.1, 0.0, 2.0)
            )
        # `is False` is deliberately avoided throughout: arm 1 returns a numpy
        # bool, so identity against the Python singleton does not hold.
        assert not population_similar
        assert not fraction_similar
        assert not similar
        assert cohensD is None, "a non-estimable effect size must not be a number"
        assert any(
            "DEGENERATE_SIZE_POPULATION" in record.message for record in caplog.records
        ), "the degenerate population must be reported, not swallowed"

    def test_the_warning_is_rate_limited(self, caplog):
        """A whole-genome run hits this at every locus while F1 is unfixed."""
        with caplog.at_level(
            "WARNING", logger="svirlpool.svcalling.svcomposite_merging"
        ):
            for _ in range(svcomposite_merging.DEGENERATE_POPULATION_WARN_LIMIT + 20):
                a, b = self._pair(100, 200, {"r1": 0.0}, {"r1": 0.0})
                svcomposite_merging._similar_size(a, b, 0.1, 0.0, 2.0)
        per_locus = [
            r for r in caplog.records if "DEGENERATE_SIZE_POPULATION::" in r.message
        ]
        assert len(per_locus) == svcomposite_merging.DEGENERATE_POPULATION_WARN_LIMIT
        assert any("suppressed" in r.message for r in caplog.records)

    def test_a_degenerate_pair_inside_the_tolerance_never_reaches_the_effect_size(
        self, caplog
    ):
        """The short-circuit above the call: no warning, and the arm accepts.

        With a homopolymer sequence the complexity allowance is the full event
        size, so the two population means are already within it and arm 2 never
        asks for an effect size. Degeneracy here is not an error -- nothing was
        computed on it.
        """
        a = _make_insertion_composite(
            size=100,
            size_distortions={"r1": 0.0},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=180,
            size_distortions={"r1": 0.0},
            samplename="sample2",
            consensusID="2.0",
        )
        with caplog.at_level(
            "WARNING", logger="svirlpool.svcalling.svcomposite_merging"
        ):
            _, _, population_similar, cohensD = svcomposite_merging._similar_size(
                a, b, 0.1, 1.0, 2.0
            )
        assert population_similar
        assert isinstance(cohensD, float) and np.isnan(cohensD), (
            "nan still means 'not computed', and must stay distinguishable from "
            "'computed and not estimable'"
        )
        assert not [
            r for r in caplog.records if "DEGENERATE_SIZE_POPULATION" in r.message
        ]

    def test_empty_populations_are_not_reported_as_degenerate(self, caplog):
        """No population at all is missing data, not a degenerate one."""
        a, b = self._pair(100, 200, None, None)
        with caplog.at_level(
            "WARNING", logger="svirlpool.svcalling.svcomposite_merging"
        ):
            _, _, population_similar, cohensD = svcomposite_merging._similar_size(
                a, b, 0.1, 0.0, 2.0
            )
        assert not population_similar
        assert isinstance(cohensD, float) and np.isnan(cohensD)
        assert not [
            r for r in caplog.records if "DEGENERATE_SIZE_POPULATION" in r.message
        ]

    def test_a_real_population_still_yields_a_real_effect_size(self, caplog):
        """Regression guard: the arm is guarded, not disabled.

        Populations with genuine spread must still be measured, and must still
        be able to carry a merge that the fractional arm rejects.
        """
        spread = {
            f"r{i}": v
            for i, v in enumerate([-80.0, -50.0, -20.0, 0.0, 20.0, 50.0, 80.0])
        }
        a, b = self._pair(500, 600, spread, spread)
        with caplog.at_level(
            "WARNING", logger="svirlpool.svcalling.svcomposite_merging"
        ):
            similar, fraction_similar, population_similar, cohensD = (
                svcomposite_merging._similar_size(a, b, 0.1, 0.0, 2.0)
            )
        assert not fraction_similar, "the fractional arm must not be what passes here"
        assert population_similar
        assert similar
        assert cohensD is not None and np.isfinite(cohensD)
        assert not [
            r for r in caplog.records if "DEGENERATE_SIZE_POPULATION" in r.message
        ]

    def test_a_real_population_still_rejects_when_the_effect_size_is_large(self):
        """The same populations, a size gap they cannot absorb."""
        spread = {
            f"r{i}": v for i, v in enumerate([-8.0, -5.0, -2.0, 0.0, 2.0, 5.0, 8.0])
        }
        a, b = self._pair(500, 900, spread, spread)
        _, _, population_similar, cohensD = svcomposite_merging._similar_size(
            a, b, 0.1, 0.0, 2.0
        )
        assert cohensD is not None and np.isfinite(cohensD)
        assert abs(cohensD) > 2.0
        assert not population_similar

    def test_max_cohens_d_still_moves_the_decision(self):
        """`--max_cohens_d` must actually be a knob again.

        On `main` it had no influence on any merge decision at all: with
        degenerate populations the comparison was `inf <= d`, False for every
        finite `d`.
        """
        spread = {
            f"r{i}": v
            for i, v in enumerate([-40.0, -25.0, -10.0, 0.0, 10.0, 25.0, 40.0])
        }
        a, b = self._pair(500, 600, spread, spread)
        _, _, lenient, _ = svcomposite_merging._similar_size(a, b, 0.1, 0.0, 5.0)
        _, _, strict, _ = svcomposite_merging._similar_size(a, b, 0.1, 0.0, 0.5)
        assert lenient
        assert not strict


class TestVerboseReportsWhyThereIsNoEffectSize:
    """The three `verbose` consumers of the returned Cohen's D.

    They were the only place the value surfaced at all, and they used `nan` as
    "not computed". They now have to distinguish that from "computed, and not
    estimable", which is the whole point of F2.
    """

    @pytest.fixture(autouse=True)
    def _reset_rate_limiter(self):
        svcomposite_merging._reset_degenerate_population_warnings()
        yield

    @staticmethod
    def _merge(kind, size_a, size_b, distortions, *, scale, capsys):
        a = _make_composite(
            kind,
            size=size_a,
            size_distortions=distortions,
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(size_a, seed=21),
        )
        b = _make_composite(
            kind,
            size=size_b,
            size_distortions=distortions,
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(size_b, seed=22),
        )
        _MERGERS[kind](
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=10 * max(size_a, size_b),
            min_kmer_overlap=0.0,
            scale_by_complexity_factor=scale,
            verbose=True,
        )
        return capsys.readouterr().out

    def test_degenerate_case_is_named_in_the_report(self, capsys):
        for kind in _KINDS:
            out = self._merge(
                kind, 100, 200, {"r1": 0.0, "r2": 0.0}, scale=0.0, capsys=capsys
            )
            assert "Cohen's D" in out, kind
            assert "not estimable" in out, kind

    def test_computed_case_still_prints_the_number(self, capsys):
        spread = {
            f"r{i}": v
            for i, v in enumerate([-80.0, -50.0, -20.0, 0.0, 20.0, 50.0, 80.0])
        }
        for kind in _KINDS:
            out = self._merge(kind, 500, 600, spread, scale=0.0, capsys=capsys)
            assert "Cohen's D:" in out, kind
            assert "not estimable" not in out, kind
            assert "threshold: 2.0" in out, kind

    def test_not_computed_case_says_so(self, capsys):
        """Populations absent: nothing was computed, and the report says that
        rather than printing nothing at all."""
        for kind in _KINDS:
            out = self._merge(kind, 500, 520, None, scale=0.0, capsys=capsys)
            assert "not computed" in out, kind


class TestSizePopulationsKeepSubBasepairResolution:
    """N2 -- the distortion values must reach Cohen's *d* unquantised.

    `get_size_populations()` used to wrap every value in `int()` and
    `_similar_size` cast the result to `np.int32` on top. Both were invisible
    while F1 pinned every value at exactly 0.0. With real distortion estimates
    they are not: truncation toward zero shrinks the within-group spread, which
    inflates `|d|` and biases the arm toward rejection -- and, at the values the
    noise model actually produces, can collapse a perfectly informative
    population into a *constant* one, manufacturing the very degeneracy this
    section is about.
    """

    def test_populations_are_not_truncated(self):
        composite = _make_insertion_composite(
            size=500, size_distortions={"r1": 35.95, "r2": 129.96, "r3": 16.46}
        )
        assert composite.get_size_populations() == pytest.approx([35.95, 129.96, 16.46])

    def test_sub_basepair_distortions_do_not_become_a_degenerate_population(self):
        """Three distinct sub-bp estimates truncate to three identical zeros.

        Both populations below have a real spread of 0.8 bp. Truncated, both
        become constant vectors and the pair reached `cohens_d` as a degenerate
        one -- `inf` pre-fix, "not estimable" under the new contract. Untruncated
        it is an ordinary, perfectly measurable effect size.
        """
        a = _make_insertion_composite(
            size=100,
            size_distortions={"r1": 0.1, "r2": 0.5, "r3": 0.9},
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(100, seed=31),
        )
        b = _make_insertion_composite(
            size=102,
            size_distortions={"r1": -0.9, "r2": -0.5, "r3": -0.1},
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(102, seed=32),
        )
        _, _, _, cohensD = svcomposite_merging._similar_size(a, b, 0.0, 0.0, 2.0)
        assert cohensD is not None, (
            "truncation turned two spread populations into constant vectors"
        )
        assert np.isfinite(cohensD)
