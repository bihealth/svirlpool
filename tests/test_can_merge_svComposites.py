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
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.SVcomposite import SVcomposite
from svirlpool.svcalling.svcomposite_merging import (
    can_merge_svComposites_deletions,
    can_merge_svComposites_insertions,
    can_merge_svComposites_inversions,
    sizetolerance_from_SVcomposite,
)

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
        """Documents that the (1+tol)*max_size fraction formula always passes for positive sizes.

        Even with very different sizes (100 vs 200, 100% difference), the check
        |a_adj - b_adj| <= (1 + tol) * max(a_adj, b_adj) is trivially True because
        the left side equals max - min <= max <= (1+tol)*max. Merging is only blocked
        by the k-mer or proximity checks, not the fraction formula.
        """
        a = _make_insertion_composite(
            size=100,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=200,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
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
        """Documents that empty populations with differing sizes still merge via the fraction formula.

        With no size_distortions, population_similar=False. However fraction_similar is still
        True (the formula is trivially satisfied), so similar_size=True and merge proceeds.
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
        assert can_merge_svComposites_insertions(
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
        """With the complexity-adjusted size formula, 5% size differences merge regardless of fraction tolerance."""
        # 500 vs 525 → diff/max = 25/525 ≈ 4.8%
        # The (1 + tolerance) * max_size formula accepts any positive-size pair;
        # both assertions merge because similar_size is True via the fraction check.
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
        # With 10%: merges (fraction passes)
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
            scale_by_complexity_factor=0.0,
        )
        # With 1%: also merges because the adjusted size formula still passes
        assert can_merge_svComposites_insertions(
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
        """Documents that the (1+tol)*max_size fraction formula always passes for positive sizes.

        Even with 150% size difference (100 vs 250), the check is trivially True. Population
        check fails (high Cohen's D), but fraction_similar=True wins via 'or', so merge proceeds.
        """
        a = _make_deletion_composite(
            size=100,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=250,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=350,
            min_kmer_overlap=0.0,
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
        """Documents that empty populations with differing sizes still merge via the fraction formula.

        With no size_distortions, population_similar=False. However fraction_similar is still
        True (the formula is trivially satisfied), so similar_size=True and merge proceeds.
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
        assert can_merge_svComposites_deletions(
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
        """With a very strict Cohen's D threshold and borderline sizes, behavior changes."""
        # 500 vs 560 → diff/max = 60/560 ≈ 10.7% → fails 10% fraction
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
        # With lenient d=2.0: Cohen's D ~1.1 → merge via population test
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )
        # With strict d=0.5 and strict fraction: fraction check still passes for large sequences
        # (the (1 + tolerance) * max_size formula accepts any positive-size pair)
        assert can_merge_svComposites_deletions(
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


class TestSizeGateTautologyCharacterisation:
    """Characterisation of the size gate as it behaves on `main` (F3).

    On `main` the fractional arm reads

        |a_adj - b_adj| <= (1.0 + tol) * max(a_adj, b_adj)

    and `get_size()` returns a magnitude, so `|x - y| <= max(x, y)` already
    holds for every non-negative pair -- before the `1.0 +` even adds its bonus.
    The arm is therefore true whenever `max_size > 0`, for every tolerance,
    which makes the whole size gate inert and the Cohen's-D arm it is OR-ed with
    dead code.

    These tests PASS on `main` and MUST FAIL once the `1.0 +` is dropped. They
    exist so that the defect is pinned by an executable statement rather than by
    a prose description, and so that the fix has something to break.
    """

    def test_fraction_arm_accepts_every_size_pair(self):
        """Even a 50 bp and a 5000 bp event merge, at any tolerance."""
        failures = []
        for kind in _KINDS:
            for size_a, size_b in _SIZE_PAIRS:
                for tolerance in (0.0, 0.01, 0.06, 0.1):
                    if not _size_gate(kind, size_a, size_b, tolerance=tolerance):
                        failures.append((kind, size_a, size_b, tolerance))
        assert failures == [], (
            f"the gate rejected {len(failures)} pair(s) it accepts on main: {failures}"
        )


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


@pytest.mark.xfail(
    strict=True,
    reason="F3: on main the fractional arm is |a-b| <= (1+tol)*max(a,b), which is "
    "tautologically true for non-negative sizes. Passes once the '1.0 +' is dropped.",
)
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

    @pytest.mark.xfail(
        strict=True,
        reason="F4: on main the function returns the bare mean complexity, so the "
        "relation is inverted (high complexity returns the larger number).",
    )
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

    @pytest.mark.xfail(
        strict=True,
        reason="F4: on main an all-zero track and an absent track are conflated, and "
        "both return 0.0; an absent track additionally raises TypeError for "
        "insertions and inversions.",
    )
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

    @pytest.mark.xfail(
        strict=True,
        reason="F3/F4: on main the gate is inert, so both regimes merge.",
    )
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

    @pytest.mark.xfail(
        strict=True,
        reason="F3/F4: on main the gate is inert, so the allowance cannot be withdrawn.",
    )
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
