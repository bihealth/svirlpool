"""Unit tests for can_merge_svComposites_insertions and can_merge_svComposites_deletions.

Tests the two-test size comparison logic:
  1) Fractional size difference check (default 10%)
  2) Population-driven Cohen's D on background signals
If either test passes, the variants are considered similar in size.
"""

import random

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.SVcomposite import SVcomposite
from svirlpool.svcalling.svcomposite_merging import (
    can_merge_svComposites_deletions,
    can_merge_svComposites_insertions,
    sizetolerance_from_SVcomposite,
)

# ---------------------------------------------------------------------------
# Helper factories
# ---------------------------------------------------------------------------


def _random_dna(length: int, seed: int) -> str:
    """A high-complexity DNA string, reproducible from `seed`.

    The merge criterion is parameterised by sequence complexity, so a test that
    leaves the sequence at the default homopolymer is not testing the threshold
    it names — it is testing the maximum complexity allowance. Use this wherever
    the intent is "well-resolved sequence"; use an explicit "A" * n wherever the
    intent is a repeat.

    Note complexity is only computed for sequences up to
    `sequence_complexity_max_length` (300 bp, SVpatterns.set_sequence). Above
    that no track exists and the complexity allowance is zero regardless of
    sequence, so keep sizes at or below 300 in tests that exercise complexity.
    """
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


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


# ===========================================================================
# INSERTION TESTS
# ===========================================================================


class TestComplexityTolerance:
    """The complexity -> tolerance relation the size gate is built on.

    Kept separate and first, so that if `complexity_local_track` ever changes
    scale or normalisation the failure lands here rather than surfacing as an
    inexplicable merge/reject flip in the tests below.
    """

    def test_low_complexity_grants_more_tolerance_than_high(self):
        """Lower complexity must grant a larger size tolerance. Monotone, in bp.

        High complexity means well-resolved sequence: the aligner has little
        freedom in where it places an indel or how large it calls it, so two
        sizes should have to agree closely. Low complexity means the opposite.
        """
        homopolymer = _make_insertion_composite(size=200, sequence="A" * 200)
        random_dna = _make_insertion_composite(
            size=200, sequence=_random_dna(200, seed=42)
        )

        tol_low = sizetolerance_from_SVcomposite(homopolymer)
        tol_high = sizetolerance_from_SVcomposite(random_dna)

        assert tol_low > tol_high, (
            f"homopolymer tolerance {tol_low:.1f} bp should exceed "
            f"random-sequence tolerance {tol_high:.1f} bp"
        )
        # Bounded by the event's own size, and non-negative.
        assert 0.0 <= tol_high <= 200.0
        assert 0.0 <= tol_low <= 200.0

    def test_absent_complexity_track_grants_no_tolerance(self):
        """No complexity evidence must fail CLOSED, not open.

        `set_sequence` only computes a complexity track up to
        `sequence_complexity_max_length` (300 bp), so a larger insertion has no
        track at all. That is missing data, and it must grant no allowance — as
        opposed to a track that is present and all zeros, which is genuine
        minimum complexity and does grant the maximum (covered above).
        """
        no_track = _make_insertion_composite(size=1000, sequence="A" * 1000)
        assert sizetolerance_from_SVcomposite(no_track) == 0.0


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
        kwargs = dict(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=300,
            min_kmer_overlap=0.0,
        )
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
        """With the corrected fraction formula, a 150% size difference is rejected.

        Sizes 100 vs 250 fail the corrected fraction test (the relative size
        difference far exceeds the 10% tolerance). The population check also
        fails (high Cohen's D), so the merge is correctly rejected.
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
        assert not can_merge_svComposites_deletions(
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
