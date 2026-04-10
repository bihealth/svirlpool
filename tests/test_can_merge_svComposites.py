"""Unit tests for can_merge_svComposites_insertions and can_merge_svComposites_deletions.

Tests the two-test size comparison logic:
  1) Fractional size difference check (default 10%)
  2) Population-driven Cohen's D on background signals
If either test passes, the variants are considered similar in size.
"""

from svirlpool.localassembly import SVprimitives, SVpatterns
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.SVcomposite import SVcomposite
from svirlpool.svcalling.multisample_sv_calling import (
    can_merge_svComposites_deletions,
    can_merge_svComposites_insertions,
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
        )

    def test_very_different_sizes_reject(self):
        """Insertions with very different sizes (>10% and non-overlapping populations) should not merge."""
        # 300 vs 600 → diff/max = 300/600 = 50% → fraction fails
        # tight populations around distinct means → Cohen's D very large
        a = _make_insertion_composite(
            size=300,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=600,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=600,
            min_kmer_overlap=0.0,
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
        )

    def test_empty_populations_fraction_fail(self):
        """When populations are empty and sizes differ >10%, should reject."""
        a = _make_insertion_composite(
            size=300,
            size_distortions=None,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=600,
            size_distortions=None,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=600,
            min_kmer_overlap=0.0,
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
        )

    def test_strict_tolerance_rejects_borderline(self):
        """With a very strict tolerance (1%), a 5% size difference should fail via fraction test."""
        # 500 vs 525 → diff/max = 25/525 ≈ 4.8% → fails 1% but passes 10%
        # tight populations → large Cohen's D
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
        # With default 10%: should merge (fraction passes)
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
        )
        # With 1%: fraction fails; tight populations → Cohen's D >2 → reject
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.01,
            d=2.0,
            near=150,
            min_kmer_overlap=0.7,
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
        """Deletions with very different sizes should not merge."""
        # 200 vs 800 → diff/max = 75% → fails fraction; tight populations → large Cohen's D
        a = _make_deletion_composite(
            size=200,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=800,
            size_distortions={"r1": 1, "r2": -1, "r3": 2},
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=1000,
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
        """Empty populations and different sizes → reject."""
        a = _make_deletion_composite(
            size=300,
            size_distortions=None,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=600,
            size_distortions=None,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            d=2.0,
            near=600,
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
        # With strict d=0.5 and strict fraction: both tests fail → reject
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.05,
            d=0.5,
            near=150,
            min_kmer_overlap=0.7,
        )
