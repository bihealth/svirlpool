"""Unit tests for the vertical-merge size-similarity gate.

The gate is one floored relative bound on the two sizes,
|a - b| <= max(floor, tol * ref(a, b)) (see `svcomposite_merging.size_tolerance`).
A Cohen's d arm on the per-read size-distortion populations used to be OR-ed
with it; it was removed, and sizes are all that decide.

The same gate is reached from three entry points -- insertions, deletions and
inversions -- and the tests below deliberately exercise all three with the same
inputs, because the gate was once triplicated by copy-paste and drifted.
"""

import json
import random
from gzip import open as gzopen
from pathlib import Path

import cattrs
import numpy as np
import pytest

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.signalprocessing.alignments_to_rafs import (
    get_start_end,
    parse_SVsignals_from_alignment,
)
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.SVcomposite import SVcomposite
from svirlpool.svcalling.svcomposite_merging import (
    can_merge_svComposites_deletions,
    can_merge_svComposites_insertions,
    can_merge_svComposites_inversions,
)
from svirlpool.svcalling.svcomposite_utils import cohens_d
from svirlpool.util.datatypes import Alignment

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
    chr: str = "chr1",
    ref_start: int = 1000,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    sequence: str | None = None,
    reads: list[str] | None = None,
) -> SVcomposite:
    """Create an insertion SVcomposite with controllable size."""
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
    )
    if sequence is None:
        sequence = "A" * size
    pattern.set_sequence(sequence)
    return SVcomposite.from_SVpattern(pattern)


def _make_deletion_composite(
    *,
    size: int = 500,
    chr: str = "chr1",
    ref_start: int = 1000,
    samplename: str = "sample1",
    consensusID: str = "1.0",
    sequence: str | None = None,
    reads: list[str] | None = None,
) -> SVcomposite:
    """Create a deletion SVcomposite with controllable size."""
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
    sequence_a: str | None = None,
    sequence_b: str | None = None,
) -> bool:
    """Merge decision with every gate other than the size gate neutralised.

    `near` is set far beyond both events and `min_kmer_overlap` to 0.0, so the
    proximity and k-mer arms always pass and the return value *is* the size
    gate's verdict.
    """
    a = _make_composite(
        kind,
        size=size_a,
        samplename="sample1",
        consensusID="1.0",
        sequence=sequence_a,
    )
    b = _make_composite(
        kind,
        size=size_b,
        samplename="sample2",
        consensusID="2.0",
        sequence=sequence_b,
    )
    return _MERGERS[kind](
        a=a,
        b=b,
        apriori_size_difference_fraction_tolerance=tolerance,
        near=10 * max(size_a, size_b, 100),
        min_kmer_overlap=0.0,
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
        """Two identical-size insertions at the same locus should merge."""
        a = _make_insertion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=500,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_similar_size_within_fraction_tolerance(self):
        """Two insertions within 10% size difference should merge via fractional test."""
        # 500 vs 540 → diff=40, 10% of 540=54 → within tolerance
        a = _make_insertion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=540,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_different_sizes_beyond_fraction_reject(self):
        """Sizes differ by >10%: rejected.

        A Cohen's d arm on per-read size-distortion populations used to merge
        such pairs when the populations overlapped. It was removed; only the
        sizes decide.
        """
        # sizes: 500 vs 600 → diff/max = 100/600 ≈ 16.7% → fraction test fails
        a = _make_insertion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=600,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_very_different_sizes_reject(self):
        """Very different sizes are rejected when the sequence is well resolved.

        Sizes 100 vs 200 (a 2x difference) fail the fraction test, which requires
        the relative size difference to be within `tol` (here 10%), and nothing
        else can accept a pair.
        """
        a = _make_insertion_composite(
            size=100,
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(100, seed=1),
        )
        b = _make_insertion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(200, seed=2),
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=300,
            min_kmer_overlap=0.0,
        )

    def test_very_different_sizes_do_not_merge_in_low_complexity(self):
        """A homopolymer grants no size allowance any more.

        The complexity gate accepted a 100 bp and a 200 bp insertion inside a
        poly-A run as one allele family. It was removed: it never improved the
        trio benchmark, and in a trio differently sized repeat alleles are
        mostly distinct alleles. Sequence complexity no longer decides a merge.
        """
        a = _make_insertion_composite(
            size=100,
            samplename="sample1",
            consensusID="1.0",
            sequence="A" * 100,
        )
        b = _make_insertion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
            sequence="A" * 200,
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=300,
            min_kmer_overlap=0.0,
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
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(200, seed=3),
        )
        b = _make_insertion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
            sequence="A" * 200,
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=300,
            min_kmer_overlap=0.0,
        )

    def test_small_size_difference_passes_the_fraction_test(self):
        """500 vs 520 bp is within the 10% fraction tolerance and merges."""
        a = _make_insertion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=520,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_large_size_difference_fails_the_fraction_test(self):
        """A large size difference fails the fraction test and is rejected.

        The fraction test fails because a 100 vs 200 size difference (100%)
        exceeds the 10% tolerance, so similar_size=False and the merge is correctly
        rejected.
        """
        a = _make_insertion_composite(
            size=100,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=300,
            min_kmer_overlap=0.0,
        )

    def test_not_near_rejects(self):
        """Insertions on different chromosomes should not merge."""
        a = _make_insertion_composite(
            size=500,
            chr="chr1",
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=500,
            chr="chr2",
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_kmer_similarity_rejects(self):
        """Same-size insertions with different sequence content should be rejected by k-mer check."""
        a = _make_insertion_composite(
            size=200,
            samplename="sample1",
            consensusID="1.0",
            sequence="ATCGATCG" * 25,  # 200bp AT-rich
        )
        b = _make_insertion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
            sequence="GCGCGCGC" * 25,  # 200bp GC-rich
        )
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_strict_tolerance_rejects_borderline(self):
        """A borderline 4.8% size difference merges under a lenient tolerance but is rejected under a strict one."""
        # 500 vs 525 → diff/max = 25/525 ≈ 4.8%
        # The corrected check requires diff <= tol * max_size, so whether this pair
        # merges now genuinely depends on the tolerance value, as the test name implies.
        a = _make_insertion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_insertion_composite(
            size=525,
            samplename="sample2",
            consensusID="2.0",
        )
        # With 10% tolerance: 4.8% difference is within tolerance -> merges
        assert can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )
        # With 1% tolerance: 4.8% difference exceeds tolerance -> rejected
        assert not can_merge_svComposites_insertions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.01,
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
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=500,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_similar_size_within_fraction_tolerance(self):
        """Deletions within 10% size difference should merge via fractional test."""
        # 1000 vs 1080 → diff=80, 10% of 1080=108 → within tolerance
        a = _make_deletion_composite(
            size=1000,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=1080,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_different_sizes_beyond_fraction_reject(self):
        """Sizes differ by >10%: rejected."""
        # 500 vs 600 → 16.7% → fraction fails
        a = _make_deletion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=600,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_very_different_sizes_reject(self):
        """A 150% size difference is rejected.

        Sizes 100 vs 250 fail the fraction test (the relative size difference far
        exceeds the 10% tolerance), and nothing else can accept a pair.
        """
        a = _make_deletion_composite(
            size=100,
            samplename="sample1",
            consensusID="1.0",
            sequence=_random_dna(100, seed=11),
        )
        b = _make_deletion_composite(
            size=250,
            samplename="sample2",
            consensusID="2.0",
            sequence=_random_dna(250, seed=12),
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=350,
            min_kmer_overlap=0.0,
        )

    def test_very_different_sizes_do_not_merge_in_low_complexity(self):
        """Deletion counterpart: a homopolymer grants no size allowance."""
        a = _make_deletion_composite(
            size=100,
            samplename="sample1",
            consensusID="1.0",
            sequence="A" * 100,
        )
        b = _make_deletion_composite(
            size=250,
            samplename="sample2",
            consensusID="2.0",
            sequence="A" * 250,
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=350,
            min_kmer_overlap=0.0,
        )

    def test_small_size_difference_passes_the_fraction_test(self):
        """500 vs 530 bp is within the 10% fraction tolerance and merges."""
        a = _make_deletion_composite(
            size=500,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=530,
            samplename="sample2",
            consensusID="2.0",
        )
        assert can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_large_size_difference_fails_the_fraction_test(self):
        """A large size difference fails the fraction test and is rejected.

        The fraction test fails because a 100 vs 200 size difference (100%)
        exceeds the 10% tolerance, so similar_size=False and the merge is correctly
        rejected.
        """
        a = _make_deletion_composite(
            size=100,
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=300,
            min_kmer_overlap=0.0,
        )

    def test_not_near_rejects(self):
        """Deletions on different chromosomes should not merge."""
        a = _make_deletion_composite(
            size=500,
            chr="chr1",
            samplename="sample1",
            consensusID="1.0",
        )
        b = _make_deletion_composite(
            size=500,
            chr="chr2",
            samplename="sample2",
            consensusID="2.0",
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
        )

    def test_kmer_similarity_rejects_deletions(self):
        """Same-size deletions with different reference sequence content should be rejected."""
        a = _make_deletion_composite(
            size=200,
            samplename="sample1",
            consensusID="1.0",
            sequence="ATCGATCG" * 25,  # 200bp AT-rich
        )
        b = _make_deletion_composite(
            size=200,
            samplename="sample2",
            consensusID="2.0",
            sequence="GCGCGCGC" * 25,  # 200bp GC-rich
        )
        assert not can_merge_svComposites_deletions(
            a=a,
            b=b,
            apriori_size_difference_fraction_tolerance=0.1,
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

    def test_all_three_paths_agree_in_low_and_high_complexity(self):
        """Same, in low- and high-complexity sequence."""
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
                            sequence_a=seq_a,
                            sequence_b=seq_b,
                        )
                        for kind in _KINDS
                    }
                    if len(set(verdicts.values())) != 1:
                        disagreements.append(
                            (
                                size_a,
                                size_b,
                                tolerance,
                                sequence,
                                verdicts,
                            )
                        )
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


class TestSequenceComplexityDecidesNothing:
    """The complexity size allowance was removed from the vertical merge.

    It granted `(1 - mean_complexity) * |size|` bp and short-circuited the
    size test. It is gone: the size test sees sizes only, so the same sizes decide the same way whatever the sequence.
    """

    def test_low_and_high_complexity_decide_alike(self):
        for size_a, size_b in [(100, 200), (200, 205), (100, 250), (300, 330)]:
            for kind in _KINDS:
                low = _size_gate(
                    kind,
                    size_a,
                    size_b,
                    tolerance=0.1,
                    sequence_a="A" * size_a,
                    sequence_b="A" * size_b,
                )
                high = _size_gate(
                    kind,
                    size_a,
                    size_b,
                    tolerance=0.1,
                    sequence_a=_random_dna(size_a, seed=1),
                    sequence_b=_random_dna(size_b, seed=2),
                )
                assert low == high, (kind, size_a, size_b)

    def test_a_complexity_difference_is_not_a_size_difference(self):
        """Regression guard: two events of the SAME size must always merge."""
        for kind in _KINDS:
            assert _size_gate(
                kind,
                200,
                200,
                tolerance=0.0,
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

    def test_a_constant_pair_is_detected_when_the_pooled_std_is_not_exactly_zero(
        self,
    ):
        """The degeneracy test cannot be `pooled_std == 0`, and this is why.

        `np.mean` of an array of *identical* elements is not exactly that
        element: numpy sums pairwise, and the value need not be representable.
        19 copies of 6164.5678 have a range of exactly 0 but a sample standard
        deviation of 9.3e-13, so the old guard did not fire and the quotient
        came back around 1e16 -- finite, and indistinguishable in a log from a
        real effect size.

        This is not a contrived array. It is what arm 2 hands to `cohens_d` at
        every locus where the distortion estimates are all 0.0: a constant
        population, shifted by a fractional complexity tolerance. On the real
        chr6 VNTR pair (4,930 bp vs 18,028 bp) the value observed was
        -2.07e16, never `inf`.
        """
        x = np.full(19, 4930.0) + 1234.5678
        y = np.full(19, 18028.0) - 1234.5678
        assert np.ptp(x) == 0.0 and np.ptp(y) == 0.0, "both samples are constant"
        assert np.std(x, ddof=1) != 0.0, (
            "the premise of this test: the computed spread is not exactly zero"
        )
        assert cohens_d(x, y) is None

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


# ===========================================================================
# F5: THE log_size DISJUNCT
# ===========================================================================


def _fraction_arm_with_log_disjunct(
    size_a: float, size_b: float, tolerance: float
) -> bool:
    """Arm 1 exactly as it stood before F5, disjunct included.

    Kept verbatim as the reference the post-F5 form is compared against, so the
    equivalence below is asserted against the removed code and not merely
    against a restatement of what replaced it.
    """
    max_size = max(abs(size_a), abs(size_b))
    log_size = np.log2(abs(size_a - size_b) + 1)
    return bool(
        (max_size > 0 and abs(size_a - size_b) <= tolerance * max_size)
        or abs(size_a - size_b) < log_size
    )


def _fraction_arm_without_log_disjunct(
    size_a: float, size_b: float, tolerance: float
) -> bool:
    """Arm 1 as it stands after F5."""
    max_size = max(abs(size_a), abs(size_b))
    return bool(max_size > 0 and abs(size_a - size_b) <= tolerance * max_size)


class TestLogSizeDisjunctIsUnreachable:
    """`|a - b| < log2(|a - b| + 1)` never fired, and could not have.

    The disjunct read as a floor on the admissible absolute size difference --
    a small-size safety net -- and was not one. F5 removes it. These tests are
    the proof that removing it is behaviour-preserving; there is deliberately
    no fails-before/passes-after test, because a pure dead-code deletion cannot
    have one, and manufacturing one would misrepresent the change.
    """

    def test_no_integer_difference_satisfies_it(self):
        """Exhaustive over d in [0, 100000), as the defect catalogue claims."""
        d = np.arange(0, 100_000, dtype=np.float64)
        satisfied = d < np.log2(d + 1.0)
        assert not satisfied.any(), (
            f"the disjunct fires for integer differences {d[satisfied][:10]}"
        )

    def test_over_the_reals_the_solution_set_is_exactly_the_open_unit_interval(self):
        """The catalogue's "no non-negative solution" holds over Z, not over R.

        Let f(d) = log2(1 + d) - d on d >= 0. Then

            f(0) = 0,  f(1) = log2(2) - 1 = 0,
            f'(d) = 1 / ((1 + d) * ln 2) - 1,

        so f'(d) = 0 at d = 1/ln2 - 1 ~= 0.4427 and f' > 0 to the left of it,
        f' < 0 to the right. f is therefore strictly concave with a single
        interior maximum between its two roots 0 and 1, positive on (0, 1) and
        strictly negative on (1, inf). The inequality d < log2(d + 1) is
        f(d) > 0, so its non-negative solution set is the *open* interval
        (0, 1): empty at both endpoints, and empty over the integers, but not
        empty over the reals.

        That distinction matters for how the disjunct is retired. It was not
        merely unreachable, it was unreachable *because* sizes happen to be
        integral (see the next test) -- so it would have woken up silently, and
        admitted sub-base-pair differences, had size ever become fractional.
        Deleting it removes that latent coupling; the catalogue's stated reason
        for deleting it does not survive contact with the real line.
        """

        def f(d: float) -> float:
            return float(np.log2(1.0 + d) - d)

        assert f(0.0) == 0.0
        assert f(1.0) == 0.0
        # Positive strictly between the two roots -- the disjunct *is*
        # satisfiable over the reals.
        for d in (1e-9, 0.01, 0.25, 1.0 / np.log(2.0) - 1.0, 0.5, 0.9, 1.0 - 1e-9):
            assert f(d) > 0.0, f"expected f({d}) > 0"
            assert d < np.log2(d + 1.0)
        # ... and strictly negative everywhere above 1.
        for d in (1.0 + 1e-9, 1.5, 2.0, 12.0, 1e3, 1e9):
            assert f(d) < 0.0, f"expected f({d}) < 0"
            assert not d < np.log2(d + 1.0)
        # The interior maximum sits where the derivative vanishes.
        peak = 1.0 / np.log(2.0) - 1.0
        grid = np.linspace(0.0, 1.0, 100_001)
        assert abs(grid[np.argmax(np.log2(1.0 + grid) - grid)] - peak) < 1e-4

    def test_svcomposite_sizes_are_integral_so_the_disjunct_was_dead_in_production(
        self,
    ):
        """`get_size()` is int-valued for all three SV types.

        Sizes are differences of alignment coordinates, so `|size_a - size_b|`
        is a non-negative integer and never lands in (0, 1).
        """
        for kind in _KINDS:
            for size in (1, 12, 57, 500, 18028):
                composite = _make_composite(kind, size=size)
                got = composite.get_size()
                assert isinstance(got, int), f"{kind}: get_size() returned {type(got)}"
                assert got == size


class TestRemovingTheLogSizeDisjunctChangesNoVerdict:
    """The F5 deletion is behaviour-preserving, on the arm and on the gate."""

    _TOLERANCES = (0.0, 0.001, 0.01, 0.06, 0.1, 0.25, 0.5, 0.9, 1.0)

    def test_the_two_arm_forms_agree_on_every_integer_size_pair(self):
        """Direct before/after comparison of arm 1, with and without the disjunct."""
        sizes = list(range(0, 60)) + [
            12,
            30,
            50,
            57,
            58,
            100,
            182,
            193,
            199,
            200,
            201,
            365,
            390,
            675,
            792,
            2554,
            2782,
            4930,
            18028,
        ]
        for tolerance in self._TOLERANCES:
            for size_a in sizes:
                for size_b in sizes:
                    with_disjunct = _fraction_arm_with_log_disjunct(
                        size_a, size_b, tolerance
                    )
                    without = _fraction_arm_without_log_disjunct(
                        size_a, size_b, tolerance
                    )
                    assert with_disjunct == without, (
                        f"verdict moved at tol={tolerance} sizes=({size_a}, {size_b}): "
                        f"{with_disjunct} -> {without}"
                    )

    def test_the_gate_matches_the_pure_fractional_bound_for_all_three_sv_types(self):
        """The live gate's verdict is the fractional bound, with no floor under it.

        The gate's return value *is* arm 1. Sweeping sizes, tolerances and all three entry
        points, it agrees with the disjunct-free bound everywhere -- and, by the
        test above, therefore with the pre-F5 form as well.
        """
        sizes = (1, 2, 11, 12, 13, 24, 25, 30, 50, 57, 58, 182, 193, 199, 200, 201, 390)
        for kind in _KINDS:
            for tolerance in self._TOLERANCES:
                for size_a in sizes:
                    for size_b in sizes:
                        expected = _fraction_arm_without_log_disjunct(
                            size_a, size_b, tolerance
                        )
                        got = _size_gate(kind, size_a, size_b, tolerance=tolerance)
                        assert got == expected, (
                            f"{kind}: tol={tolerance} sizes=({size_a}, {size_b}) "
                            f"gave {got}, expected {expected}"
                        )

    def test_the_gate_still_agrees_across_the_catalogue_size_pairs(self):
        """The pairs the other size-gate tests use, re-checked against the bound."""
        for kind in _KINDS:
            for size_a, size_b in _SIZE_PAIRS:
                for tolerance in self._TOLERANCES:
                    expected = _fraction_arm_without_log_disjunct(
                        size_a, size_b, tolerance
                    )
                    assert (
                        _size_gate(kind, size_a, size_b, tolerance=tolerance)
                        == expected
                    ), f"{kind}: ({size_a}, {size_b}) at tol={tolerance}"

    def test_no_floor_is_granted_below_the_fractional_bound(self):
        """By default there is no absolute-difference floor, in particular none at 12 bp.

        A floor option (--size-tolerance-floor) was tried and removed: floors of
        3-20 bp lowered precision in the svp_merging size study (2026-09-30).

        The dissertation text proposes F = 12 bp -- the consensus-level indel
        parse threshold (`consensus_align.py --min-signal-size`, default 12).
        F5 deliberately does *not* introduce it: a floor is a loosening of the
        gate and belongs in its own, separately benchmarkable change, not behind
        a dead-code removal. This test pins the absence, so that adding a floor
        later is a visible, deliberate edit rather than a silent one.

        The window in which a 12 bp floor would differ from the default
        fractional bound is max(size_a, size_b) < 12 / 0.06 = 200 bp; above that
        the fractional arm already admits +/-12 bp on its own.
        """
        tolerance = 0.06
        for kind in _KINDS:
            # Inside the window: a 12 bp difference on a 50 bp event is 24% and
            # is rejected, floor or no floor.
            assert not _size_gate(kind, 50, 62, tolerance=tolerance)
            assert not _size_gate(kind, 100, 112, tolerance=tolerance)
            # The tightest miss: 12 bp apart, just under the crossover.
            assert not _size_gate(kind, 185, 197, tolerance=tolerance)
            # At and above the crossover (0.06 * max >= 12, i.e. max >= 200)
            # the fractional arm admits 12 bp on its own.
            assert _size_gate(kind, 200, 212, tolerance=tolerance)
            assert _size_gate(kind, 1000, 1012, tolerance=tolerance)


# ===========================================================================
# INDEL SIGNAL SIGN CONVENTION (F8)
# ===========================================================================


def _load_simulated_alignments(name: str) -> list[Alignment]:
    """Load the simulated single-indel alignments used as fixtures."""
    path = Path(__file__).parent / "data" / "consensus_class" / name
    with gzopen(path, "rt") as f:
        data = json.load(f)
    return cattrs.structure(data["alignments"], list[Alignment])


class TestIndelSignalSignConvention:
    """Pin the sign convention of parsed indel signals.

    Deletion sizes are magnitudes from the moment they are parsed
    (``size=int(abs(delr - dell))`` in
    ``alignments_to_rafs.parse_SVsignals_from_alignment``); only
    ``SVsignal.sv_type`` tells insertions and deletions apart.
    """

    def test_indel_signals_are_unsigned_magnitudes_at_the_source(self):
        """A 20 bp deletion and a 20 bp insertion both parse to ``size == +20``."""
        sizes_by_type: dict[int, list[int]] = {}
        for name in (
            "simulated.with_deletion.json.gz",
            "simulated.with_insertion.json.gz",
        ):
            for aln in _load_simulated_alignments(name):
                pysam_aln = aln.to_pysam()
                ref_start, ref_end, read_start, read_end = get_start_end(pysam_aln)
                for signal in parse_SVsignals_from_alignment(
                    alignment=pysam_aln,
                    ref_start=ref_start,
                    ref_end=ref_end,
                    read_start=read_start,
                    read_end=read_end,
                    min_signal_size=10,
                    min_bnd_size=50,
                ):
                    sizes_by_type.setdefault(signal.sv_type, []).append(signal.size)

        assert 0 in sizes_by_type and 1 in sizes_by_type, sizes_by_type
        # sv_type 1 == deletion: the size is +20, NOT -20
        assert [20] == sizes_by_type[1]
        # sv_type 0 == insertion
        assert [20] == sizes_by_type[0]
        # only the separate sv_type field distinguishes them
        assert sizes_by_type[0] == sizes_by_type[1]
