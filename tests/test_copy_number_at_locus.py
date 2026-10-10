"""The copy number of an insertion's locus comes from the CN track too.

An insertion's reference region is a point (start == end). ``IntervalTree[s:s]``
is empty, so ``copy_number_at_locus`` used to fall back to the diploid default
for every insertion: on a male chrX (CN 1) deletions were written haploid and
insertions diploid (2026-10-10, HG002 30x: 383 INS "1/1" against 165 DEL "1").
"""

from intervaltree import IntervalTree
from xxhash import xxh64

from svirlpool.svcalling.multisample_sv_calling import (
    copy_number_at_locus,
    genotype_of_sample,
)


def _cn_tracks(cn: int) -> dict[str, dict[str, IntervalTree]]:
    tree = IntervalTree()
    tree[0:1_000_000] = cn
    return {"HG002": {"chrX": tree}}


def test_point_locus_reads_the_track():
    assert copy_number_at_locus("HG002", "chrX", 500_000, 500_000, _cn_tracks(1)) == 1


def test_span_locus_is_unchanged():
    assert copy_number_at_locus("HG002", "chrX", 500_000, 500_200, _cn_tracks(1)) == 1


def test_point_locus_outside_the_track_stays_diploid():
    assert copy_number_at_locus("HG002", "chrX", 2_000_000, 2_000_000, _cn_tracks(1)) == 2


def test_insertion_on_a_haploid_track_gets_a_haploid_genotype():
    reads = [xxh64(f"read{i}").intdigest() for i in range(20)]
    cov = IntervalTree()
    for r in reads:
        cov[499_000:501_000] = r
    gt = genotype_of_sample(
        samplename="HG002",
        chrname="chrX",
        start=500_000,
        end=500_000,
        raw_alt_reads=set(reads),
        covtrees={"HG002": {"chrX": cov}},
        cn_tracks=_cn_tracks(1),
        apply_breakpoint_margin=True,
    )
    assert gt.genotype == "1"
