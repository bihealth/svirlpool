"""--core-alignment-merge: same allele judged by aligning the consensus cores.

Synthetic haplotypes built from a random reference; the cores of the samples
cover different, overlapping stretches of their haplotype, as the cores of two
samples' candidate regions do. Positions are on the cores (core_start 0).
"""

import random
import shutil
import sys
from pathlib import Path

import pytest

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.SVcomposite import SVcomposite

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None
    and not (Path(sys.executable).parent / "minimap2").exists(),
    reason="minimap2 is needed",
)


def rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def with_snvs(seq: str, positions: list[int]) -> str:
    s = list(seq)
    for p in positions:
        s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
    return "".join(s)


def _primitive(sample, cid, sv_type, read_start, read_end, ref_start, ref_end, size, reverse):
    return SVprimitives.SVprimitive(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        size=size,
        sv_type=sv_type,
        chr="chr1",
        repeatIDs=[],
        original_alt_sequences=[],
        original_ref_sequences=[],
        samplename=sample,
        consensusID=cid,
        alignmentID=0,
        svID=0,
        aln_is_reverse=reverse,
        consensus_aln_interval=("chr1", 0, 10_000),
        genotypeMeasurement=genotyping.GenotypeMeasurement(
            start_on_consensus=read_start,
            supporting_reads_start=[f"{sample}_r{i}" for i in range(5)],
        ),
    )


def ins(sample, cid, core, start, size, ref_pos, reverse=False):
    """An insertion occupying core[start:start+size] (oriented core coordinates)."""
    seq = core[start : start + size]
    rs, re = start, start + size
    if reverse:  # the stored consensus is the reverse complement of `core`
        rs, re = len(core) - re, len(core) - rs
    pattern = SVpatterns.SVpatternInsertion(
        SVprimitives=[_primitive(sample, cid, 0, rs, re, ref_pos, ref_pos + 1, size, reverse)]
    )
    pattern.set_sequence(seq)
    return SVcomposite.from_SVpattern(pattern)


def dele(sample, cid, core, pos, deleted, ref_pos, reverse=False):
    """A deletion of `deleted` at the point `pos` of the (oriented) core."""
    rs = len(core) - pos if reverse else pos
    pattern = SVpatterns.SVpatternDeletion(
        SVprimitives=[
            _primitive(sample, cid, 1, rs, rs, ref_pos, ref_pos + len(deleted), len(deleted), reverse)
        ]
    )
    pattern.set_sequence(deleted)
    return SVcomposite.from_SVpattern(pattern)


def cores_of(**named):
    """{(sample, cid): (0, stored core)} from name -> (sample, cid, oriented core, reverse)."""
    out = {}
    for sample, cid, core, reverse in named.values():
        out[(sample, cid)] = (0, cam.revcomp(core) if reverse else core)
    return out


def never(a, b):
    return False


def always(a, b):
    return True


def run(group, cores, mode="align", params=None, old=never):
    stats = cam.GroupStats()
    merged, rest = cam.merge_group_with_core_alignments(
        group, cores, old, old, mode, params or cam.CoreAlnParams(), stats
    )
    assert rest == []
    return merged, stats


def records(merged):
    return sorted(
        sorted(f"{p.samplename}:{p.read_start}" for p in svc.svPatterns) for svc in merged
    )


# --------------------------------------------------------------------------- #


@pytest.fixture
def rng():
    return random.Random(7)


@pytest.mark.parametrize("mother_reverse", [False, True])
def test_one_haplotype_decomposed_differently_becomes_one_record(rng, mother_reverse):
    left, right, inserted = rand_seq(2000, rng), rand_seq(2000, rng), rand_seq(944, rng)
    hap = left + inserted + right
    child = hap[200 : 2000 + 944 + 1700]
    mother = with_snvs(hap[500 : 2000 + 944 + 1800], [100, 2700])  # other flank SNVs
    # child: INS 149 + 199 + 596 (the consensus-to-reference alignment split it)
    c0 = 2000 - 200
    group = [
        ins("HG002", "1.0", child, c0, 149, 2000),
        ins("HG002", "1.0", child, c0 + 149, 199, 2000),
        ins("HG002", "1.0", child, c0 + 348, 596, 2000),
        ins("HG004", "5.0", mother, 2000 - 500, 944, 2000, reverse=mother_reverse),
    ]
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, mother_reverse)
    )
    merged, stats = run(group, cores)
    assert len(merged) == 1
    assert len(merged[0].svPatterns) == 4
    assert stats.links_alignment == 3
    # the record's size is the sum of the pieces of its best-supported consensus
    assert merged[0].get_size() in (944, 944)


def _tandem_repeat(rng):
    """left flank, a 26-bp unit, the reference repeat (8 units), right flank"""
    unit = rand_seq(26, rng)
    return rand_seq(1500, rng), unit, unit * 8, rand_seq(1500, rng)


def test_tandem_repeat_length_variants_kept_apart_at_default_tolerance(rng):
    # child +104 (4 units), mother +121 (4 units and 17 bp): 17 bp apart
    left, unit, ref_tr, right = _tandem_repeat(rng)
    child = left + ref_tr + unit * 4 + right
    mother = left + ref_tr + unit * 4 + unit[:17] + right
    end_tr = len(left) + len(ref_tr)
    group = [
        ins("HG002", "1.0", child, end_tr, 104, 1700),
        ins("HG004", "5.0", mother, end_tr, 121, 1700),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False))
    merged, stats = run(group, cores)
    assert len(merged) == 2
    assert stats.distinct_pairs == 1
    # a tolerance of 20% of the SV size accepts the 17 bp difference
    lenient = cam.CoreAlnParams(frac_tol=0.2)
    merged, _ = run(group, cores, params=lenient)
    assert len(merged) == 1


def test_same_repeat_allele_with_shifted_breakpoints_is_linked(rng):
    left, unit, ref_tr, right = _tandem_repeat(rng)
    child = left + ref_tr + unit * 4 + right
    mother = child  # the same allele
    group = [
        # placed at the start of the repeat on the child, at its end on the mother
        ins("HG002", "1.0", child, len(left), 104, 1500),
        ins("HG004", "5.0", mother, len(left) + len(ref_tr), 104, 1708),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False))
    merged, stats = run(group, cores)
    assert len(merged) == 1
    assert stats.links_alignment == 1


def test_distinct_insertions_at_one_locus_are_not_linked(rng):
    left, right = rand_seq(1500, rng), rand_seq(1500, rng)
    child = left + rand_seq(300, rng) + right
    mother = left + rand_seq(300, rng) + right
    group = [
        ins("HG002", "1.0", child, 1500, 300, 1500),
        ins("HG004", "5.0", mother, 1500, 300, 1500),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False))
    # even when the old rule would merge them, the alignment decides
    merged, stats = run(group, cores, old=always)
    assert len(merged) == 2
    assert stats.links_old_only == 1


def test_two_nearby_deletions_are_not_fused(rng):
    ref = rand_seq(3000, rng)
    # child D56 at 1425, mother D60 at 1570 (145 bp apart)
    child = ref[:1425] + ref[1481:]
    mother = ref[:1570] + ref[1630:]
    group = [
        dele("HG002", "1.0", child, 1425 - 100, ref[1425:1481], 1425),
        dele("HG004", "5.0", mother, 1570 - 300, ref[1570:1630], 1570),
    ]
    cores = cores_of(
        c=("HG002", "1.0", child[100:], False), m=("HG004", "5.0", mother[300:], False)
    )
    merged, _ = run(group, cores, old=always)
    assert len(merged) == 2
    # "or" keeps every link of the old rule
    merged, _ = run(group, cores, mode="or", old=always)
    assert len(merged) == 1


def test_same_deletion_is_linked(rng):
    ref = rand_seq(3000, rng)
    hap = ref[:1500] + ref[1800:]
    child, mother = hap[100:], with_snvs(hap[400:], [50, 2000])
    group = [
        dele("HG002", "1.0", child, 1400, ref[1500:1800], 1500),
        dele("HG004", "5.0", mother, 1100, ref[1500:1800], 1500, reverse=True),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, True))
    merged, stats = run(group, cores)
    assert len(merged) == 1
    assert stats.links_alignment == 1


def test_cores_that_do_not_align_fall_back_to_the_old_rule(rng):
    a, b = rand_seq(3000, rng), rand_seq(3000, rng)
    group = [
        ins("HG002", "1.0", a, 1500, 200, 1500),
        ins("HG004", "5.0", b, 1500, 200, 1500),
    ]
    cores = cores_of(c=("HG002", "1.0", a, False), m=("HG004", "5.0", b, False))
    merged, stats = run(group, cores, old=always)
    assert len(merged) == 1
    assert stats.evaluable_pairs == 0 and stats.links_fallback == 1
    merged, _ = run(group, cores, old=never)
    assert len(merged) == 2


def test_two_alleles_of_one_sample_are_never_joined_through_partners(rng):
    left, unit, ref_tr, right = _tandem_repeat(rng)
    h104 = left + ref_tr + unit * 4 + right
    h130 = left + ref_tr + unit * 5 + right
    end_tr = len(left) + len(ref_tr)
    group = [
        ins("HG002", "7.0", h104, end_tr, 104, 1708),  # child, allele 1
        ins("HG002", "7.1", h130, end_tr, 130, 1708),  # child, allele 2
        ins("HG004", "5.0", h130, end_tr, 130, 1708),  # mother
        ins("HG003", "9.0", h104, end_tr, 104, 1708),  # father
    ]
    cores = cores_of(
        a1=("HG002", "7.0", h104, False),
        a2=("HG002", "7.1", h130, False),
        m=("HG004", "5.0", h130, False),
        f=("HG003", "9.0", h104, False),
    )
    # the old rule says yes to every pair; the alignment says which allele is which
    merged, stats = run(group, cores, mode="or", old=always)
    assert len(merged) == 2
    assert stats.cannot_link_blocked > 0
    samples = sorted(sorted(p.samplename for p in svc.svPatterns) for svc in merged)
    assert samples == [["HG002", "HG003"], ["HG002", "HG004"]]


def test_insertion_shift_range_spans_the_repeat():
    unit = "ACGTTGCAAGT"
    seq = "GGGGCCCC" + unit * 6 + "TTTTAAAA"
    s = 8 + 2 * len(unit)
    lo, hi = cam.insertion_shift_range(seq, s, s + 2 * len(unit))
    assert lo == 8 and hi == 8 + 6 * len(unit)


def test_deletion_shift_range_spans_the_repeat():
    unit = "ACGTTGCAAGT"
    seq = "GGGGCCCC" + unit * 4 + "TTTTAAAA"  # unit * 2 deleted somewhere in here
    lo, hi = cam.deletion_shift_range(seq, 8 + len(unit), unit * 2)
    assert lo == 8 and hi == 8 + 4 * len(unit)


def test_pair_alignment_inversion_round_trip():
    chain = [
        cam.Block(qs=10, qe=110, ts=0, te=105, score=90, cigar=cam.parse_cigar("40=5D30=1X29=")),
        cam.Block(qs=300, qe=400, ts=105, te=205, score=100, cigar=cam.parse_cigar("100=")),
    ]
    aln = cam.PairAlignment.from_chain(chain)
    assert aln.map_pos(10) == 0
    assert aln.map_pos(50) == 45
    assert aln.window_differences(0, 200) == (5 + 190, 1)
    inv = aln.inverted()
    assert inv.map_pos(45) == 50
    assert inv.inverted() == aln


def test_one_deletion_placed_at_two_copies_of_an_imperfect_repeat(rng):
    # five 40-bp copies that differ by one base each; the child's consensus
    # deletes the first copy, the mother's the last: the same haplotype (only
    # mismatches between them), the two deletions 160 bp apart
    unit = rand_seq(40, rng)
    copies = [with_snvs(unit, [5 + 7 * k]) for k in range(5)]
    left, right = rand_seq(1500, rng), rand_seq(1500, rng)
    child = left + "".join(copies[1:]) + right
    mother = left + "".join(copies[:4]) + right
    group = [
        dele("HG002", "1.0", child, len(left), copies[0], 1500),
        dele("HG004", "5.0", mother, len(left) + 160, copies[4], 1660),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False))
    # the haplotypes agree across both deletions and neither lands on the
    # other or on another deletion of the other consensus: one event, linked
    # whatever the pairwise rule says
    merged, stats = run(group, cores, old=never)
    assert len(merged) == 1 and stats.links_twin == 1


def test_two_events_both_haplotypes_carry_are_not_twins(rng):
    # child and mother share a haplotype with two deletions 600 bp apart; each
    # deletion lands on its counterpart, so the cross pairs are not linked
    ref = rand_seq(4000, rng)
    hap = ref[:1500] + ref[1580:2200] + ref[2270:]
    child, mother = hap, hap[200:]
    group = [
        dele("HG002", "1.0", child, 1500, ref[1500:1580], 1500),
        dele("HG002", "1.0", child, 2120, ref[2200:2270], 2200),
        dele("HG004", "5.0", mother, 1300, ref[1500:1580], 1500),
        dele("HG004", "5.0", mother, 1920, ref[2200:2270], 2200),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False))
    merged, stats = run(group, cores, old=always)
    assert records(merged) == [["HG002:1500", "HG004:1300"], ["HG002:2120", "HG004:1920"]]
    assert stats.links_twin == 0


def test_short_and_long_allele_of_an_imperfect_vntr_stay_apart(rng):
    # 18-bp units with one random substitution each: no exact shifting, and a
    # window around the child's +180 alone looks identical in the father's +1080
    unit = rand_seq(18, rng)

    def units(n):
        return "".join(with_snvs(unit, [rng.randrange(18)]) for _ in range(n))

    left, ref_units, extra, right = rand_seq(1500, rng), units(20), units(60), rand_seq(1500, rng)
    child = left + extra[:180] + ref_units + right
    father = left + extra + ref_units + right
    group = [
        ins("HG002", "1.0", child, len(left), 180, 1500),
        ins("HG003", "9.0", father, len(left), 1080, 1500),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), f=("HG003", "9.0", father, False))
    merged, stats = run(group, cores, old=always)
    assert len(merged) == 2
    assert stats.distinct_pairs == 1
