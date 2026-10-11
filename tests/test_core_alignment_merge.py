"""--core-alignment-merge: same allele judged by aligning the consensus cores.

Synthetic haplotypes built from a random reference; the cores of the samples
cover different, overlapping stretches of their haplotype, as the cores of two
samples' candidate regions do. Positions are on the cores (core_start 0).
"""

import random
from types import SimpleNamespace

import pytest

from svirlpool.localassembly import SVpatterns, SVprimitives
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import genotyping
from svirlpool.svcalling.SVcomposite import SVcomposite


def rand_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def with_snvs(seq: str, positions: list[int]) -> str:
    s = list(seq)
    for p in positions:
        s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
    return "".join(s)


def _primitive(
    sample, cid, sv_type, read_start, read_end, ref_start, ref_end, size, reverse
):
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
        SVprimitives=[
            _primitive(sample, cid, 0, rs, re, ref_pos, ref_pos + 1, size, reverse)
        ]
    )
    pattern.set_sequence(seq)
    return SVcomposite.from_SVpattern(pattern)


def dele(sample, cid, core, pos, deleted, ref_pos, reverse=False):
    """A deletion of `deleted` at the point `pos` of the (oriented) core."""
    rs = len(core) - pos if reverse else pos
    pattern = SVpatterns.SVpatternDeletion(
        SVprimitives=[
            _primitive(
                sample,
                cid,
                1,
                rs,
                rs,
                ref_pos,
                ref_pos + len(deleted),
                len(deleted),
                reverse,
            )
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


def run(group, cores, params=None, old=never):
    stats = cam.GroupStats()
    merged, rest = cam.merge_group_with_core_alignments(
        group, cores, old, old, params or cam.CoreAlnParams(), stats
    )
    assert rest == []
    return merged, stats


def records(merged):
    return sorted(
        sorted(f"{p.samplename}:{p.read_start}" for p in svc.svPatterns)
        for svc in merged
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
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
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
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
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
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
    # even when the old rule would merge them, the alignment decides
    merged, stats = run(group, cores, old=always)
    assert len(merged) == 2
    assert stats.evaluable_pairs == 1 and stats.distinct_pairs == 1
    assert stats.links_fallback == 0


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
    merged, stats = run(group, cores, old=always)
    assert len(merged) == 2
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
        cam.Block(
            qs=10,
            qe=110,
            ts=0,
            te=105,
            score=90,
            cigar=cam.parse_cigar("40=5D30=1X29="),
        ),
        cam.Block(
            qs=300, qe=400, ts=105, te=205, score=100, cigar=cam.parse_cigar("100=")
        ),
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
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
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
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
    merged, stats = run(group, cores, old=always)
    assert records(merged) == [
        ["HG002:1500", "HG004:1300"],
        ["HG002:2120", "HG004:1920"],
    ]
    assert stats.links_twin == 0


def test_short_and_long_allele_of_an_imperfect_vntr_stay_apart(rng):
    # 18-bp units with one random substitution each: no exact shifting, and a
    # window around the child's +180 alone looks identical in the father's +1080
    unit = rand_seq(18, rng)

    def units(n):
        return "".join(with_snvs(unit, [rng.randrange(18)]) for _ in range(n))

    left, ref_units, extra, right = (
        rand_seq(1500, rng),
        units(20),
        units(60),
        rand_seq(1500, rng),
    )
    child = left + extra[:180] + ref_units + right
    father = left + extra + ref_units + right
    group = [
        ins("HG002", "1.0", child, len(left), 180, 1500),
        ins("HG003", "9.0", father, len(left), 1080, 1500),
    ]
    cores = cores_of(
        c=("HG002", "1.0", child, False), f=("HG003", "9.0", father, False)
    )
    merged, stats = run(group, cores, old=always)
    assert len(merged) == 2
    assert stats.distinct_pairs == 1


def test_alternating_small_indels_in_a_repeat_count_by_their_net_length():
    # one +13 insertion followed by alternating 2-bp deletions and insertions,
    # as the aligner writes a shifted repeat: 37 indel bases, net +11
    cigar = (
        "50=13I20="
        + "".join("2D9=" for _ in range(6))
        + "".join("2I9=" for _ in range(5))
        + "50="
    )
    ops = cam.parse_cigar(cigar)
    qlen = sum(n for n, op in ops if op in "=I")
    tlen = sum(n for n, op in ops if op in "=D")
    aln = cam.PairAlignment.from_chain([cam.Block(0, qlen, 0, tlen, 100, ops)])
    assert aln.window_differences(0, qlen) == (13 + 12 + 10, 0)
    assert aln.window_differences(0, qlen, cancel_dist=50) == (11, 0)
    # two differences far apart do not cancel
    far = cam.parse_cigar("50=20I300=20D50=")
    aln = cam.PairAlignment.from_chain([cam.Block(0, 390, 0, 390, 100, far)])
    assert aln.window_differences(0, 390, cancel_dist=50) == (40, 0)


def test_two_events_of_one_consensus_are_not_bridged_by_partners(rng):
    # the child's consensus carries two different insertions 400 bp apart; the
    # mother's core covers only the first, the father's only the second, so the
    # parents' cores do not align and the pairwise rule (here: always) links
    # them; without the same-consensus guard the chain child.x1 - mother -
    # father - child.x2 puts the child's two events into one record
    left, mid, right = rand_seq(1500, rng), rand_seq(400, rng), rand_seq(1500, rng)
    i1, i2 = rand_seq(100, rng), rand_seq(120, rng)
    child = left + i1 + mid + i2 + right
    mother = left + i1 + mid[:200]
    father = mid[200:] + i2 + right
    group = [
        ins("HG002", "1.0", child, 1500, 100, 1500),
        ins("HG002", "1.0", child, 2000, 120, 1900),
        ins("HG004", "5.0", mother, 1500, 100, 1500),
        ins("HG003", "9.0", father, 200, 120, 1900),
    ]
    cores = cores_of(
        c=("HG002", "1.0", child, False),
        m=("HG004", "5.0", mother, False),
        f=("HG003", "9.0", father, False),
    )
    merged, stats = run(group, cores, old=always)
    assert stats.same_consensus_cannot == 1
    child_counts = sorted(
        sum(p.samplename == "HG002" for p in svc.svPatterns) for svc in merged
    )
    assert max(child_counts) == 1


def test_pieces_linked_to_one_partner_stay_together(rng):
    # the guard allows two pieces of the child's consensus in one record
    # because the alignment links both to the mother's single insertion
    left, right, inserted = rand_seq(2000, rng), rand_seq(2000, rng), rand_seq(600, rng)
    hap = left + inserted + right
    child, mother = hap[100:], hap[300:]
    group = [
        ins("HG002", "1.0", child, 1900, 250, 2000),
        ins("HG002", "1.0", child, 2150, 350, 2000),
        ins("HG004", "5.0", mother, 1700, 600, 2000),
    ]
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
    merged, stats = run(group, cores, old=never)
    assert len(merged) == 1
    assert stats.same_consensus_pieces == 1 and stats.same_consensus_cannot == 0


def test_guarded_union_keeps_rejected_pairs_of_one_sample_apart(monkeypatch):
    # the rule accepts x-y and y-z but rejects x-z (x and z of one sample):
    # x-y (closer in size) is joined, z stays alone
    class Fake:
        def __init__(self, sample, size):
            self.svPatterns = [SimpleNamespace(samplename=sample)]
            self.size = size

        def get_size(self):
            return self.size

    x, y, z = Fake("HG002", 100), Fake("HG003", 104), Fake("HG002", 115)
    accepted = {(id(x), id(y)), (id(y), id(z))}
    monkeypatch.setattr(
        cam.SVcomposite, "from_SVpatterns", staticmethod(lambda svPatterns: svPatterns)
    )
    out = cam.guarded_union([x, y, z], lambda a, b: (id(a), id(b)) in accepted)
    assert sorted(len(r) for r in out) == [1, 2]
    joined = next(r for r in out if len(r) == 2)
    assert {p.samplename for p in joined} == {"HG002", "HG003"}
    assert all(p is not z.svPatterns[0] for p in joined)


# --------------------------------------------------------------------------- #
# aligner, switches, scheduling
# --------------------------------------------------------------------------- #


def test_cs_string_gives_operations_and_asm20_score():
    cigar, score = cam.cigar_and_score_from_cs(":10*ag+acg:5-tt:3")
    assert cigar == [(10, "="), (1, "X"), (3, "I"), (5, "="), (2, "D"), (3, "=")]
    # 18 matches, one mismatch (-4), gaps of 3 and 2: -min(6 + 2l, 26 + l)
    assert score == 18 - 4 - 12 - 10


def test_align_cores_finds_the_indel_between_two_haplotypes(rng):
    left, right = rand_seq(1500, rng), rand_seq(1500, rng)
    a = left + rand_seq(200, rng) + right
    b = with_snvs(left + right, [700])
    blocks = cam.align_cores({"a": a, "b": b}, cam.CoreAlnParams())
    assert set(blocks) == {("b", "a")}  # query b on target a
    aln = cam.PairAlignment.from_chain(cam.best_colinear_chain(blocks[("b", "a")]))
    # b lacks a's 200 inserted bases: a deletion on the query, one mismatch
    assert aln.window_differences(0, len(b)) == (200, 1)
    for q, t, n in aln.runs:
        mismatches = sum(1 for k in range(n) if b[q + k] != a[t + k])
        assert mismatches == sum(1 for m in aln.mm_q if q <= m < q + n)
    # deterministic
    assert cam.align_cores({"a": a, "b": b}, cam.CoreAlnParams()) == blocks


def test_twins_are_not_linked_when_switched_off(rng):
    unit = rand_seq(40, rng)
    copies = [with_snvs(unit, [5 + 7 * k]) for k in range(5)]
    left, right = rand_seq(1500, rng), rand_seq(1500, rng)
    child = left + "".join(copies[1:]) + right
    mother = left + "".join(copies[:4]) + right
    group = [
        dele("HG002", "1.0", child, len(left), copies[0], 1500),
        dele("HG004", "5.0", mother, len(left) + 160, copies[4], 1660),
    ]
    cores = cores_of(
        c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False)
    )
    merged, stats = run(group, cores, params=cam.CoreAlnParams(link_twins=False))
    assert len(merged) == 2 and stats.links_twin == 0


def test_group_scheduling_does_not_change_the_result(rng):
    from svirlpool.svcalling import svcomposite_merging as M

    left, right, inserted = rand_seq(2000, rng), rand_seq(2000, rng), rand_seq(944, rng)
    hap = left + inserted + right
    child = hap[200 : 2000 + 944 + 1700]
    mother = with_snvs(hap[500 : 2000 + 944 + 1800], [100, 2700])
    father = left[300:] + right[:1500]
    c0 = 2000 - 200
    composites = [
        ins("HG002", "1.0", child, c0, 149, 2000),
        ins("HG002", "1.0", child, c0 + 149, 199, 2000),
        ins("HG002", "1.0", child, c0 + 348, 596, 2000),
        ins("HG004", "5.0", mother, 2000 - 500, 944, 2000),
        dele("HG003", "9.0", father, 1000, right[1500:1600], 4000),
    ]
    core_alignment = cam.CoreAlignmentMerge(
        cores=cores_of(
            c=("HG002", "1.0", child, False),
            m=("HG004", "5.0", mother, False),
            f=("HG003", "9.0", father, False),
        )
    )

    def merged_records(threads):
        merged = M.merge_svComposites(
            svComposites=composites,
            apriori_size_difference_fraction_tolerance=0.1,
            near=150,
            min_kmer_overlap=0.7,
            threads=threads,
            core_alignment=core_alignment,
        )
        return [
            sorted(f"{p.samplename}:{p.read_start}" for p in svc.svPatterns)
            for svc in merged
        ]

    one = merged_records(1)
    assert one == merged_records(3)
    assert sorted(len(r) for r in one) == [1, 4]


def test_cores_are_read_from_the_db_and_old_headers_are_skipped(tmp_path):
    import pickle
    import sqlite3

    db = tmp_path / "svirltile.db"
    conn = sqlite3.connect(db)
    conn.execute(
        "CREATE TABLE consensus_sequences "
        "(consensusID TEXT PRIMARY KEY, sequence BLOB, description TEXT)"
    )
    rows = [
        ("1.0", "aaCCGGTTaa", "1.0 core=2-8"),
        ("2.0", "ACGTACGT", "2.0"),  # written before the headers named the core
        ("3.0", "ACGT", "3.0 core=0-4"),
    ]
    conn.executemany(
        "INSERT INTO consensus_sequences VALUES (?, ?, ?)",
        [(cid, pickle.dumps(seq), desc) for cid, seq, desc in rows],
    )
    conn.commit()
    conn.close()
    cores, without_core = cam.load_core_sequences(db, {"1.0", "2.0"})
    assert cores == {"1.0": (2, "CCGGTT")}
    assert without_core == 1


@pytest.mark.parametrize("crop_margin", [None, 300])
def test_cropped_cores_give_the_same_decisions(rng, crop_margin):
    # the decomposed insertion of the first test, with SNVs far from the SVs
    # that only the uncropped alignment sees
    left, right, inserted = rand_seq(2000, rng), rand_seq(2000, rng), rand_seq(944, rng)
    hap = left + inserted + right
    child = with_snvs(hap[200 : 2000 + 944 + 1700], [50, 4200])
    mother = with_snvs(hap[500 : 2000 + 944 + 1800], [100, 2700])
    c0 = 2000 - 200
    group = [
        ins("HG002", "1.0", child, c0, 149, 2000),
        ins("HG002", "1.0", child, c0 + 149, 199, 2000),
        ins("HG002", "1.0", child, c0 + 348, 596, 2000),
        ins("HG004", "5.0", mother, 2000 - 500, 944, 2000, reverse=True),
    ]
    cores = cores_of(c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, True))
    params = cam.CoreAlnParams(crop_margin=crop_margin)
    merged, stats = run(group, cores, params=params)
    assert len(merged) == 1 and len(merged[0].svPatterns) == 4
    assert stats.links_alignment == 3
    # the blocks come back in whole-core coordinates
    placed, oriented = cam.place_composites(group, cores)
    keys = {p.key for p in placed}
    windows = {
        k: (
            min(p.amb_start for p in placed if p.key == k),
            max(p.amb_end for p in placed if p.key == k),
        )
        for k in keys
    }
    alns = cam.align_group_cores(
        oriented, keys, params, cam.GroupStats(), windows=windows
    )
    ck, mk = ("HG002", "1.0"), ("HG004", "5.0")
    aln = alns[(ck, mk)]
    for q, t, n in aln.runs:
        assert q + n <= len(oriented[ck]) and t + n <= len(oriented[mk])
        assert sum(
            a != b
            for a, b in zip(
                oriented[ck][q : q + n], oriented[mk][t : t + n], strict=True
            )
        ) == sum(1 for m in aln.mm_q if q <= m < q + n)
    if crop_margin is not None:
        assert aln.q0 >= windows[ck][0] - crop_margin
        assert aln.q1 <= windows[ck][1] + crop_margin
