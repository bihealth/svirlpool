"""Read retrieval and assembly at an inversion, on simulated ONT-like reads.

A read through an inversion aligns in at least three fragments: left flank,
inverted segment on the other strand, right flank (+/-/+, or -/+/- for a read
from the reverse strand). The consensus code fetches the read from the BAM
(read cache), cuts it to the candidate regions from the extents of all its
fragments, and assembles the cut reads. This checks that the cut read is one
co-linear piece of the inverted haplotype and that both lamassemble methods
assemble the inversion.
"""

import difflib
import random
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from svirlpool.localassembly import consensus
from svirlpool.localassembly.read_cache import ReadSequenceCache
from svirlpool.util import datatypes

pytestmark = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("minimap2", "lastdb", "lastal", "mafft")),
    reason="minimap2, LAST and MAFFT are required",
)

MAT = Path(__file__).parent.parent / "data" / "lamassemble-mats" / "promethion.mat"
METHODS = ("lamassemble", "lamassemble-onestrand")
COMP = str.maketrans("ACGT", "TGCA")
INV_START, INV_LEN = 12_000, 3_000
BP1, BP2 = INV_START, INV_START + INV_LEN
BUFFER_CLIPPED = 500  # svirlpool run's default


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def mutate(rng: random.Random, s: str, rate: float = 0.04) -> str:
    out = []
    for c in s:
        r = rng.random()
        if r < rate * 0.375:
            continue  # deletion
        if r < rate * 0.75:
            out.append(rng.choice("ACGT".replace(c, "")))  # substitution
            continue
        out.append(c)
        if r < rate:
            out.append(rng.choice("ACGT"))  # insertion
    return "".join(out)


def haplotypes(seed: int = 1, repeat_len: int = 0, total: int = 30_000):
    """(reference, inverted haplotype); optionally with inverted repeats at the breakpoints."""
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(total)]
    if repeat_len:
        rep = "".join(rng.choice("ACGT") for _ in range(repeat_len))
        ref[BP1 - repeat_len : BP1] = rep
        ref[BP2 : BP2 + repeat_len] = rc(rep)
    ref = "".join(ref)
    return ref, ref[:BP1] + rc(ref[BP1:BP2]) + ref[BP2:]


def simulate(rng, hap, prefix, n, left, right, overhang=(1500, 4000)):
    """n reads over [left, right) with flanking overhangs, every other one reverse."""
    reads = []
    for i in range(n):
        start = left - rng.randint(*overhang)
        end = right + rng.randint(*overhang)
        s = mutate(rng, hap[start:end])
        reads.append((f"{prefix}{i}", rc(s) if i % 2 else s))
    return reads


def align(tmp: Path, ref: str, reads) -> Path:
    (tmp / "ref.fa").write_text(f">chrT\n{ref}\n")
    (tmp / "reads.fq").write_text(
        "".join(f"@{n}\n{s}\n+\n{'5' * len(s)}\n" for n, s in reads)
    )
    with open(tmp / "reads.sam", "w") as f:
        subprocess.run(
            ["minimap2", "-ax", "map-ont", "-Y", str(tmp / "ref.fa"), str(tmp / "reads.fq")],
            check=True, stdout=f, stderr=subprocess.DEVNULL,
        )  # fmt: skip
    bam = tmp / "reads.bam"
    pysam.sort("-o", str(bam), str(tmp / "reads.sam"))
    pysam.index(str(bam))
    return bam


def hits(tmp: Path, target: str, query: str) -> list[dict]:
    """minimap2 alignments (>= 100 bp) of query on target."""
    (tmp / "t.fa").write_text(f">t\n{target}\n")
    (tmp / "q.fa").write_text(f">q\n{query}\n")
    out = subprocess.run(
        ["minimap2", "-c", "-x", "asm10", str(tmp / "t.fa"), str(tmp / "q.fa")],
        check=True, capture_output=True, text=True,
    ).stdout  # fmt: skip
    result = []
    for line in out.splitlines():
        x = line.split("\t")
        h = {"qlen": int(x[1]), "qs": int(x[2]), "qe": int(x[3]), "strand": x[4],
             "ts": int(x[7]), "te": int(x[8]), "ident": int(x[9]) / int(x[10])}  # fmt: skip
        if h["qe"] - h["qs"] >= 100:
            result.append(h)
    return result


def longest(hs: list[dict]) -> dict | None:
    return max(hs, key=lambda h: h["qe"] - h["qs"], default=None)


def crs_for(layout: str) -> list[datatypes.CandidateRegion]:
    if layout == "one_cr":
        return [datatypes.CandidateRegion(0, "chrT", 0, BP1 - 100, BP2 + 100, [])]
    return [
        datatypes.CandidateRegion(0, "chrT", 0, BP1 - 100, BP1 + 100, []),
        datatypes.CandidateRegion(1, "chrT", 0, BP2 - 100, BP2 + 100, []),
    ]


def retrieve(bam: Path, crs):
    """The read retrieval of process_consensus_container: fetch, interval, cut."""
    alns, records = {}, {}
    with ReadSequenceCache(bam) as cache:
        for cr in crs:
            cr_alns, cr_seqs = cache.fetch_for_cr(cr)
            alns[cr.crID] = cr_alns
            records.update(cr_seqs)
    intervals = consensus.get_read_alignment_intervals_in_cr(
        crs=crs, dict_alignments=alns, buffer_clipped_length=BUFFER_CLIPPED
    )
    cut = consensus.trim_reads(
        dict_alignments=alns,
        intervals=consensus.get_max_extents_of_read_alignments_on_cr(intervals),
        read_records=records,
    )
    return alns, cut


def assemble(tmp: Path, reads: dict, method: str) -> str:
    fa = tmp / f"cluster.{method}.fa"
    fa.write_text("".join(f">{n}\n{r.seq}\n" for n, r in reads.items()))
    seq = consensus.assemble_consensus(MAT, "inv", fa, tmp / "cons.fa", 1, 120, method)
    assert seq, f"{method} produced no consensus"
    return seq


def assert_contains_inversion(tmp: Path, ref: str, alt: str, cons: str):
    # one co-linear alignment to the inverted haplotype across both breakpoints
    h = longest(hits(tmp, alt, cons))
    assert h is not None
    assert (h["qe"] - h["qs"]) >= 0.95 * len(cons)  # ends may be ragged
    assert h["ts"] <= BP1 - 150 and h["te"] >= BP2 + 150
    assert h["ident"] >= 0.97
    # but not to the reference: there the inverted segment breaks it up
    r = longest(hits(tmp, ref, cons))
    assert r is None or (r["qe"] - r["qs"]) < 0.9 * len(cons)


@pytest.fixture(scope="module")
def inversion(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("inversion")
    ref, alt = haplotypes()
    rng = random.Random(2)
    reads = simulate(rng, alt, "alt", 12, BP1, BP2) + simulate(
        rng, ref, "ref", 6, BP1, BP2
    )
    bam = align(tmp, ref, reads)
    fragments: dict[str, list[str]] = {}
    with pysam.AlignmentFile(str(bam)) as f:
        for a in f.fetch():
            fragments.setdefault(a.query_name, []).append("-" if a.is_reverse else "+")
    return tmp, ref, alt, bam, fragments


def test_simulated_reads_align_in_three_fragments(inversion):
    _tmp, _ref, _alt, _bam, fragments = inversion
    alt_frags = [f for n, f in fragments.items() if n.startswith("alt")]
    assert len(alt_frags) == 12
    # the case under test: >= 3 fragments on both strands, on reads of both strands
    assert all(len(f) >= 3 and {"+", "-"} <= set(f) for f in alt_frags)


@pytest.mark.parametrize("layout", ["one_cr", "two_crs"])
def test_cut_reads_are_colinear_pieces_of_their_haplotype(inversion, layout):
    tmp, ref, alt, bam, _ = inversion
    _alns, cut = retrieve(bam, crs_for(layout))
    assert len(cut) == 18
    for name, rec in cut.items():
        hap = alt if name.startswith("alt") else ref
        h = longest(hits(tmp, hap, str(rec.seq)))
        assert h is not None, name
        assert (h["qe"] - h["qs"]) >= 0.95 * len(rec.seq), name
        assert h["ts"] <= BP1 - 150 and h["te"] >= BP2 + 150, name


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("layout", ["one_cr", "two_crs"])
def test_consensus_of_inversion_reads_contains_the_inversion(inversion, layout, method):
    tmp, ref, alt, bam, _ = inversion
    _alns, cut = retrieve(bam, crs_for(layout))
    alt_reads = {n: r for n, r in cut.items() if n.startswith("alt")}
    assert_contains_inversion(tmp, ref, alt, assemble(tmp, alt_reads, method))


def test_one_strand_assembles_inversion_reads_like_both_strands(inversion):
    # not always identical (cross-strand pairs are scored with the forward
    # matrix), but at most a few bases apart
    tmp, _ref, _alt, bam, _ = inversion
    _alns, cut = retrieve(bam, crs_for("one_cr"))
    alt_reads = {n: r for n, r in cut.items() if n.startswith("alt")}
    both = assemble(tmp, alt_reads, METHODS[0])
    one = assemble(tmp, alt_reads, METHODS[1])
    assert abs(len(one) - len(both)) <= 5
    assert difflib.SequenceMatcher(None, one, both, autojunk=False).ratio() > 0.995


@pytest.mark.parametrize("method", METHODS)
def test_unphased_cluster_follows_its_majority_haplotype(inversion, method):
    # 12 inverted + 6 reference reads in one cluster (e.g. the phasing failed)
    tmp, ref, alt, bam, _ = inversion
    _alns, cut = retrieve(bam, crs_for("one_cr"))
    assert_contains_inversion(tmp, ref, alt, assemble(tmp, cut, method))


def oriented_strands(tmp: Path, hap: str, reads: dict) -> dict[str, str]:
    """Strand of each read's longest hit on its haplotype."""
    return {n: longest(hits(tmp, hap, str(r.seq)))["strand"] for n, r in reads.items()}


@pytest.mark.parametrize("seed", [3, 5])
def test_orient_reads_to_reference_puts_inversion_reads_on_one_strand(tmp_path, seed):
    ref, alt = haplotypes()
    rng = random.Random(seed)
    reads = simulate(rng, alt, "alt", 8, BP1, BP2)
    # two reads whose flanks (300-400 bp) are too short to align on their own:
    # their only alignment is the inverted fragment
    reads += [
        (f"short{i}", s)
        for i, (_, s) in enumerate(simulate(rng, alt, "x", 2, BP1, BP2, (300, 400)))
    ]
    bam = align(tmp_path, ref, reads)
    alns, cut = retrieve(bam, crs_for("one_cr"))
    assert all(
        len({a.is_reverse for c in alns for a in alns[c] if a.query_name == n}) == 1
        for n in ("short0", "short1")
    )
    oriented = consensus.orient_reads_to_reference(cut, alns)
    # one strand, and that of the flanks on the reference
    assert set(oriented_strands(tmp_path, alt, oriented).values()) == {"+"}


def test_orient_reads_to_reference_by_alignment_without_inversion(tmp_path):
    # no read aligns on both strands: the strand of the alignment, unchanged
    ref, _alt = haplotypes()
    reads = simulate(random.Random(6), ref, "ref", 6, BP1, BP2)
    bam = align(tmp_path, ref, reads)
    alns, cut = retrieve(bam, crs_for("one_cr"))
    oriented = consensus.orient_reads_to_reference(cut, alns)
    assert set(oriented_strands(tmp_path, ref, oriented).values()) == {"+"}
    for n, r in oriented.items():
        is_reverse = next(
            a.is_reverse for c in alns for a in alns[c] if a.query_name == n
        )
        expected = cut[n].seq.reverse_complement() if is_reverse else cut[n].seq
        assert str(r.seq) == str(expected), n


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="with inverted repeats at the breakpoints longer than the cut flank, the "
    "cut inverted reads equal the reverse complement of the reference window, so "
    "the assembled consensus alone cannot show the inversion (in HG002, e.g. "
    "chr17:5.98 Mb and chr3:187.41 Mb, the padding with read flanks added after "
    "the assembly restores it and the inversion is called)",
)
@pytest.mark.parametrize("method", METHODS)
def test_inversion_between_long_inverted_repeats(tmp_path, method):
    ref, alt = haplotypes(repeat_len=1000)
    rng = random.Random(4)
    reads = simulate(rng, alt, "alt", 10, BP1 - 1000, BP2 + 1000)
    bam = align(tmp_path, ref, reads)
    _alns, cut = retrieve(bam, crs_for("two_crs"))
    assert_contains_inversion(tmp_path, ref, alt, assemble(tmp_path, cut, method))
