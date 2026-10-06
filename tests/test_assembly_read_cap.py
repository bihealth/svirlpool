"""--assembly-max-reads: the reads an allele is assembled from."""

from types import SimpleNamespace

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus


def test_no_cap_or_few_reads_keeps_all():
    names = ["b", "a", "c"]
    assert consensus.assembly_reads(names, {}, 0) == names
    assert consensus.assembly_reads(names, {}, 3) == names


def test_spanning_reads_are_capped():
    spans = {f"r{i:02d}": [("chr1", 0, 1000)] for i in range(30)}
    spans["short"] = [("chr1", 400, 600)]
    kept = consensus.assembly_reads(list(spans), spans, 10)
    assert len(kept) == 10
    assert "short" not in kept  # spanning reads first, the short one adds nothing
    assert kept == sorted(kept)  # input order kept, ties by name


def test_a_long_region_is_tiled_at_the_cap():
    # reads of 3 kb, 10 starting every 1 kb along a 20 kb region
    spans = {
        f"r{s:02d}_{i}": [("chr1", s * 1000, s * 1000 + 3000)]
        for s in range(18)
        for i in range(10)
    }
    kept = consensus.assembly_reads(list(spans), spans, 6)
    depth: dict[int, int] = {}
    for rn in kept:
        for _, start, end in spans[rn]:
            for b in range(start // 100, end // 100):
                depth[b] = depth.get(b, 0) + 1
    assert len(kept) < len(spans)
    assert set(depth) == set(range(200))
    assert min(depth.values()) >= 6  # every bin of the region keeps the cap


def test_reads_without_spans_are_left_out_unless_none_has_one():
    spans = {f"r{i}": [("chr1", 0, 500)] for i in range(5)}
    assert "x" not in consensus.assembly_reads([*spans, "x"], spans, 3)
    assert consensus.assembly_reads(["x", "y", "z"], {}, 2) == ["x", "y", "z"]


def test_spans_are_clipped_to_the_candidate_regions():
    def aln(name, chrom, start, end):
        return SimpleNamespace(
            query_name=name, reference_name=chrom, reference_start=start,
            reference_end=end,
        )  # fmt: skip

    crs = [
        SimpleNamespace(chr="chr1", referenceStart=10_000, referenceEnd=11_000),
        SimpleNamespace(chr="chr1", referenceStart=12_000, referenceEnd=12_500),
    ]
    alns = {
        0: [aln("a", "chr1", 0, 50_000), aln("b", "chr2", 0, 50_000)],
        1: [aln("c", "chr1", 12_400, 13_000), aln("d", "chr1", 20_000, 21_000)],
    }
    spans = consensus.read_reference_spans(alns, crs, flank=200)
    assert spans == {"a": [("chr1", 9_800, 12_700)], "c": [("chr1", 12_400, 12_700)]}


def test_write_reads_writes_the_subset_next_to_all_reads(tmp_path):
    pool = {
        f"r{i}": SeqRecord(Seq("ACGT" * 10), id=f"r{i}", description="")
        for i in range(8)
    }
    spans = {rn: [("chr1", 0, 1000)] for rn in pool}
    reads_fasta = tmp_path / "reads.phased.0.fasta"
    cap = consensus.AssemblyReadCap(max_reads=5, spans=spans)
    path = cap.write_reads(reads_fasta, list(pool), pool, "1.0")
    assert path == tmp_path / "reads.phased.0.assembly.fasta"
    assert path.read_text().count(">") == 5
    assert cap.write_reads(reads_fasta, ["r0", "r1"], pool, "1.0") == reads_fasta
