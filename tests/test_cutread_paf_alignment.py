"""Cut reads are aligned to their consensus in PAF mode.

Only the extent, read name and strand of each alignment enter
``Consensus.intervals_cutread_alignments``, so no SAM / base-level alignment is
parsed any more.
"""

import os
import random
import stat
import subprocess
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus, consensus_class
from svirlpool.util import util

PAF_LINE = (
    "read1\t1000\t10\t990\t-\tcons\t2000\t500\t1480\t950\t985\t60\t"
    "tp:A:P\tcm:i:80\ts1:i:900\ts2:i:0\tdv:f:0.01\trl:i:0"
)


def test_parse_paf_line():
    record = util.parse_paf_line(PAF_LINE + "\n")
    assert record == util.PafRecord(
        query_name="read1",
        query_length=1000,
        query_start=10,
        query_end=990,
        is_forward=False,
        target_name="cons",
        target_length=2000,
        target_start=500,
        target_end=1480,
        n_matches=950,
        block_length=985,
        mapq=60,
        alignment_type="P",
    )


def test_parse_paf_line_without_tp_tag_is_primary():
    line = "\t".join(PAF_LINE.split("\t")[:12])
    assert "P" == util.parse_paf_line(line).alignment_type


def test_parse_paf_line_rejects_short_lines():
    with pytest.raises(ValueError, match="PAF"):
        util.parse_paf_line("read1\t1000\t10")


def _record(target, start, end, name, forward, tp="P"):
    return util.PafRecord(
        name, 100, 0, 100, forward, target, 1000, start, end, 90, 100, 60, tp
    )


def test_intervals_are_in_sorted_bam_order():
    """(target, start, forward before reverse), ties in minimap2 order."""
    records = [
        _record("b", 5, 50, "r1", True),
        _record("a", 30, 90, "r2", False),
        _record("a", 30, 80, "r3", True),
        _record("a", 10, 60, "r4", False),
        _record("a", 10, 70, "r5", False),
    ]
    assert [
        (10, 60, "r4", False),
        (10, 70, "r5", False),
        (30, 80, "r3", True),
        (30, 90, "r2", False),
        (5, 50, "r1", True),
    ] == consensus.cutread_intervals_from_paf(records)
    # an explicit target order (FASTA order) wins over the name order
    assert (
        "r1"
        == consensus.cutread_intervals_from_paf(records, target_order={"b": 0, "a": 1})[
            0
        ][2]
    )


# ---------------------------------------------------------------------------
# end to end with minimap2 on synthetic sequences
# ---------------------------------------------------------------------------


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _mutate(seq: str, rng: random.Random, rate: float = 0.02) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < rate / 3:
            continue  # deletion
        if r < 2 * rate / 3:
            out.append(rng.choice("ACGT"))  # substitution
            continue
        out.append(base)
        if r < rate:
            out.append(rng.choice("ACGT"))  # insertion
    return "".join(out)


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    with open(path, "w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n{seq}\n")


@pytest.fixture
def synthetic(tmp_path):
    rng = random.Random(7)
    cons_seq = _random_dna(3000, rng)
    # read name -> (start, end, forward) of the consensus part it carries
    truth = {
        "fwd_full": (0, 3000, True),
        "rev_full": (0, 3000, False),
        "fwd_left": (0, 1800, True),
        "rev_right": (1200, 3000, False),
        "fwd_mid": (700, 2300, True),
    }
    reads = {}
    for name, (start, end, forward) in truth.items():
        seq = _mutate(cons_seq[start:end], rng)
        reads[name] = seq if forward else str(Seq(seq).reverse_complement())
    cons_fa = tmp_path / "consensus.fasta"
    reads_fa = tmp_path / "reads.fasta"
    _write_fasta(cons_fa, {"7.0": cons_seq})
    _write_fasta(reads_fa, reads)
    return cons_fa, reads_fa, cons_seq, truth, reads


def test_minimap_paf_finds_every_read_on_its_strand(synthetic):
    cons_fa, reads_fa, _cons_seq, truth, _reads = synthetic
    records = util.align_reads_with_minimap_paf(
        reference=cons_fa,
        reads=reads_fa,
        aln_args=consensus.CUTREAD_TO_CONSENSUS_ALN_ARGS,
        threads=1,
    )
    assert {r.query_name for r in records} == set(truth)
    for r in records:
        start, end, forward = truth[r.query_name]
        assert r.target_name == "7.0"
        assert r.is_forward is forward
        # the alignment ends lie close to the true extent (reads carry 2% errors)
        assert abs(r.target_start - start) <= 50, r
        assert abs(r.target_end - end) <= 50, r
        assert r.alignment_type != "S"


def test_final_consensus_builds_core_relative_intervals(synthetic):
    cons_fa, reads_fa, cons_seq, truth, _reads = synthetic
    result = consensus.final_consensus(
        reads_fasta=reads_fa,
        consensus_fasta_path=cons_fa,
        consensus_sequence=cons_seq,
        ID="7.0",
        crIDs=[7],
        original_regions=[("chr1", 10_000, 13_000)],
        threads=1,
        verbose=False,
    )
    assert isinstance(result, consensus_class.Consensus)
    assert result.get_used_readnames() == set(truth)
    intervals = result.intervals_cutread_alignments
    assert intervals == sorted(intervals, key=lambda x: (x[0], not x[3]))
    for start, end, name, forward in intervals:
        assert forward is truth[name][2]
        assert 0 <= start < end <= len(cons_seq)
    # reads spanning the middle of the consensus
    spanning = {name for start, end, name, _ in intervals if start <= 1500 <= end}
    assert spanning == set(truth)


def test_final_consensus_without_alignments_returns_none(tmp_path):
    rng = random.Random(3)
    cons_fa = tmp_path / "consensus.fasta"
    reads_fa = tmp_path / "reads.fasta"
    _write_fasta(cons_fa, {"1.0": _random_dna(2000, rng)})
    _write_fasta(reads_fa, {"unrelated": _random_dna(2000, rng)})
    assert (
        consensus.final_consensus(
            reads_fasta=reads_fa,
            consensus_fasta_path=cons_fa,
            consensus_sequence="A",
            ID="1.0",
            crIDs=[1],
            original_regions=[("chr1", 0, 2000)],
            threads=1,
            verbose=False,
        )
        is None
    )


def test_unused_reads_go_to_the_consensus_they_align_to():
    rng = random.Random(11)
    seq_a, seq_b = _random_dna(2500, rng), _random_dna(2500, rng)
    consensus_objects = {
        cid: consensus_class.Consensus(
            ID=cid,
            crIDs=[int(cid.split(".")[0])],
            original_regions=[("chr1", 0, 2500)],
            consensus_sequence=seq,
            intervals_cutread_alignments=[(0, 2500, f"used_{cid}", True)],
        )
        for cid, seq in (("3.0", seq_a), ("3.1", seq_b))
    }
    pool = {
        "to_a": SeqRecord(Seq(_mutate(seq_a[100:2000], rng)), id="to_a", name="to_a"),
        "to_b_rev": SeqRecord(
            Seq(_mutate(seq_b[400:2400], rng)).reverse_complement(),
            id="to_b_rev",
            name="to_b_rev",
        ),
        "nowhere": SeqRecord(Seq(_random_dna(2000, rng)), id="nowhere", name="nowhere"),
    }
    consensus.add_unaligned_reads_to_consensuses_inplace(
        consensus_objects=consensus_objects, pool=pool
    )
    a = consensus_objects["3.0"].intervals_cutread_alignments
    b = consensus_objects["3.1"].intervals_cutread_alignments
    assert ["used_3.0", "to_a"] == [x[2] for x in a]
    assert ["used_3.1", "to_b_rev"] == [x[2] for x in b]
    assert a[1][3] is True and b[1][3] is False
    assert abs(a[1][0] - 100) <= 50 and abs(a[1][1] - 2000) <= 50
    assert abs(b[1][0] - 400) <= 50 and abs(b[1][1] - 2400) <= 50


def test_minimap_paf_timeout_kills_minimap2(tmp_path, monkeypatch):
    """The PAF helper keeps the timeout behaviour of the SAM helper."""
    fake = tmp_path / "bin" / "minimap2"
    fake.parent.mkdir()
    fake.write_text("#!/bin/sh\nsleep 30\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{fake.parent}{os.pathsep}{os.environ['PATH']}")
    with pytest.raises(TimeoutError):
        util.align_reads_with_minimap_paf(
            reference=tmp_path / "ref.fa", reads=tmp_path / "reads.fa", timeout=0.5
        )


def test_minimap_paf_failure_raises(tmp_path, monkeypatch):
    fake = tmp_path / "bin" / "minimap2"
    fake.parent.mkdir()
    fake.write_text("#!/bin/sh\nexit 3\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{fake.parent}{os.pathsep}{os.environ['PATH']}")
    with pytest.raises(subprocess.CalledProcessError):
        util.align_reads_with_minimap_paf(
            reference=tmp_path / "ref.fa", reads=tmp_path / "reads.fa"
        )
