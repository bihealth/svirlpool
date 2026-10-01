"""The in-process lamassemble port (svirlpool.localassembly.lamassemble)."""

import dataclasses
import difflib
import shutil
import subprocess
from pathlib import Path

import pytest

from svirlpool.localassembly import consensus
from svirlpool.localassembly import lamassemble as lm

ROOT = Path(__file__).parent.parent
MAT = ROOT / "data" / "lamassemble-mats" / "promethion.mat"
TANDEM_REPEAT = Path(__file__).parent / "data" / "consensus" / "lamassemble_tandem_repeat.fasta"
PARAMS = consensus.lamassemble_params()

needs_tools = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("lastdb", "lastal", "mafft")),
    reason="LAST and MAFFT are required",
)


def _write(path: Path, seqs) -> Path:
    path.write_text("".join(f">{n}\n{s}\n" for n, s in seqs))
    return path


@needs_tools
def test_direct_disttbfast_reproduces_mafft_driver():
    assert lm.direct_mafft_available()
    seqs = lm._self_check_sequences()
    direct = lm.multiple_alignment(seqs, str(MAT), PARAMS, direct=True)
    driver = lm.multiple_alignment(seqs, str(MAT), PARAMS, direct=False)
    assert direct == driver
    assert len(direct) == len(seqs)


@needs_tools
@pytest.mark.parametrize("n", [2, 3, 6])
def test_direct_and_driver_agree_for_few_sequences(n):
    # the mafft driver runs one cycle instead of two for exactly two sequences
    seqs = lm._self_check_sequences()[:n]
    direct = lm.assemble(seqs, MAT, PARAMS, direct=True)
    driver = lm.assemble(seqs, MAT, PARAMS, direct=False)
    assert direct == driver
    assert len(direct) > 300


@needs_tools
@pytest.mark.skipif(shutil.which("lamassemble") is None, reason="lamassemble command not installed")
@pytest.mark.parametrize("which", ["fixture", "tandem_repeat"])
def test_same_consensus_as_lamassemble_command(tmp_path, which):
    if which == "fixture":
        reads = _write(tmp_path / "reads.fa", lm._self_check_sequences())
    else:
        reads = TANDEM_REPEAT
    out = subprocess.run(
        ["lamassemble", "-P", "1", "-f", "fa", "-s", "2", "-g", "67", "-m", "50", str(MAT), str(reads)],
        check=True, capture_output=True, text=True,
    ).stdout  # fmt: skip
    expected = "".join(line for line in out.splitlines() if not line.startswith(">"))
    assert lm.assemble(lm.read_sequences(reads), MAT, PARAMS) == expected


def test_orient_by_kmers_puts_reverse_reads_on_one_strand():
    seqs = [s for _, s in lm._self_check_sequences()]
    flipped = lm.orient_by_kmers(seqs)
    # _self_check_sequences reverse-complements reads 2 and 5
    reference = flipped[0]
    assert [f != reference for f in flipped] == [False, False, True, False, False, True]


@needs_tools
def test_one_strand_consensus_is_close_to_both_strands():
    # not identical: cross-strand pairs are scored with the forward matrix
    seqs = lm._self_check_sequences()
    both = lm.assemble(seqs, MAT, PARAMS)
    one = lm.assemble(seqs, MAT, dataclasses.replace(PARAMS, both_strands=False))
    assert abs(len(one) - len(both)) <= 2
    assert difflib.SequenceMatcher(None, one, both, autojunk=False).ratio() > 0.99


@needs_tools
def test_one_strand_falls_back_to_both_strands_in_tandem_repeat():
    # Oriented alike, these repeat reads give every seed more than -m 50
    # matches and nothing aligns; the fallback must recover the consensus.
    seqs = lm.read_sequences(TANDEM_REPEAT)
    both = lm.assemble(seqs, MAT, PARAMS)
    one = lm.assemble(seqs, MAT, dataclasses.replace(PARAMS, both_strands=False))
    assert len(both) > 1000
    assert one == both


@needs_tools
def test_timeout_raises():
    seqs = lm.read_sequences(TANDEM_REPEAT)
    with pytest.raises(lm.LamassembleTimeout):
        lm.assemble(seqs, MAT, PARAMS, timeout=1e-6)


@needs_tools
def test_make_consensus_with_lamassemble_writes_fasta(tmp_path):
    reads = _write(tmp_path / "reads.fa", lm._self_check_sequences())
    out = tmp_path / "cons.fa"
    seq = consensus.make_consensus_with_lamassemble(
        lamassemble_mat=MAT, reads_file=reads, output=out,
        consensus_name="7.0", threads=1, timeout=60,
    )  # fmt: skip
    assert seq is not None
    assert out.read_text() == f">7.0\n{seq}\n"


def test_make_consensus_with_lamassemble_records_timeout(tmp_path, monkeypatch):
    def slow(*a, **kw):
        raise lm.LamassembleTimeout("timed out")

    monkeypatch.setattr(lm, "assemble", slow)
    reads = _write(tmp_path / "reads.fa", lm._self_check_sequences())
    with consensus.tool_timeouts.watch() as hits:
        seq = consensus.make_consensus_with_lamassemble(
            lamassemble_mat=MAT, reads_file=reads, output=tmp_path / "c.fa",
            consensus_name="x", threads=1, timeout=1,
        )  # fmt: skip
    assert seq is None
    assert hits == ["lamassemble"]


def test_read_sequences_fasta_and_fastq(tmp_path):
    fa = tmp_path / "r.fa"
    fa.write_text(">a desc\nACG\nTT\n\n>b\nGG\n")
    fq = tmp_path / "r.fq"
    fq.write_text("@a\nACGT\n+\nIIII\n@b\nGG\n+\nII\n")
    assert lm.read_sequences(fa) == [("a", "ACGTT"), ("b", "GG")]
    assert lm.read_sequences(fq) == [("a", "ACGT"), ("b", "GG")]


def test_disttbfast_args_follow_the_mafft_driver():
    gap = lm.mafft_gap_options(lm.alignment_scores(str(MAT)))
    args = lm.disttbfast_args(gap, num_seqs=5)
    assert args[args.index("-E") + 1] == "2"
    assert lm.disttbfast_args(gap, num_seqs=2)[args.index("-E") + 1] == "1"
    op = gap[gap.index("--op") + 1]
    assert args[args.index("-f") + 1] == "-" + op
    assert args[args.index("-V") + 1] == "-" + op
