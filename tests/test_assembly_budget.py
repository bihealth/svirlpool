"""Deterministic bounds on the consensus cost: --max-assembly-bp, lamassemble's
-m fallback (--lamassemble-max-initial-matches) and the consensus processes'
private temporary directory (--tmp-dir)."""

import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

import svirlpool
from svirlpool.localassembly import consensus, tool_timeouts
from svirlpool.localassembly import lamassemble as lm

ROOT = Path(__file__).parent.parent
MAT = ROOT / "data" / "lamassemble-mats" / "promethion.mat"
TANDEM_REPEAT = (
    Path(__file__).parent / "data" / "consensus" / "lamassemble_tandem_repeat.fasta"
)
needs_tools = pytest.mark.skipif(
    not all(shutil.which(t) for t in ("lastdb", "lastal", "mafft")),
    reason="LAST and MAFFT are required",
)


# --------------------------------------------------------------------------
# --max-assembly-bp


def test_assembly_size_is_unlimited_outside_the_block():
    tool_timeouts.check_assembly_size(10**12)
    with tool_timeouts.assembly_size_limit(0):
        tool_timeouts.check_assembly_size(10**12)


def test_assembly_over_the_limit_raises_and_passes_broad_handlers():
    with tool_timeouts.assembly_size_limit(100):
        tool_timeouts.check_assembly_size(100)  # at the limit: fine
        with pytest.raises(tool_timeouts.AssemblyTooLarge):
            try:
                tool_timeouts.check_assembly_size(101)
            except Exception:  # the consensus code's fallbacks
                pytest.fail("AssemblyTooLarge was caught by `except Exception`")
    tool_timeouts.check_assembly_size(101)  # the limit ends with the block


def test_assemble_consensus_checks_the_size_before_assembling(tmp_path, monkeypatch):
    reads = tmp_path / "reads.fasta"
    reads.write_text(">a\n" + "ACGT" * 30 + "\n>b\n" + "ACGT" * 30 + "\n")  # 240 bp

    def must_not_run(**kwargs):
        raise AssertionError("assembled although over the limit")

    monkeypatch.setattr(consensus, "make_consensus_with_lamassemble", must_not_run)
    with (
        tool_timeouts.assembly_size_limit(239),
        pytest.raises(tool_timeouts.AssemblyTooLarge),
    ):
        consensus.assemble_consensus(
            lamassemble_mat=MAT, name="x", reads_fasta=reads,
            consensus_fasta_path=tmp_path / "c.fasta", threads=1, timeout=10,
            method="lamassemble-onestrand",
        )  # fmt: skip
    assert consensus.fasta_bp(reads) == 240


class _Cache:
    fetch_calls = cache_misses_reads = cache_hits_reads = 0

    def __init__(self, **kwargs):
        pass

    def advance(self, window_start):
        pass

    def close(self):
        pass


def test_container_over_the_assembly_size_is_dropped(tmp_path, monkeypatch):
    containers = {
        crID: {"crs": [SimpleNamespace(crID=crID, chr="chr1", referenceStart=start)]}
        for crID, start in ((1, 100), (2, 5000))
    }
    monkeypatch.setattr(
        consensus, "load_crs_containers_from_db", lambda path_db, crIDs: containers
    )
    monkeypatch.setattr(consensus.read_cache_mod, "ReadSequenceCache", _Cache)
    monkeypatch.setattr(consensus, "available_cpus", lambda: 24)
    calls = []

    def process(crs_dict, threads, timeout, **kwargs):
        crID = next(iter(crs_dict))
        calls.append((crID, threads, timeout))
        # container 1 assembles 5000 bp, container 2 500 bp
        tool_timeouts.check_assembly_size(5000 if crID == 1 else 500)
        return {"c": SimpleNamespace(original_regions=[1])}, {}

    monkeypatch.setattr(consensus, "process_consensus_container", process)
    monkeypatch.setattr(
        consensus.consensus_class,
        "CrsContainerResult",
        lambda consensus_dicts, unused_reads: SimpleNamespace(
            unstructure=lambda: {"n": len(consensus_dicts)}
        ),
    )
    db = tmp_path / "containers.db"
    db.write_text("")
    out = tmp_path / "consensus.jsonl"
    consensus.crs_containers_to_consensus(
        samplename="s",
        input=db,
        output=out,
        lamassemble_mat=None,
        path_alignments=tmp_path / "reads.bam",
        threads=16,
        buffer_clipped_sequence=500,
        consensus_method="lamassemble",
        reference=None,
        escalation=[(1, 20), (4, 60)],
        max_assembly_bp=1000,
    )
    # dropped at the first level, without escalating
    assert calls == [(1, 1, 20), (2, 1, 20)]
    assert [json.loads(line) for line in out.read_text().splitlines()] == [
        {"n": 0},
        {"n": 1},
    ]
    tool_timeouts.check_assembly_size(10**12)  # no limit left behind


# --------------------------------------------------------------------------
# lamassemble's -m fallback


def _attempts(params):
    return [(p.m, p.both_strands, p.m_fallback) for p in lm.layout_attempts(params)]


def test_layout_attempts():
    one = lm.LamassembleParams(m=50, both_strands=False)
    assert _attempts(one) == [(50, False, ()), (50, True, ())]  # as before
    both = lm.LamassembleParams(m=50)
    assert _attempts(both) == [(50, True, ())]  # lamassemble itself
    adaptive = lm.LamassembleParams(m=10, m_fallback=(50,), both_strands=False)
    assert _attempts(adaptive) == [
        (10, False, ()), (10, True, ()), (50, False, ()), (50, True, ()),
    ]  # fmt: skip


def _run_attempts(monkeypatch, params, complete_at):
    """multiple_alignment with LAST, the layout and MAFFT replaced: the layout
    links every read from attempt `complete_at` on. Returns the (m, both
    strands) of the LAST runs."""
    runs = []

    def pairwise(p, scores, sequences, tmpdir, threads, deadline):
        runs.append((p.m, p.both_strands))
        return []

    def layout(p, num_seqs, alignments):
        complete = len(runs) - 1 >= complete_at
        order = [(0, i) for i in range(1, num_seqs)] if complete else []
        return [(0, False)] * num_seqs, order, []

    monkeypatch.setattr(lm, "_pairwise_alignments", pairwise)
    monkeypatch.setattr(lm, "_layout_of_seqs", layout)
    monkeypatch.setattr(lm, "_mafft_alignment", lambda *a, **k: [])
    seqs = [("a", "ACGTACGTAC"), ("b", "ACGTACGTAC"), ("c", "ACGTACGTAC")]
    lm.multiple_alignment(seqs, str(MAT), params, direct=False)
    return runs


def test_multiple_alignment_stops_at_the_first_complete_layout(monkeypatch):
    adaptive = lm.LamassembleParams(m=10, m_fallback=(50,), both_strands=False)
    assert _run_attempts(monkeypatch, adaptive, 0) == [(10, False)]
    assert _run_attempts(monkeypatch, adaptive, 2) == [
        (10, False),
        (10, True),
        (50, False),
    ]
    # never complete: every attempt, the last result is used (as before)
    assert len(_run_attempts(monkeypatch, adaptive, 99)) == 4


def test_lamassemble_params_follow_the_setting():
    assert consensus.parse_max_initial_matches("10,50") == (10, 50)
    with pytest.raises(ValueError):
        consensus.parse_max_initial_matches("10,x")
    with pytest.raises(ValueError):
        consensus.parse_max_initial_matches("0")
    default = consensus.lamassemble_params(both_strands=False)
    assert (default.m, default.m_fallback) == (50, ())
    with consensus.max_initial_matches((10, 50)):
        p = consensus.lamassemble_params(both_strands=False)
        assert (p.m, p.m_fallback, p.both_strands) == (10, (50,), False)
    assert consensus.lamassemble_params().m == 50


@needs_tools
def test_small_m_falls_back_in_a_tandem_repeat():
    # the reads of this short tandem repeat need -m 50 to be linked at all
    seqs = lm.read_sequences(TANDEM_REPEAT)
    params = consensus.lamassemble_params(both_strands=False)
    m50 = lm.assemble(seqs, MAT, params)
    adaptive = lm.assemble(seqs, MAT, lm.LamassembleParams(**{
        **params.__dict__, "m": 10, "m_fallback": (50,)}))  # fmt: skip
    assert len(m50) > 1000
    assert adaptive == m50


# --------------------------------------------------------------------------
# --tmp-dir


def test_private_tmp_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("TMPDIR", "/somewhere/else")
    before = tempfile.tempdir
    with consensus.private_tmp_dir(tmp_path) as private:
        assert private.parent == tmp_path and private.is_dir()
        assert tempfile.gettempdir() == str(private)
        assert os.environ["TMPDIR"] == str(private)
        # the tools started inside inherit it
        out = subprocess.run(
            ["sh", "-c", "echo $TMPDIR"], capture_output=True, text=True
        )
        assert out.stdout.strip() == str(private)
        with tempfile.TemporaryDirectory() as t:
            assert Path(t).parent == private
        (private / "left_over").write_text("x")
    assert not private.exists()
    assert os.environ["TMPDIR"] == "/somewhere/else"
    assert tempfile.tempdir == before


# --------------------------------------------------------------------------
# command lines and the workflow


RUN_ARGS = [
    "run", "--samplename", "s", "--workdir", "w", "--alignments", "a.bam",
    "--reference", "r.fa", "--trf", "t.bed", "--mononucleotides", "m.bed",
    "--threads", "1",
]  # fmt: skip
CONSENSUS_ARGS = [
    "-s", "s", "-i", "c.db", "-a", "a.bam", "-o", "o.jsonl",
    "-r", "r.fa",
]  # fmt: skip


def test_options_and_their_defaults():
    from svirlpool.__main__ import get_parser

    # the defaults since the trio tests of README sections 7-8
    run = get_parser().parse_args(RUN_ARGS)
    assert run.max_assembly_bp == 0
    assert run.lamassemble_max_initial_matches == "10,50"
    assert run.assembly_max_reads == 30
    assert run.container_time_limit == 0
    assert run.consensus_escalation_bp == "100000,300000"
    assert run.consensus_tmp_dir is None
    cons = consensus.get_consensus_parser().parse_args(CONSENSUS_ARGS)
    assert cons.max_assembly_bp == 0
    assert cons.lamassemble_max_initial_matches == "10,50"
    assert cons.assembly_max_reads == 30
    assert cons.container_time_limit == 0
    assert cons.escalation_bp == "100000,300000"
    assert cons.tmp_dir is None


def test_workflow_passes_the_options():
    smk = (Path(svirlpool.__file__).parent / "workflows" / "main.smk").read_text()
    assert re.search(r'max_assembly_bp = config\.get\("max_assembly_bp", 0\)', smk)
    assert "--max-assembly-bp {params.max_assembly_bp}" in smk
    assert (
        "--lamassemble-max-initial-matches {params.lamassemble_max_initial_matches}"
        in smk
    )
    expr = re.search(r"tmp_dir_arg=(.*?),\n", smk).group(1)
    assert eval(expr, {"consensus_tmp_dir": ""}) == ""  # noqa: S307
    assert eval(expr, {"consensus_tmp_dir": "/tmp"}) == "--tmp-dir /tmp"  # noqa: S307

    assert "--escalation-bp {params.escalation_bp}" in smk
    assert 'config.get("consensus_escalation_bp", "100000,300000")' in smk
    assert 'config.get("assembly_max_reads", 30)' in smk
    assert 'config.get("container_time_limit", 0)' in smk


# --------------------------------------------------------------------------
# --escalation-bp: the starting level by assembly size


def test_parse_escalation_bp():
    assert consensus.parse_escalation_bp("100000,300000") == (100000, 300000)
    assert consensus.parse_escalation_bp("0") == ()
    assert consensus.parse_escalation_bp("") == ()
    for bad in ("300000,100000", "-5", "x"):
        with pytest.raises(ValueError):
            consensus.parse_escalation_bp(bad)


def test_required_level():
    with tool_timeouts.at_level(0, 3, (100, 300)):
        assert [tool_timeouts.required_level(bp) for bp in (100, 101, 300, 301)] == [
            0,
            1,
            1,
            2,
        ]
    with tool_timeouts.at_level(0, 2, (100, 300)):
        assert tool_timeouts.required_level(10**6) == 1  # capped at the last level
    with tool_timeouts.at_level(0, 3, ()):
        assert tool_timeouts.required_level(10**6) == 0  # no size rule


def test_require_level_only_escalates_below_the_needed_level():
    with (
        tool_timeouts.at_level(1, 3, (100, 300)),
        tool_timeouts.watch(escalate=True) as hits,
    ):
        tool_timeouts.require_level(200)  # level 2 is the current one
        assert hits == []
        with pytest.raises(tool_timeouts.EscalateTo) as e:
            try:
                tool_timeouts.require_level(400)
            except Exception:  # the consensus code's fallbacks
                pytest.fail("EscalateTo was caught by `except Exception`")
        assert e.value.level == 2
    # outside watch() and at the last level (no escalation): no-op
    with tool_timeouts.at_level(0, 3, (100, 300)):
        tool_timeouts.require_level(10**6)
        with tool_timeouts.watch(escalate=False):
            tool_timeouts.require_level(10**6)


def test_large_assembly_starts_at_its_level(tmp_path, monkeypatch):
    containers = {
        crID: {"crs": [SimpleNamespace(crID=crID, chr="chr1", referenceStart=start)]}
        for crID, start in ((1, 100), (2, 5000), (3, 9000))
    }
    monkeypatch.setattr(
        consensus, "load_crs_containers_from_db", lambda path_db, crIDs: containers
    )
    monkeypatch.setattr(consensus.read_cache_mod, "ReadSequenceCache", _Cache)
    monkeypatch.setattr(consensus, "available_cpus", lambda: 24)
    calls = []
    sizes = {1: 50_000, 2: 200_000, 3: 500_000}

    def process(crs_dict, threads, timeout, **kwargs):
        crID = next(iter(crs_dict))
        calls.append((crID, threads, timeout))
        tool_timeouts.require_level(sizes[crID])  # as assemble_consensus does
        return {"c": SimpleNamespace(original_regions=[1])}, {}

    monkeypatch.setattr(consensus, "process_consensus_container", process)
    monkeypatch.setattr(
        consensus.consensus_class,
        "CrsContainerResult",
        lambda consensus_dicts, unused_reads: SimpleNamespace(
            unstructure=lambda: {"n": len(consensus_dicts)}
        ),
    )
    db = tmp_path / "containers.db"
    db.write_text("")
    out = tmp_path / "consensus.jsonl"
    consensus.crs_containers_to_consensus(
        samplename="s",
        input=db,
        output=out,
        lamassemble_mat=None,
        path_alignments=tmp_path / "reads.bam",
        threads=16,
        buffer_clipped_sequence=500,
        consensus_method="lamassemble",
        reference=None,
        escalation=[(1, 20), (4, 60), (12, 120)],
        escalation_bp=(100_000, 300_000),
    )
    assert calls == [
        (1, 1, 20),  # 50 kb: level 1
        (2, 1, 20), (2, 4, 60),  # 200 kb: straight to level 2
        (3, 1, 20), (3, 12, 120),  # 500 kb: straight to level 3
    ]  # fmt: skip
    assert [json.loads(line) for line in out.read_text().splitlines()] == [{"n": 1}] * 3
