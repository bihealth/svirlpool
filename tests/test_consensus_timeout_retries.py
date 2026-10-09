"""A container in which a tool timed out is processed again inside its batch,
with more threads and time."""

import json
import subprocess
import time
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus, read_phasing, tool_timeouts


def test_escalation_is_parsed():
    assert consensus.parse_escalation("1:20,4:60,12:120") == [
        (1, 20),
        (4, 60),
        (12, 120),
    ]
    for bad in ("", "1:20,4", "0:20", "1:0", "a:b"):
        with pytest.raises(ValueError):
            consensus.parse_escalation(bad)


def test_escalation_threads_are_capped(monkeypatch):
    monkeypatch.setattr(consensus, "available_cpus", lambda: 8)
    levels = [(1, 20), (4, 60), (12, 120)]
    assert consensus.container_attempts(levels) == [(1, 20), (4, 60), (8, 120)]
    assert consensus.container_attempts(levels, max_threads=2) == [
        (1, 20),
        (2, 60),
        (2, 120),
    ]


def test_timeouts_are_recorded_only_inside_watch():
    tool_timeouts.record("lamassemble")  # no-op
    with tool_timeouts.watch() as hits:
        tool_timeouts.record("lamassemble")
    assert hits == ["lamassemble"]
    with tool_timeouts.watch() as hits:
        pass
    assert hits == []


def test_phasing_timeout_is_recorded(tmp_path, monkeypatch):
    def timing_out(*args, **kwargs):
        raise subprocess.TimeoutExpired("minimap2", 20)

    monkeypatch.setattr(read_phasing, "run_ava", timing_out)
    reads = {
        f"r{i}": SeqRecord(Seq("ACGT" * 50), id=f"r{i}", description="")
        for i in range(10)
    }
    with tool_timeouts.watch() as hits:
        res = read_phasing.phase_reads(reads, tmp_dir_path=tmp_path)
    assert res.status == "failed"
    assert hits == ["phasing all-vs-all"]


class _Cache:
    fetch_calls = cache_misses_reads = cache_hits_reads = 0

    def __init__(self, **kwargs):
        pass

    def advance(self, window_start):
        pass

    def close(self):
        pass


def _run_batch(tmp_path, monkeypatch, timing_out_attempts, escalation=None):
    """Run the batch driver on two containers; container 1 times out on its
    first `timing_out_attempts` attempts. Returns the (crID, threads, timeout)
    of every call and the output lines."""
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
    finished = []  # attempts that went on after their timeout (the fallbacks)

    def process(crs_dict, threads, timeout, **kwargs):
        crID = next(iter(crs_dict))
        calls.append((crID, threads, timeout))
        if crID == 1 and sum(c[0] == 1 for c in calls) <= timing_out_attempts:
            tool_timeouts.record("lamassemble")
            finished.append((crID, threads, timeout))
        return {}, {}

    monkeypatch.setattr(consensus, "process_consensus_container", process)
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
        escalation=escalation or [(1, 20), (4, 60), (12, 120), (32, 300)],
    )
    lines = [json.loads(line) for line in out.read_text().splitlines()]
    _run_batch.finished = finished
    return calls, lines


def test_timed_out_container_escalates_one_level(tmp_path, monkeypatch):
    calls, lines = _run_batch(tmp_path, monkeypatch, timing_out_attempts=1)
    assert calls == [(1, 1, 20), (1, 4, 60), (2, 1, 20)]
    assert len(lines) == 2  # one result per container, escalations replace


def test_escalation_stops_at_the_ceiling(tmp_path, monkeypatch):
    # the last level asks for 32 threads, the batch is capped at 16
    calls, lines = _run_batch(tmp_path, monkeypatch, timing_out_attempts=99)
    assert calls == [(1, 1, 20), (1, 4, 60), (1, 12, 120), (1, 16, 300), (2, 1, 20)]
    assert len(lines) == 2
    # abandoned at the timeout below the last level, the fallbacks run only there
    assert _run_batch.finished == [(1, 16, 300)]


def test_escalate_passes_broad_exception_handlers():
    with tool_timeouts.watch(escalate=True) as hits:
        with pytest.raises(tool_timeouts.Escalate):
            try:
                tool_timeouts.record("lamassemble")
            except Exception:  # the consensus code's fallbacks
                pass
    assert hits == ["lamassemble"]


def test_single_level_does_not_escalate(tmp_path, monkeypatch):
    calls, _ = _run_batch(
        tmp_path, monkeypatch, timing_out_attempts=99, escalation=[(1, 20)]
    )
    assert calls == [(1, 1, 20), (2, 1, 20)]


def test_heavy_container_starts_at_the_last_level(tmp_path, monkeypatch):
    containers = {
        crID: {"crs": [SimpleNamespace(crID=crID, chr="chr1", referenceStart=start)]}
        for crID, start in ((1, 100), (2, 5000))
    }
    monkeypatch.setattr(
        consensus, "load_crs_containers_from_db", lambda path_db, crIDs: containers
    )
    monkeypatch.setattr(consensus.read_cache_mod, "ReadSequenceCache", _Cache)
    monkeypatch.setattr(consensus, "available_cpus", lambda: 24)
    calls, caches = [], []

    def process(crs_dict, threads, timeout, heavy_container_bp, phasing_cache, **kw):
        crID = next(iter(crs_dict))
        calls.append((crID, threads, timeout))
        caches.append(id(phasing_cache))
        if crID == 1 and heavy_container_bp:
            tool_timeouts.skip_to_last("200000 bp of cut reads")
        return {}, {}

    monkeypatch.setattr(consensus, "process_consensus_container", process)
    db = tmp_path / "containers.db"
    db.write_text("")
    consensus.crs_containers_to_consensus(
        samplename="s",
        input=db,
        output=tmp_path / "consensus.jsonl",
        lamassemble_mat=None,
        path_alignments=tmp_path / "reads.bam",
        threads=16,
        buffer_clipped_sequence=500,
        consensus_method="lamassemble",
        reference=None,
        escalation=[(1, 20), (4, 60), (12, 120)],
        heavy_container_bp=100_000,
    )
    # skip_to_last is a no-op at the last level, so container 1 finishes there
    assert calls == [(1, 1, 20), (1, 12, 120), (2, 1, 20)]
    # one phasing cache per container, kept across its levels
    assert caches[0] == caches[1] != caches[2]


def test_skip_to_last_is_a_no_op_outside_escalation():
    tool_timeouts.skip_to_last("big")
    with tool_timeouts.watch(escalate=False) as hits:
        tool_timeouts.skip_to_last("big")
    assert hits == []


def test_time_limit_interrupts_and_passes_broad_handlers():
    t0 = time.monotonic()
    with pytest.raises(tool_timeouts.ContainerTimeLimit):
        with tool_timeouts.time_limit(0.2):
            try:
                time.sleep(5)
            except Exception:  # the consensus code's fallbacks
                pass
    assert time.monotonic() - t0 < 2
    time.sleep(0.3)  # disarmed after the block


def test_time_limit_kills_a_running_tool():
    t0 = time.monotonic()
    with pytest.raises(tool_timeouts.ContainerTimeLimit):
        with tool_timeouts.time_limit(0.2):
            subprocess.run(["sleep", "10"], check=True)
    assert time.monotonic() - t0 < 2


def test_time_limit_waits_for_a_shielded_block():
    done = []
    with pytest.raises(tool_timeouts.ContainerTimeLimit):
        with tool_timeouts.time_limit(0.1):
            with tool_timeouts.shielded():
                time.sleep(0.3)
                done.append(True)
            time.sleep(5)
    assert done == [True]


def test_no_time_limit():
    with tool_timeouts.time_limit(0):
        time.sleep(0.1)


def test_container_over_the_time_limit_is_dropped(tmp_path, monkeypatch):
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
        if crID == 1:
            if len(calls) == 1:
                tool_timeouts.record("lamassemble")  # escalates
            time.sleep(5)  # the next level runs into the limit
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
        container_time_limit=0.3,
    )
    assert calls == [(1, 1, 20), (1, 4, 60), (2, 1, 20)]
    # container 1 dropped (empty result), container 2 kept
    assert [json.loads(line) for line in out.read_text().splitlines()] == [
        {"n": 0},
        {"n": 1},
    ]
    with tool_timeouts.watch(escalate=True) as hits:
        with pytest.raises(tool_timeouts.EscalateToLast):
            tool_timeouts.skip_to_last("big")
    assert hits == ["big"]
