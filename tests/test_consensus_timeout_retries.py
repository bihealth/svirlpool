"""A container in which a tool timed out is processed again inside its batch,
with more threads and time."""

import json
import subprocess
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

    def process(crs_dict, threads, timeout, **kwargs):
        crID = next(iter(crs_dict))
        calls.append((crID, threads, timeout))
        if crID == 1 and sum(c[0] == 1 for c in calls) <= timing_out_attempts:
            tool_timeouts.record("lamassemble")
        return {}, {}

    monkeypatch.setattr(consensus, "process_consensus_container", process)
    db = tmp_path / "containers.db"
    db.write_text("")
    out = tmp_path / "consensus.jsonl"
    consensus.crs_containers_to_consensus(
        samplename="s",
        input=db,
        copy_number_tracks=tmp_path / "cn.bed.gz",
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


def test_single_level_does_not_escalate(tmp_path, monkeypatch):
    calls, _ = _run_batch(
        tmp_path, monkeypatch, timing_out_attempts=99, escalation=[(1, 20)]
    )
    assert calls == [(1, 1, 20), (2, 1, 20)]
