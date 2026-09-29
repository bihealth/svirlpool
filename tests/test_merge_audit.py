"""sv-calling --merge-audit: per-process part files and their concatenation."""

import json
import multiprocessing as mp
from types import SimpleNamespace

import numpy as np
import pytest

from svirlpool.svcalling import merge_audit
from svirlpool.svcalling.svcomposite_merging import _audit_value


@pytest.fixture(autouse=True)
def _reset_audit():
    yield
    merge_audit.AUDIT_PATH = None
    merge_audit._handle = None
    merge_audit._handle_pid = None


def _worker(i: int) -> None:
    merge_audit.record(stage="vertical", i=i)


def test_disabled_is_a_noop(tmp_path):
    merge_audit.record(stage="vertical")
    assert not merge_audit.enabled()
    assert merge_audit.finalize() is None
    assert list(tmp_path.iterdir()) == []


def test_parts_of_forked_workers_are_concatenated(tmp_path):
    out = tmp_path / "audit.jsonl"
    merge_audit.start(out)
    merge_audit.record(stage="horizontal", i=-1)
    with mp.get_context("fork").Pool(3) as pool:
        pool.map(_worker, range(20))
    assert merge_audit.finalize() == out
    rows = [json.loads(line) for line in out.read_text().splitlines()]
    assert sorted(r["i"] for r in rows) == list(range(-1, 20))
    assert list(tmp_path.glob("audit.jsonl.part.*")) == []


def test_start_drops_stale_parts(tmp_path):
    out = tmp_path / "audit.jsonl"
    (tmp_path / "audit.jsonl.part.1").write_text('{"stale":true}\n')
    merge_audit.start(out)
    merge_audit.record(stage="vertical")
    merge_audit.finalize()
    assert "stale" not in out.read_text()


def test_pattern_key():
    svp = SimpleNamespace(
        samplenamed_consensusID="HG002:12.1",
        read_start=100,
        read_end=150,
        get_sv_type=lambda: "DEL",
    )
    assert merge_audit.pattern_key(svp) == "HG002:12.1:100-150:DEL"


def test_audit_values_are_json_safe():
    values = [np.bool_(True), np.int64(3), np.float64(0.5), float("nan"), None]
    assert json.dumps([_audit_value(v) for v in values]) == "[true, 3, 0.5, null, null]"
