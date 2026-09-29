"""Structured record of the merge decisions of `sv-calling`, for grading them.

Enabled by `sv-calling --merge-audit <path>`. Every pair the horizontal or the
vertical merge evaluates is written as one JSON object per line: the stage, the
decision and its reason, the features the decision was taken on, and the
pattern keys of both sides. A pattern key (see `pattern_key`) is also written
to the VCF as INFO/PATTERNIDS, so a run without merging
(`--repeat-collapse-mode none --no-vertical-merge`) benchmarked against a truth
set labels every pattern with the truth allele it matches, and each audited
pair can then be graded: same truth allele (should merge), different allele
(should not), or unmatched.

The vertical merge runs in one worker process per chromosome, so every process
appends to its own part file (`<path>.part.<pid>`); `finalize` concatenates the
parts into `<path>` once calling is done.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import IO, Any

#: Set once from the CLI before any worker process is forked; None = off.
AUDIT_PATH: Path | None = None

#: Horizontal candidates further apart than this on the consensus are not
#: recorded: no mode merges them, and there are quadratically many.
HORIZONTAL_WINDOW: int = 1000

_handle: IO[str] | None = None
_handle_pid: int | None = None


def start(path: Path | str) -> None:
    """Enable the audit, dropping part files a previous run left behind."""
    global AUDIT_PATH
    AUDIT_PATH = Path(path)
    for stale in AUDIT_PATH.parent.glob(AUDIT_PATH.name + ".part.*"):
        stale.unlink()


def enabled() -> bool:
    return AUDIT_PATH is not None


def pattern_key(svp: Any) -> str:
    """`sample:crID.subID:start-end:TYPE`, start/end on the consensus sequence.

    Unique per pattern within one svirltile DB, and identical across sv-calling
    runs on that DB, whatever they merge.
    """
    return f"{svp.samplenamed_consensusID}:{svp.read_start}-{svp.read_end}:{svp.get_sv_type()}"


def pattern_keys(svc: Any) -> list[str]:
    return sorted(pattern_key(p) for p in svc.svPatterns)


def record(**fields: Any) -> None:
    """Append one decision. A no-op unless the audit is enabled."""
    global _handle, _handle_pid
    if AUDIT_PATH is None:
        return
    pid = os.getpid()
    if _handle is None or _handle_pid != pid:
        # A forked worker inherits the parent's handle; it must not share it.
        _handle = open(f"{AUDIT_PATH}.part.{pid}", "a")
        _handle_pid = pid
    _handle.write(json.dumps(fields, separators=(",", ":")) + "\n")
    _handle.flush()


def finalize() -> Path | None:
    """Concatenate the per-process parts into AUDIT_PATH; returns it."""
    global _handle, _handle_pid
    if AUDIT_PATH is None:
        return None
    if _handle is not None:
        _handle.close()
        _handle, _handle_pid = None, None
    out = Path(AUDIT_PATH)
    parts = sorted(out.parent.glob(out.name + ".part.*"))
    with open(out, "w") as fout:
        for part in parts:
            with open(part) as fin:
                for line in fin:
                    fout.write(line)
            part.unlink()
    return out
