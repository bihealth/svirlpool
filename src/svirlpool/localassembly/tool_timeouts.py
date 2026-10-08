"""Record the external-tool timeouts hit while one container is processed.

A tool that runs into its timeout degrades a container's result instead of
failing it: the read phasing falls back to one consensus, a cluster's
lamassemble consensus is dropped, the spectral all-vs-all subsamples the
reads. The batch driver (`consensus.crs_containers_to_consensus`) watches for
these and processes such a container again with more threads and time.

`time_limit` bounds the wall clock of one container over all its levels: a
container still unfinished then is dropped (`--container-time-limit`).

`assembly_size_limit` bounds the bp of any one assembly of a container: a
container about to assemble more is dropped (`--max-assembly-bp`). Unlike the
wall clock, this does not depend on the machine or its load.
"""

from __future__ import annotations

import logging
import signal
import threading
from collections.abc import Iterator
from contextlib import contextmanager

log = logging.getLogger(__name__)

_hits: list[str] | None = None
_escalate = False
_shield_depth = 0  # > 0 inside `shielded()`
_limit_pending = False  # the time limit ran out inside `shielded()`


class Escalate(BaseException):
    """Raised by `record` inside `watch(escalate=True)`: the container is
    abandoned at its first timeout, since it is processed again at the next
    level anyway. A BaseException, so that the broad `except Exception`
    fallbacks of the consensus code, which would finish a degraded result
    first, let it through."""


class EscalateToLast(Escalate):
    """Raised by `skip_to_last` inside `watch(escalate=True)`: the container
    is predicted to time out below the last level, so it goes there directly."""


def skip_to_last(reason: str) -> None:
    """Send the container straight to the last level (a no-op at the last
    level and outside `watch()`)."""
    if _hits is not None and _escalate:
        _hits.append(reason)
        raise EscalateToLast(reason)


def record(tool: str) -> None:
    """Note that `tool` timed out (a no-op outside `watch()`)."""
    if _hits is not None:
        _hits.append(tool)
        if _escalate:
            raise Escalate(tool)


@contextmanager
def watch(escalate: bool = False) -> Iterator[list[str]]:
    """Collect the timeouts recorded inside the block into the yielded list;
    with `escalate`, the first one raises `Escalate`."""
    global _hits, _escalate
    outer = _hits, _escalate
    _hits, _escalate = [], escalate
    try:
        yield _hits
    finally:
        _hits, _escalate = outer


class ContainerTimeLimit(BaseException):
    """Raised inside `time_limit()` once its wall-clock budget is spent. A
    BaseException for the same reason as `Escalate`."""


def _on_alarm(signum, frame) -> None:
    global _limit_pending
    if _shield_depth > 0:
        _limit_pending = True
    else:
        raise ContainerTimeLimit()


@contextmanager
def shielded() -> Iterator[None]:
    """Defer a `ContainerTimeLimit` to the end of the block, for state shared
    beyond one container (the batch's read cache) that must not be left half
    updated."""
    global _shield_depth, _limit_pending
    _shield_depth += 1
    try:
        yield
    finally:
        _shield_depth -= 1
        if _shield_depth == 0 and _limit_pending:
            _limit_pending = False
            raise ContainerTimeLimit()


@contextmanager
def time_limit(seconds: float) -> Iterator[None]:
    """Raise `ContainerTimeLimit` in the block after `seconds` of wall clock
    (0: no limit). The SIGALRM interrupts the Python code wherever it is: a
    running `subprocess.run` / `check_call` kills its child on the way out,
    and the tools' temporary directories are removed by their `with` blocks.
    Code in C (pysam, numpy) is interrupted when it returns. A no-op, with a
    warning, outside the main thread or without SIGALRM."""
    global _limit_pending
    if seconds <= 0:
        yield
        return
    if (
        not hasattr(signal, "setitimer")
        or threading.current_thread() is not threading.main_thread()
    ):
        log.warning("container time limit needs SIGALRM in the main thread; not set")
        yield
        return
    _limit_pending = False
    previous = signal.signal(signal.SIGALRM, _on_alarm)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        _limit_pending = False


class AssemblyTooLarge(BaseException):
    """Raised by `check_assembly_size` inside `assembly_size_limit()` when an
    assembly's reads exceed its bp budget. A BaseException for the same reason
    as `Escalate`: no fallback may assemble the same reads some other way."""


_max_assembly_bp = 0  # 0: no limit


def check_assembly_size(bp: int) -> None:
    """Raise `AssemblyTooLarge` if an assembly of `bp` bp of reads is over the
    budget of the enclosing `assembly_size_limit()` (a no-op outside it)."""
    if 0 < _max_assembly_bp < bp:
        raise AssemblyTooLarge(f"{bp} bp in one assembly (> {_max_assembly_bp})")


@contextmanager
def assembly_size_limit(max_bp: int) -> Iterator[None]:
    """Inside the block, an assembly of more than `max_bp` bp raises
    `AssemblyTooLarge` before it starts (0: no limit). lamassemble's time grows
    with the square of its input, so this bounds a container's cost the way a
    time limit does, but deterministically."""
    global _max_assembly_bp
    outer = _max_assembly_bp
    _max_assembly_bp = max(0, int(max_bp))
    try:
        yield
    finally:
        _max_assembly_bp = outer
