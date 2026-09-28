"""Record the external-tool timeouts hit while one container is processed.

A tool that runs into its timeout degrades a container's result instead of
failing it: the read phasing falls back to one consensus, a cluster's
lamassemble consensus is dropped, the spectral all-vs-all subsamples the
reads. The batch driver (`consensus.crs_containers_to_consensus`) watches for
these and processes such a container again with more threads and time.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager

_hits: list[str] | None = None
_escalate = False


class Escalate(BaseException):
    """Raised by `record` inside `watch(escalate=True)`: the container is
    abandoned at its first timeout, since it is processed again at the next
    level anyway. A BaseException, so that the broad `except Exception`
    fallbacks of the consensus code, which would finish a degraded result
    first, let it through."""


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
