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


def record(tool: str) -> None:
    """Note that `tool` timed out (a no-op outside `watch()`)."""
    if _hits is not None:
        _hits.append(tool)


@contextmanager
def watch() -> Iterator[list[str]]:
    """Collect the timeouts recorded inside the block into the yielded list."""
    global _hits
    outer, _hits = _hits, []
    try:
        yield _hits
    finally:
        _hits = outer
