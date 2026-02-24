"""
Signal Loss Logger

Provides a dedicated logging facility for tracking data loss throughout the
SV signal processing pipeline. Every time signals, primitives, or patterns
are filtered, merged, skipped, or otherwise discarded, an entry is written
to a dedicated log file so users can audit exactly where and why data was
lost or transformed.

Usage in other modules::

    from ..util.signal_loss_logger import get_signal_loss_logger

    loss_logger = get_signal_loss_logger()
    loss_logger.log_filtered(
        stage="parse_sv_signals_from_consensus",
        consensusID="12.0",
        reason="signal outside interval_core",
        details={"sv_type": 0, "ref_start": 100, "ref_end": 200, ...},
    )

At the end of a run, call ``flush()`` or let the process exit normally;
all entries are written to the configured log file.
"""

from __future__ import annotations

import json
import logging
import threading
from pathlib import Path

# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------
_signal_loss_logger: SignalLossLogger | None = None
_lock = threading.Lock()

log = logging.getLogger(__name__)


class SignalLossLogger:
    """Collects and persists records of data loss / filtering events."""

    # Categories of data loss
    FILTERED = "FILTERED"
    MERGED = "MERGED"
    SKIPPED = "SKIPPED"
    EMPTY_INPUT = "EMPTY_INPUT"
    NO_SUPPORT = "NO_SUPPORT"

    def __init__(self, output_path: Path | str | None = None) -> None:
        """
        Parameters
        ----------
        output_path : Path or str, optional
            Path to the signal-loss log file.  If ``None``, a Python logger
            named ``"svirlpool.signal_loss"`` is used instead (writes to
            stderr or wherever the root logger is configured).
        """
        self._output_path = Path(output_path) if output_path is not None else None
        self._file_handle = None
        self._logger = logging.getLogger("svirlpool.signal_loss")
        self._lock = threading.Lock()
        self._entry_count = 0

        if self._output_path is not None:
            self._output_path.parent.mkdir(parents=True, exist_ok=True)
            self._file_handle = open(self._output_path, "w")
            # Write header
            self._file_handle.write(
                "category\tstage\tconsensusID\treason\tcount\tdetails\n"
            )
            self._file_handle.flush()

    # ----- public API -----

    def log_filtered(
        self,
        stage: str,
        reason: str,
        count: int = 1,
        consensusID: str = "",
        details: dict | list | str | None = None,
    ) -> None:
        """Record that data was filtered / discarded."""
        self._write(self.FILTERED, stage, consensusID, reason, count, details)

    def log_merged(
        self,
        stage: str,
        reason: str,
        count: int = 1,
        consensusID: str = "",
        details: dict | list | str | None = None,
    ) -> None:
        """Record that data was merged (individual signals collapsed)."""
        self._write(self.MERGED, stage, consensusID, reason, count, details)

    def log_skipped(
        self,
        stage: str,
        reason: str,
        count: int = 1,
        consensusID: str = "",
        details: dict | list | str | None = None,
    ) -> None:
        """Record that data was skipped (e.g. missing prerequisite)."""
        self._write(self.SKIPPED, stage, consensusID, reason, count, details)

    def log_empty_input(
        self,
        stage: str,
        reason: str,
        count: int = 0,
        consensusID: str = "",
        details: dict | list | str | None = None,
    ) -> None:
        """Record that an empty input was received (nothing to process)."""
        self._write(self.EMPTY_INPUT, stage, consensusID, reason, count, details)

    def log_no_support(
        self,
        stage: str,
        reason: str,
        count: int = 1,
        consensusID: str = "",
        details: dict | list | str | None = None,
    ) -> None:
        """Record that data had no supporting evidence and was dropped."""
        self._write(self.NO_SUPPORT, stage, consensusID, reason, count, details)

    # ----- summary -----

    def get_entry_count(self) -> int:
        return self._entry_count

    def flush(self) -> None:
        with self._lock:
            if self._file_handle is not None:
                self._file_handle.flush()

    def close(self) -> None:
        with self._lock:
            if self._file_handle is not None:
                self._file_handle.close()
                self._file_handle = None

    # ----- internals -----

    def _write(
        self,
        category: str,
        stage: str,
        consensusID: str,
        reason: str,
        count: int,
        details: dict | list | str | None,
    ) -> None:
        details_str = ""
        if details is not None:
            if isinstance(details, (dict, list)):
                try:
                    details_str = json.dumps(details, default=str)
                except (TypeError, ValueError):
                    details_str = str(details)
            else:
                details_str = str(details)

        with self._lock:
            self._entry_count += 1
            line = f"{category}\t{stage}\t{consensusID}\t{reason}\t{count}\t{details_str}\n"
            if self._file_handle is not None:
                self._file_handle.write(line)
                # Flush every write — worker processes spawned by mp.Pool
                # may be terminated without proper cleanup (os._exit), so
                # buffered data would be silently lost.
                self._file_handle.flush()
            else:
                # Fall back to standard Python logging
                self._logger.info(line.rstrip("\n"))


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def init_signal_loss_logger(output_path: Path | str | None = None) -> SignalLossLogger:
    """Initialise (or re-initialise) the global signal-loss logger.

    This should be called once at the start of a pipeline run, before any
    worker processes are spawned.  Worker processes should call
    ``init_signal_loss_logger()`` with their own output path if they want
    per-worker log files, or use the default (Python logger) fallback.
    """
    global _signal_loss_logger
    with _lock:
        if _signal_loss_logger is not None:
            _signal_loss_logger.close()
        _signal_loss_logger = SignalLossLogger(output_path=output_path)
    return _signal_loss_logger


def get_signal_loss_logger() -> SignalLossLogger:
    """Return the global signal-loss logger, creating one if needed."""
    global _signal_loss_logger
    if _signal_loss_logger is None:
        with _lock:
            if _signal_loss_logger is None:
                _signal_loss_logger = SignalLossLogger()  # fallback to Python logger
    return _signal_loss_logger
