"""
Streaming read-sequence cache for the consensus stage.

A :class:`ReadSequenceCache` keeps one open :class:`pysam.AlignmentFile`
handle and caches full read sequences (as :class:`Bio.SeqRecord.SeqRecord`)
plus the alignments that have been seen, keyed by readname. It is intended
to be used by a driver that iterates over a chromosomally adjacent list of
candidate regions (CRs): every CR triggers at most one ``fetch`` for the
new reads in that interval, and :meth:`advance` evicts cached reads whose
alignments have ended before the next CR's start.

The cache deliberately scopes a single ``pysam`` file handle across all
CRs in a batch so that BAM ``fetch`` calls are amortised across many CRs
instead of being repeated for every container as in the old per-job
design.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

log = logging.getLogger(__name__)


def _seqrecord_from_aln(aln: pysam.AlignedSegment) -> SeqRecord:
    """Build a SeqRecord in original read orientation from an alignment.

    The caller must ensure ``aln`` carries the full read sequence (i.e.
    no hard-clip, ``query_sequence`` length equals ``infer_read_length``).
    """
    seq = (
        Seq(aln.query_sequence).reverse_complement()
        if aln.is_reverse
        else Seq(aln.query_sequence)
    )
    if aln.query_qualities:
        qualities = (
            {"phred_quality": aln.query_qualities[::-1]}
            if aln.is_reverse
            else {"phred_quality": aln.query_qualities}
        )
    else:
        qualities = {"phred_quality": None}
    return SeqRecord(
        seq=seq,
        letter_annotations=qualities,
        id=aln.query_name,
        name=aln.query_name,
    )


class ReadSequenceCache:
    """Sliding read-sequence cache backed by a single open BAM handle.

    Cache state is scoped to a single chromosome: switching to a new
    chromosome triggers a full flush. Within a chromosome the cache
    evicts entries lazily via :meth:`advance` based on the rightmost
    reference end seen for each readname.
    """

    def __init__(self, path_alignments: Path):
        self.path_alignments = Path(path_alignments)
        self._samfile: pysam.AlignmentFile = pysam.AlignmentFile(
            str(self.path_alignments), "rb"
        )
        # readname -> full SeqRecord (original read orientation)
        self._seqs: dict[str, SeqRecord] = {}
        # readname -> list of alignments encountered for that read on the
        # current chromosome (used for downstream supplementary handling
        # and for re-returning alignments belonging to a CR's interval)
        self._alns: dict[str, list[pysam.AlignedSegment]] = {}
        # readname -> rightmost reference_end on current chromosome
        self._max_ref_end: dict[str, int] = {}
        # for dedup across overlapping CR fetches:
        # readname -> set of (flag, reference_start) tuples already stored
        self._aln_keys: dict[str, set[tuple[int, int]]] = {}
        self._current_chr: str | None = None
        # fetch counter for tests / debug
        self.fetch_calls: int = 0
        self.cache_hits_reads: int = 0
        self.cache_misses_reads: int = 0

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def fetch_for_cr(
        self, cr
    ) -> tuple[list[pysam.AlignedSegment], dict[str, SeqRecord]]:
        """Return alignments and full-read sequences for a CR.

        On a chromosome switch the cache is flushed first.
        """
        if self._current_chr != cr.chr:
            if self._current_chr is not None:
                log.debug(
                    "ReadSequenceCache: chr switch %s -> %s, flushing %d cached reads",
                    self._current_chr,
                    cr.chr,
                    len(self._seqs),
                )
            self._flush()
            self._current_chr = cr.chr

        cr_alns: list[pysam.AlignedSegment] = []
        new_readnames: set[str] = set()
        new_alns_for_resolution: list[pysam.AlignedSegment] = []

        self.fetch_calls += 1
        # pysam is compiled with Cython `profile=True`, which does not handle the
        # new PyTrace_STOP_ITERATION (event 10) added in Python 3.12.  Disabling
        # the trace function for the duration of the pysam iterator avoids the
        # resulting ValueError when coverage (or any tracer) is active.
        _tracer = sys.gettrace()
        if _tracer is not None:
            sys.settrace(None)
        try:
            _fetched = list(
                self._samfile.fetch(
                    cr.chr, max(int(cr.referenceStart), 0), int(cr.referenceEnd)
                )
            )
        finally:
            if _tracer is not None:
                sys.settrace(_tracer)
        for aln in _fetched:
            if aln.is_secondary:
                continue
            qname = aln.query_name
            key = (int(aln.flag), int(aln.reference_start))
            seen_keys = self._aln_keys.setdefault(qname, set())
            if key in seen_keys:
                cr_alns.append(_pick_existing_aln(self._alns[qname], key) or aln)
                continue
            seen_keys.add(key)
            self._alns.setdefault(qname, []).append(aln)
            ref_end = aln.reference_end if aln.reference_end is not None else 0
            prev = self._max_ref_end.get(qname, -1)
            if ref_end > prev:
                self._max_ref_end[qname] = ref_end
            cr_alns.append(aln)
            if qname not in self._seqs:
                new_readnames.add(qname)
                new_alns_for_resolution.append(aln)

        if new_readnames:
            self.cache_misses_reads += len(new_readnames)
            self._resolve_new_read_sequences(new_alns_for_resolution, new_readnames)
        # count hits = readnames in cr_alns already present in _seqs at entry
        cr_readnames = {a.query_name for a in cr_alns}
        self.cache_hits_reads += len(cr_readnames - new_readnames)

        cr_seqs = {rn: self._seqs[rn] for rn in cr_readnames if rn in self._seqs}
        log.info(
            f"PROGRESS fetch_for_cr({cr.chr}:{cr.referenceStart}-{cr.referenceEnd}) "
            f"found {len(cr_alns)} alns, {len(new_readnames)} new reads, "
            f"cache now {self.current_size_reads()} reads ({self.current_size_bp()} bp)"
        )
        return cr_alns, cr_seqs

    def advance(self, window_start: int) -> None:
        """Evict any cached read whose alignments end strictly before ``window_start``.

        Call this after a container has been processed, passing the
        leftmost ``referenceStart`` of the next container on the same
        chromosome. Pass ``float('inf')`` (or a very large int) to evict
        everything (e.g. before a chromosome switch).
        """
        if not self._max_ref_end:
            return
        to_evict = [rn for rn, end in self._max_ref_end.items() if end < window_start]
        for rn in to_evict:
            self._seqs.pop(rn, None)
            self._alns.pop(rn, None)
            self._max_ref_end.pop(rn, None)
            self._aln_keys.pop(rn, None)
        if to_evict:
            log.info(
                f"PROGRESS advance(window_start={window_start}) evicted {len(to_evict)} reads, "
                f"cache now {len(self._seqs)} reads ({self.current_size_bp()} bp)"
            )

    def current_size_reads(self) -> int:
        return len(self._seqs)

    def current_size_bp(self) -> int:
        return sum(len(rec.seq) for rec in self._seqs.values())

    def close(self) -> None:
        try:
            self._samfile.close()
        except Exception:  # pragma: no cover - defensive
            pass
        self._flush()

    def __enter__(self) -> "ReadSequenceCache":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------
    def _flush(self) -> None:
        self._seqs.clear()
        self._alns.clear()
        self._max_ref_end.clear()
        self._aln_keys.clear()

    def _resolve_new_read_sequences(
        self,
        new_alns: list[pysam.AlignedSegment],
        new_readnames: set[str],
    ) -> None:
        """Populate ``self._seqs`` for the readnames in ``new_readnames``.

        Most alignments returned by an interval fetch carry the full
        read sequence already (primary, no hard clip). For the rest we
        defer to :func:`consensus.get_full_read_sequences_of_alignments`
        which knows how to follow supplementary alignments. That helper
        opens its own BAM handle, but it is invoked only for the (rare)
        reads whose local fetch returned an alignment lacking the full
        SEQ.
        """
        # Cheap path: alignment already carries the full SEQ.
        deferred: list[pysam.AlignedSegment] = []
        for aln in new_alns:
            qname = aln.query_name
            if qname in self._seqs:
                continue
            if aln.query_sequence and aln.infer_read_length() == len(
                aln.query_sequence
            ):
                self._seqs[qname] = _seqrecord_from_aln(aln)
            else:
                deferred.append(aln)

        # Expensive path: resolve via supplementary scan for the residue.
        still_missing = {
            aln.query_name for aln in deferred if aln.query_name not in self._seqs
        }
        if not still_missing:
            return

        # Lazy import to avoid circular dependency with consensus.py.
        from . import consensus

        deferred_by_qname: dict[int, list[pysam.AlignedSegment]] = {0: []}
        for aln in deferred:
            if aln.query_name in still_missing:
                deferred_by_qname[0].append(aln)
        resolved = consensus.get_full_read_sequences_of_alignments(
            dict_alignments=deferred_by_qname,
            path_alignments=self.path_alignments,
        )
        for qname, rec in resolved.items():
            self._seqs[qname] = rec


def _pick_existing_aln(
    alns: list[pysam.AlignedSegment], key: tuple[int, int]
) -> pysam.AlignedSegment | None:
    flag, ref_start = key
    for a in alns:
        if int(a.flag) == flag and int(a.reference_start) == ref_start:
            return a
    return None
