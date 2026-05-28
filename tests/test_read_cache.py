"""Tests for :mod:`svirlpool.localassembly.read_cache`."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pysam

from svirlpool.localassembly.read_cache import ReadSequenceCache

# ---------------------------------------------------------------------------
# fixtures: build a tiny BAM file with controllable read layout
# ---------------------------------------------------------------------------


def _make_bam(tmp_path: Path, reads: list[dict]) -> Path:
    """Create a small indexed BAM. ``reads`` is a list of dicts with keys
    ``name``, ``chr``, ``start``, ``seq``, optional ``flag`` (default 0).
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": 1_000_000},
            {"SN": "chr2", "LN": 1_000_000},
        ],
    }
    chr_to_tid = {"chr1": 0, "chr2": 1}
    bam_path = tmp_path / "test.bam"
    # Sort reads by (chr, start) for the coordinate-sorted BAM.
    sorted_reads = sorted(reads, key=lambda r: (chr_to_tid[r["chr"]], r["start"]))
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        for r in sorted_reads:
            a = pysam.AlignedSegment(f.header)
            a.query_name = r["name"]
            a.query_sequence = r["seq"]
            a.flag = int(r.get("flag", 0))
            a.reference_id = chr_to_tid[r["chr"]]
            a.reference_start = int(r["start"])
            a.mapping_quality = 60
            a.cigartuples = [(0, len(r["seq"]))]  # all M
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r["seq"]))
            f.write(a)
    pysam.index(str(bam_path))
    return bam_path


@dataclass
class _CR:
    """Minimal stand-in for ``datatypes.CandidateRegion`` (only used attrs)."""

    crID: int
    chr: str
    referenceStart: int
    referenceEnd: int


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------


def test_cache_dedup_across_overlapping_crs(tmp_path: Path) -> None:
    bam = _make_bam(
        tmp_path,
        [
            {"name": "r1", "chr": "chr1", "start": 100, "seq": "A" * 200},
            {"name": "r2", "chr": "chr1", "start": 150, "seq": "C" * 200},
            {"name": "r3", "chr": "chr1", "start": 800, "seq": "G" * 100},
        ],
    )
    with ReadSequenceCache(bam) as cache:
        alns1, seqs1 = cache.fetch_for_cr(_CR(1, "chr1", 100, 300))
        alns2, seqs2 = cache.fetch_for_cr(_CR(2, "chr1", 200, 400))
        # r1 and r2 are returned to both CRs but only loaded once.
        assert {a.query_name for a in alns1} == {"r1", "r2"}
        assert {a.query_name for a in alns2} == {"r1", "r2"}
        assert set(seqs1) == {"r1", "r2"}
        assert set(seqs2) == {"r1", "r2"}
        assert cache.cache_misses_reads == 2
        assert cache.cache_hits_reads == 2


def test_cache_sliding_window_eviction(tmp_path: Path) -> None:
    bam = _make_bam(
        tmp_path,
        [
            {"name": "r1", "chr": "chr1", "start": 100, "seq": "A" * 200},  # ends 300
            {"name": "r2", "chr": "chr1", "start": 150, "seq": "C" * 400},  # ends 550
            {"name": "r3", "chr": "chr1", "start": 800, "seq": "G" * 100},
        ],
    )
    with ReadSequenceCache(bam) as cache:
        cache.fetch_for_cr(_CR(1, "chr1", 100, 300))
        assert cache.current_size_reads() == 2
        # Advance past r1's end (300) but not r2's end (550).
        cache.advance(window_start=400)
        assert cache.current_size_reads() == 1
        assert "r2" in cache._seqs and "r1" not in cache._seqs
        # Advance past r2's end too.
        cache.advance(window_start=600)
        assert cache.current_size_reads() == 0


def test_cache_flushes_on_chromosome_switch(tmp_path: Path) -> None:
    bam = _make_bam(
        tmp_path,
        [
            {"name": "r1", "chr": "chr1", "start": 100, "seq": "A" * 200},
            {"name": "r2", "chr": "chr2", "start": 100, "seq": "T" * 200},
        ],
    )
    with ReadSequenceCache(bam) as cache:
        cache.fetch_for_cr(_CR(1, "chr1", 100, 300))
        assert cache.current_size_reads() == 1
        cache.fetch_for_cr(_CR(2, "chr2", 100, 300))
        # chr1 reads should be flushed.
        assert "r1" not in cache._seqs
        assert "r2" in cache._seqs
        assert cache.current_size_reads() == 1


def test_cache_reverse_complement_preserves_original_orientation(
    tmp_path: Path,
) -> None:
    # flag 16 = reverse-strand alignment. The BAM stores the reverse
    # complement of the read; the cache must return the original read
    # orientation.
    seq_on_ref = "ACGTACGTACGT"  # 12 bp; stored in BAM verbatim
    bam = _make_bam(
        tmp_path,
        [
            {
                "name": "rfwd",
                "chr": "chr1",
                "start": 100,
                "seq": seq_on_ref,
                "flag": 0,
            },
            {
                "name": "rrev",
                "chr": "chr1",
                "start": 200,
                "seq": seq_on_ref,
                "flag": 16,
            },
        ],
    )
    with ReadSequenceCache(bam) as cache:
        _alns, seqs = cache.fetch_for_cr(_CR(1, "chr1", 50, 300))
        assert str(seqs["rfwd"].seq) == seq_on_ref
        # reverse-strand read: cache should reverse-complement back.
        assert str(seqs["rrev"].seq) == "ACGTACGTACGT"[::-1].translate(
            str.maketrans("ACGT", "TGCA")
        )


def test_cache_skips_secondary_alignments(tmp_path: Path) -> None:
    # flag 256 = secondary; cache must ignore it.
    bam = _make_bam(
        tmp_path,
        [
            {"name": "r1", "chr": "chr1", "start": 100, "seq": "A" * 200},
            {
                "name": "r1",
                "chr": "chr1",
                "start": 150,
                "seq": "A" * 200,
                "flag": 256,
            },
        ],
    )
    with ReadSequenceCache(bam) as cache:
        alns, _seqs = cache.fetch_for_cr(_CR(1, "chr1", 50, 400))
        # only the primary alignment is kept
        assert len(alns) == 1
        assert int(alns[0].flag) == 0
