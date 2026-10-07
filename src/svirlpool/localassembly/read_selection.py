"""Select the reads of a crowded candidate region that most probably carry its
alleles (--read-selection-factor).

A CR with more than k = factor x the sample's median depth candidate reads is
a pile-up of fragments, collapsed copies of a repeat, or a satellite. Its
reads are ranked by how they cross the CR, from their alignment at the CR and
the alignments in its SA tag:

  crossing   one alignment covers the CR, or two collinear alignments of the
             same strand enclose it (an allele longer than the reference);
             ranked by how far the read runs on beyond the CR on its shorter
             side (`anchor`)
  one-sided  aligned beyond one end of the CR only; ranked by that reach

The k best crossing reads are kept, filled up with one-sided reads. Reads
inside the CR only (short fragments, often of other copies) are dropped. A
crowded CR without a crossing read cannot be spanned and is dropped.

On the HG002 Q100 assembly (experiments/time_limit/README.md) the top crossing
reads of collapsed copies resolved the true allele in 50% of cases, the other
crossing reads in 1%; the fragment reads of a pile-up all came from elsewhere.
"""

from __future__ import annotations

import logging
import re
from typing import NamedTuple

import pysam

log = logging.getLogger(__name__)

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
#: alignments of a read farther than this from the CR do not anchor it
NEAR = 200_000


class Placement(NamedTuple):
    chr: str
    ref_start: int
    ref_end: int
    strand: str
    read_start: int  # on the original read orientation
    read_end: int


def read_placements(aln: pysam.AlignedSegment) -> list[Placement]:
    """The read's alignments: `aln` and those of its SA tag."""
    qlen = aln.infer_read_length() or 0
    entries = [
        (
            aln.reference_name,
            aln.reference_start + 1,
            "-" if aln.is_reverse else "+",
            aln.cigarstring,
        )
    ]
    if aln.has_tag("SA"):
        for entry in str(aln.get_tag("SA")).rstrip(";").split(";"):
            f = entry.split(",")
            if len(f) >= 4:
                entries.append((f[0], int(f[1]), f[2], f[3]))
    out, seen = [], set()
    for chrom, pos, strand, cigar in entries:
        if not cigar or (chrom, pos, strand) in seen:
            continue
        seen.add((chrom, pos, strand))
        ops = [(int(n), op) for n, op in _CIGAR.findall(cigar)]
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        trail = ops[-1][0] if ops[-1][1] in "SH" else 0
        rlen = sum(n for n, op in ops if op in "MDN=X")
        qs, qe = lead, qlen - trail
        if strand == "-":
            qs, qe = qlen - qe, qlen - qs
        out.append(Placement(chrom, pos - 1, pos - 1 + rlen, strand, qs, qe))
    return out


def crossing(
    placements: list[Placement], chrom: str, start: int, end: int
) -> tuple[str, int]:
    """How a read crosses [start, end): ("crossing" | "one-sided" | "inside",
    rank value: the anchor of a crossing read, the reach of a one-sided one)."""
    near = [
        p
        for p in placements
        if p.chr == chrom and p.ref_end > start - NEAR and p.ref_start < end + NEAR
    ]
    best = -1
    for p in near:
        if p.ref_start <= start and p.ref_end >= end:
            best = max(best, min(start - p.ref_start, p.ref_end - end))
    left = [p for p in near if p.ref_start < start]
    right = [p for p in near if p.ref_end > end]
    for a in left:
        for b in right:
            if a is b or a.strand != b.strand or b.ref_end < a.ref_end:
                continue
            # collinear on the read: forward a before b, reverse b before a
            if (a.strand == "+" and a.read_end <= b.read_start + 500) or (
                a.strand == "-" and b.read_end <= a.read_start + 500
            ):
                best = max(best, min(start - a.ref_start, b.ref_end - end))
    if best >= 0:
        return "crossing", best
    reach = [start - p.ref_start for p in left] + [p.ref_end - end for p in right]
    if reach:
        return "one-sided", max(reach)
    return "inside", -1


def select_reads(
    alns: list[pysam.AlignedSegment], chrom: str, start: int, end: int, k: int
) -> tuple[set[str] | None, dict[str, int]]:
    """The reads to keep of one CR (None: all of them), and counts."""
    by_read: dict[str, list[pysam.AlignedSegment]] = {}
    for a in alns:
        by_read.setdefault(a.query_name, []).append(a)
    counts = {
        "reads": len(by_read),
        "crossing": 0,
        "one_sided": 0,
        "kept": len(by_read),
    }
    if k <= 0 or len(by_read) <= k:
        return None, counts
    crossers, one_sided = [], []
    for name, al in by_read.items():
        primary = [a for a in al if not a.is_supplementary] or al
        kind, value = crossing(read_placements(primary[0]), chrom, start, end)
        if kind == "crossing":
            crossers.append((value, name))
        elif kind == "one-sided":
            one_sided.append((value, name))
    counts["crossing"], counts["one_sided"] = len(crossers), len(one_sided)
    if not crossers:
        counts["kept"] = 0
        return set(), counts
    crossers.sort(reverse=True)
    one_sided.sort(reverse=True)
    keep = [n for _, n in crossers[:k]]
    keep += [n for _, n in one_sided[: max(0, k - len(keep))]]
    counts["kept"] = len(keep)
    return set(keep), counts
