"""Experimental: a cheap haplotype split of a container's reads from the het SNVs
in their reference alignments.

The read phasing (``read_phasing``) aligns the reads all-vs-all, which is the
expensive part of the consensus stage. For the KMeans fast path only one
question needs answering: do the reads fall into two SNV haplotypes, and if
so, do the haplotypes carry different indels? The reads' alignments to the
reference are already loaded, so the het SNVs can be read off them directly:

1. every read contributes its alignment with the largest overlap of the
   window (candidate regions +- ``flank``);
2. a *candidate site* is a reference column outside homopolymers where >= 3
   reads carry one alternative base and >= 3 reads the reference base, with
   the alternative at >= 20% of the aligned bases;
3. a site is kept only if its read split recurs at another site (>= 4 shared
   reads, >= 80% concordant in either phase): all het SNVs of a haplotype
   split the reads alike, sequencing errors do not;
4. the reads are split by the sign of their summed site alleles, with the
   site phases from the leading singular vector and refined by alternating
   majority votes; a read needs a net >= 2 votes to be assigned.

The split exists if both haplotypes have >= 3 reads.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pysam

from .read_phasing import _as_array, _expand, _homopolymer_mask

_Q_ADV = np.array([1, 1, 0, 0, 1, 0, 0, 1, 1, 0], dtype=np.int64)
_R_ADV = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1, 0], dtype=np.int64)
_MATCH_OPS = (0, 7, 8)
# A C G T -> 0..3, everything else 4
_CODE = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i


@dataclass
class RefSnvParams:
    flank: int = 5000
    min_alt: int = 3
    min_ref: int = 3
    min_alt_frac: float = 0.2
    homopolymer_len: int = 4
    min_shared: int = 4
    min_concordance: float = 0.8
    min_votes: int = 2
    min_group: int = 3


@dataclass
class RefSnvSplit:
    n_candidates: int = 0
    n_sites: int = 0
    haplotypes: dict[str, int] = field(default_factory=dict)  # read -> 0 / 1
    min_group: int = 3

    @property
    def is_split(self) -> bool:
        """Both haplotypes have at least ``min_group`` reads."""
        if not self.haplotypes:
            return False
        sizes = np.bincount(list(self.haplotypes.values()), minlength=2)
        return bool(sizes.min() >= self.min_group)


def _aligned_columns(
    aln: pysam.AlignedSegment, w0: int, w1: int
) -> tuple[np.ndarray, np.ndarray]:
    """Reference positions in [w0, w1) and the read's base codes there
    (aligned M/=/X columns only)."""
    ct = np.array(aln.cigartuples, dtype=np.int64)
    ops, lens = ct[:, 0], ct[:, 1]
    rstart = aln.reference_start + np.concatenate(([0], np.cumsum(lens * _R_ADV[ops])[:-1]))
    qstart = np.concatenate(([0], np.cumsum(lens * _Q_ADV[ops])[:-1]))
    m = np.isin(ops, _MATCH_OPS) & (rstart < w1) & (rstart + lens > w0)
    if not m.any():
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int8)
    rpos = _expand(rstart[m], lens[m])
    qpos = _expand(qstart[m], lens[m])
    keep = (rpos >= w0) & (rpos < w1)
    rpos, qpos = rpos[keep], qpos[keep]
    return rpos, _CODE[_as_array(aln.query_sequence)[qpos]]


def _best_alignments(
    alns: list[pysam.AlignedSegment], reads: set[str], chrom: str, w0: int, w1: int
) -> dict[str, pysam.AlignedSegment]:
    best: dict[str, tuple[int, pysam.AlignedSegment]] = {}
    for a in alns:
        if (
            a.query_name not in reads
            or a.reference_name != chrom
            or a.is_secondary
            or a.query_sequence is None
            or a.cigartuples is None
        ):
            continue
        ov = min(a.reference_end, w1) - max(a.reference_start, w0)
        if ov > 0 and (a.query_name not in best or ov > best[a.query_name][0]):
            best[a.query_name] = (ov, a)
    return {rn: a for rn, (_ov, a) in best.items()}


def ref_snv_haplotypes(
    alns: list[pysam.AlignedSegment],
    reads: set[str],
    ref: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    params: RefSnvParams | None = None,
) -> RefSnvSplit:
    """Split ``reads`` into two haplotypes by the het SNVs in their reference
    alignments over ``chrom:start-end`` +- ``params.flank``."""
    p = params or RefSnvParams()
    w0 = max(0, start - p.flank)
    w1 = min(ref.get_reference_length(chrom), end + p.flank)
    best = _best_alignments(alns, reads, chrom, w0, w1)
    if len(best) < 2 * p.min_group:
        return RefSnvSplit(min_group=p.min_group)
    refseq = ref.fetch(chrom, w0, w1).upper()
    ref_code = _CODE[_as_array(refseq)]
    hp = _homopolymer_mask(refseq, p.homopolymer_len)

    names = sorted(best)
    cols = [_aligned_columns(best[rn], w0, w1) for rn in names]
    # base counts per column: 0-3 ACGT, 4 other
    counts = np.zeros((w1 - w0, 5), dtype=np.int32)
    for rpos, code in cols:
        np.add.at(counts, (rpos - w0, code), 1)
    n = w1 - w0
    ref_n = counts[np.arange(n), np.minimum(ref_code, 4)] * (ref_code < 4)
    alt = counts[:, :4].copy()
    alt[np.arange(n)[ref_code < 4], ref_code[ref_code < 4]] = 0
    alt_b = alt.argmax(axis=1)
    alt_n = alt[np.arange(n), alt_b]
    frac = alt_n / np.maximum(1, alt_n + ref_n)
    cand = np.flatnonzero(
        (ref_code < 4)
        & ~hp
        & (alt_n >= p.min_alt)
        & (ref_n >= p.min_ref)
        & (frac >= p.min_alt_frac)
    )
    split = RefSnvSplit(n_candidates=len(cand), min_group=p.min_group)
    if len(cand) < 2:
        return split

    # read x site: +1 alternative, -1 reference, 0 not covered / other base
    M = np.zeros((len(names), len(cand)), dtype=np.int8)
    for i, (rpos, code) in enumerate(cols):
        j = np.searchsorted(cand, rpos - w0)
        hit = (j < len(cand)) & (cand[np.minimum(j, len(cand) - 1)] == rpos - w0)
        j, c = j[hit], code[hit]
        M[i, j[c == alt_b[cand[j]]]] = 1
        M[i, j[c == ref_code[cand[j]]]] = -1

    # recurrence: a site's split must be concordant with another site's
    Mf = M.astype(np.float32)
    agree = Mf.T @ Mf
    shared = np.abs(Mf).T @ np.abs(Mf)
    np.fill_diagonal(shared, 0)
    concordant = (shared >= p.min_shared) & (np.abs(agree) >= p.min_concordance * shared)
    keep = concordant.any(axis=1)
    split.n_sites = int(keep.sum())
    if split.n_sites < 2:
        return split
    Mk = Mf[:, keep]

    # site phases: leading right singular vector, refined by majority votes
    _u, _s, vt = np.linalg.svd(Mk, full_matrices=False)
    v = np.sign(vt[0])
    v[v == 0] = 1
    for _ in range(10):
        score = Mk @ v
        hap = np.sign(score)
        v_new = np.sign(hap @ Mk)
        v_new[v_new == 0] = v[v_new == 0]
        if np.array_equal(v_new, v):
            break
        v = v_new
    score = Mk @ v
    split.haplotypes = {
        rn: int(s < 0) for rn, s in zip(names, score, strict=True) if abs(s) >= p.min_votes
    }
    return split
