"""Experimental: the read phasing of ``read_phasing`` with its sites taken from
the reads' alignments to the REFERENCE instead of from all-vs-all alignments.

This is the ablation arm for the question whether aligning the reads to each
other is worth its cost. Everything after the sites is shared with
``read_phasing`` (``phase_from_sites``: pair weights, correlation clustering,
refinement), and the sites have the same form: one ``Site`` per target read
and variant, with the reads that carry the target's allele and the reads that
carry the other one. Only their source differs:

* SNV sites: reference columns (window = candidate regions +- flank) where a
  target read's base is carried by >= ``min_ref`` reads (the target included)
  and the most frequent other base by >= ``min_alt`` reads, the smaller side
  >= ``min_alt_frac``; columns in reference homopolymers are skipped, and a
  read does not count at a column within ``indel_window`` of an indel of its
  alignment. As in ``read_phasing``, the alleles are the bases the reads
  carry, not reference / alternative. Then the same recurrence filter.
* SV sites: windows around indels >= ``min_sv`` and clipped alignment ends
  (>= ``min_clip``) in any read, chained within ``sv_merge_dist`` and padded by
  ``sv_pad``. A read that spans a window has a net indel there; it agrees with
  the target if the two nets differ by < ``min_sv`` / 2 and disagrees if they
  differ by >= ``min_sv``, which is what their pairwise alignment would show.
  A read clipped inside the window disagrees with a spanning target and agrees
  with a clipped one.
* Low-quality reads: a read's divergence to the reference d(r) is put on the
  pairwise scale of ``read_phasing.low_quality_reads`` (a pair's divergence is
  about the sum of both reads' divergences: D(r) = d(r) + median d), then the
  same rule D(r) > max(lowq_factor * m, m + lowq_min_excess) is applied.

The read set (``select_reads``) and the parameters are those of
``read_phasing``.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pysam
from Bio.SeqRecord import SeqRecord

from .read_phasing import (
    PhasingParams,
    PhasingResult,
    Site,
    _as_array,
    _expand,
    _homopolymer_mask,
    _in_intervals,
    phase_from_sites,
    recurrent_sites,
    select_reads,
)

log = logging.getLogger(__name__)

# query / reference consumed per CIGAR op (M I D N S H P = X)
_Q_ADV = np.array([1, 1, 0, 0, 1, 0, 0, 1, 1, 0], dtype=np.int64)
_R_ADV = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1, 0], dtype=np.int64)
_CODE = np.full(256, 4, dtype=np.int8)  # A C G T -> 0..3, else 4
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i


@dataclass
class RefAln:
    """One alignment of a read to the reference, restricted to the window."""

    r_start: int
    r_end: int
    clip_left: int
    clip_right: int
    col_pos: np.ndarray  # reference positions of aligned (M/=/X) columns
    col_code: np.ndarray  # the read's base code there
    del_start: np.ndarray
    del_len: np.ndarray
    ins_after: np.ndarray  # reference position followed by an insertion
    indel_pos: np.ndarray  # sorted, + insertion / - deletion lengths
    indel_len: np.ndarray
    n_eq: int
    n_x: int
    n_small: int  # bases in indels <= 10 bp


def _parse(aln: pysam.AlignedSegment, ref_code: np.ndarray, w0: int, w1: int) -> RefAln:
    ct = np.array(aln.cigartuples, dtype=np.int64)
    ops, lens = ct[:, 0], ct[:, 1]
    rst = aln.reference_start + np.concatenate(
        ([0], np.cumsum(lens * _R_ADV[ops])[:-1])
    )
    qst = np.concatenate(([0], np.cumsum(lens * _Q_ADV[ops])[:-1]))
    m = np.isin(ops, (0, 7, 8))
    rpos = _expand(rst[m], lens[m])
    qpos = _expand(qst[m], lens[m])
    keep = (rpos >= w0) & (rpos < w1)
    rpos, qpos = rpos[keep], qpos[keep]
    code = _CODE[_as_array(aln.query_sequence)[qpos]]
    rc = ref_code[rpos - w0]
    is_eq = code == rc
    i, d = ops == 1, ops == 2
    in_w = (rst >= w0) & (rst < w1)
    small = (i | d) & (lens <= 10) & in_w
    ins_pos, ins_len = rst[i], lens[i]
    del_pos, del_len = rst[d], lens[d]
    pos = np.concatenate((ins_pos, del_pos))
    ln = np.concatenate((ins_len, -del_len))
    order = np.argsort(pos, kind="stable")
    return RefAln(
        r_start=aln.reference_start,
        r_end=aln.reference_end,
        clip_left=int(lens[0]) if ops[0] in (4, 5) else 0,
        clip_right=int(lens[-1]) if ops[-1] in (4, 5) else 0,
        col_pos=rpos,
        col_code=code,
        del_start=del_pos,
        del_len=del_len,
        ins_after=ins_pos - 1,
        indel_pos=pos[order],
        indel_len=ln[order],
        n_eq=int(is_eq.sum()),
        n_x=int((~is_eq & (code < 4)).sum()),
        n_small=int(lens[small].sum()),
    )


def _blocked(a: RefAln, pos: np.ndarray, w: int) -> np.ndarray:
    starts = np.concatenate((a.del_start - w, a.ins_after - w))
    ends = np.concatenate((a.del_start + a.del_len - 1 + w, a.ins_after + w + 1))
    return _in_intervals(pos, starts, ends)


def low_quality_reads_ref(
    by_read: dict[str, list[RefAln]], params: PhasingParams
) -> set[str]:
    d = {}
    for r, al in by_read.items():
        n_err = sum(a.n_x + a.n_small for a in al)
        n = sum(a.n_eq for a in al) + n_err
        if n > 0:
            d[r] = n_err / n
    if not d:
        return set()
    m_read = float(np.median(list(d.values())))
    pair = {r: v + m_read for r, v in d.items()}  # pairwise scale
    m = float(np.median(list(pair.values())))
    cut = max(params.lowq_factor * m, m + params.lowq_min_excess)
    return {r for r, v in pair.items() if v > cut}


def snv_sites_ref(
    by_read: dict[str, list[RefAln]],
    names: list[str],
    refseq: str,
    w0: int,
    params: PhasingParams,
) -> list[Site]:
    n_col = len(refseq)
    # C[i, c]: read i's usable base code at window column c (-1: none)
    C = np.full((len(names), n_col), -1, dtype=np.int8)
    for i, r in enumerate(names):
        seen = np.zeros(n_col, dtype=bool)
        for a in by_read.get(r, []):
            c = a.col_pos - w0
            ok = (a.col_code < 4) & ~_blocked(a, a.col_pos, params.indel_window)
            twice = seen[c]  # a column covered by two alignments: ambiguous
            C[i, c[ok & ~twice]] = a.col_code[ok & ~twice]
            C[i, c[twice]] = -1
            seen[c] = True
    counts = np.stack([(C == b).sum(axis=0) for b in range(4)], axis=1)  # col x base
    srt = np.sort(counts, axis=1)
    cand = np.flatnonzero(
        (srt[:, -1] >= params.min_alt)
        & (srt[:, -2] >= params.min_ref)
        & ~_homopolymer_mask(refseq, params.homopolymer_len)
    )
    sites: list[Site] = []
    for c in cand:
        cnt = counts[c]
        for b in range(4):
            n_b = int(cnt[b])
            if n_b < params.min_ref:
                continue
            other = cnt.copy()
            other[b] = -1
            alt = int(other.argmax())
            n_alt = int(cnt[alt])
            if (
                n_alt < params.min_alt
                or min(n_b, n_alt) / (n_b + n_alt) < params.min_alt_frac
            ):
                continue
            agree = [names[i] for i in np.flatnonzero(C[:, c] == b)]
            disagree = [names[i] for i in np.flatnonzero(C[:, c] == alt)]
            for t in agree:
                sites.append(
                    Site(
                        t,
                        w0 + int(c),
                        "snv",
                        [t] + [q for q in agree if q != t],
                        disagree,
                        n_b + n_alt,
                    )
                )
    return sites


def sv_sites_ref(
    by_read: dict[str, list[RefAln]],
    names: list[str],
    w0: int,
    w1: int,
    params: PhasingParams,
) -> list[Site]:
    ev: list[int] = []
    for r in names:
        for a in by_read.get(r, []):
            big = np.abs(a.indel_len) >= params.min_sv
            p = a.indel_pos[big]
            ev.extend(p[(p >= w0) & (p < w1)].tolist())
            if a.clip_left >= params.min_clip and w0 <= a.r_start < w1:
                ev.append(a.r_start)
            if a.clip_right >= params.min_clip and w0 <= a.r_end < w1:
                ev.append(a.r_end)
    if not ev:
        return []
    ev.sort()
    wins = [[ev[0], ev[0]]]
    for x in ev[1:]:
        if x - wins[-1][1] <= params.sv_merge_dist:
            wins[-1][1] = x
        else:
            wins.append([x, x])
    sites: list[Site] = []
    for a_, b_ in wins:
        lo, hi = a_ - params.sv_pad, b_ + params.sv_pad
        if lo < w0 or hi > w1:
            continue
        state: dict[str, tuple[str, int]] = {}
        for r in names:
            al = by_read.get(r, [])
            if any(
                (a.clip_left >= params.min_clip and lo <= a.r_start <= hi)
                or (a.clip_right >= params.min_clip and lo <= a.r_end <= hi)
                for a in al
            ):
                state[r] = ("clip", 0)
                continue
            for a in al:
                if a.r_start <= lo and a.r_end >= hi:
                    i0 = np.searchsorted(a.indel_pos, lo, side="left")
                    i1 = np.searchsorted(a.indel_pos, hi, side="right")
                    state[r] = ("span", int(a.indel_len[i0:i1].sum()))
                    break
        for t, (kt, nt) in state.items():
            agree, disagree = [t], []
            for q, (kq, nq) in state.items():
                if q == t:
                    continue
                if kt == "clip" and kq == "clip":
                    agree.append(q)
                elif kt == "clip" or kq == "clip":
                    disagree.append(q)
                elif abs(nq - nt) >= params.min_sv:
                    disagree.append(q)
                elif abs(nq - nt) < params.min_sv / 2:
                    agree.append(q)
            if len(agree) >= 2 and len(disagree) >= 2:
                sites.append(
                    Site(
                        t,
                        (a_ + b_) // 2,
                        "sv",
                        agree,
                        disagree,
                        len(agree) + len(disagree),
                    )
                )
    return sites


def phase_reads_reference(
    reads: dict[str, SeqRecord],
    alns: list[pysam.AlignedSegment],
    ref: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    flank: int,
    params: PhasingParams | None = None,
) -> PhasingResult:
    """Phase ``reads`` (the same read set as ``read_phasing.phase_reads``) by the
    variants in their reference alignments over ``chrom:start-end`` +- flank."""
    params = params or PhasingParams()
    if len(reads) < 2 * params.min_group:
        return PhasingResult(status="no_information", unassigned=sorted(reads))
    all_reads = reads
    reads = select_reads(reads, params)
    w0 = max(0, start - flank)
    w1 = min(ref.get_reference_length(chrom), end + flank)
    refseq = ref.fetch(chrom, w0, w1).upper()
    ref_code = _CODE[_as_array(refseq)]
    by_read: dict[str, list[RefAln]] = defaultdict(list)
    for a in alns:
        if (
            a.query_name not in reads
            or a.reference_name != chrom
            or a.is_secondary
            or a.is_unmapped
            or a.query_sequence is None
            or a.cigartuples is None
            or a.reference_end <= w0
            or a.reference_start >= w1
        ):
            continue
        by_read[a.query_name].append(_parse(a, ref_code, w0, w1))
    # the same alignment can be listed under several CRs of the container
    for r, al in by_read.items():
        uniq = {(a.r_start, a.r_end, a.clip_left, a.clip_right): a for a in al}
        by_read[r] = [uniq[k] for k in sorted(uniq)]
    lowq = low_quality_reads_ref(by_read, params)
    names = sorted(r for r in reads if r not in lowq)
    s_snv = snv_sites_ref(by_read, names, refseq, w0, params)
    if params.recurrence > 0:
        s_snv = recurrent_sites(s_snv, min_support=params.recurrence)
    s_sv = sv_sites_ref(by_read, names, w0, w1, params) if params.use_sv_sites else []
    return phase_from_sites(all_reads, reads, s_snv, s_sv, lowq, params)
