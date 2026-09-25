"""Experimental: phase the reads of a candidate-region container from their
all-vs-all alignments, and derive the number of alleles from the phasing.

The legacy clustering in ``consensus`` compares reads by their SV-scale
signals only and needs the local copy number as the number of clusters.  Reads
of one haplotype, however, also share every small variant -- above all SNVs --
over their whole length, and the reads are much longer than the candidate
region.  This module aligns reads cut to the candidate region +- a phasing
flank all-vs-all (minimap2 ``--eqx``; every pair is aligned once and the
other direction is derived by inverting the alignment, so every read is the
target of all others -- the caller puts the reads on one strand) and then

1. finds *informative sites* on every target read t: a column where >= 3 other
   reads share one alternative base and >= 3 reads (t included) carry t's base,
   outside homopolymers and away from indels.  A sequencing error rarely recurs
   in several reads with the same base; a heterozygous SNV does.
2. keeps a site only if its split of the reads *recurs* at another position of
   the same target (all het SNVs of a haplotype split the reads identically;
   coinciding errors form random subsets).
3. adds SV-scale sites: windows of t where a read's net indel is >= 20 bp or
   its alignment is clipped (large SV / breakend).
4. scores every read pair by agreements and disagreements over all sites (an
   agreement is weighted by how balanced the site's split is), and clusters
   the signed weights by greedy correlation clustering, which needs no number
   of clusters.
5. merges clusters that are not separated by >= 2 discriminating variant
   positions (both near-unanimous, different alleles) and assigns reads of tiny
   clusters to the cluster they have positive evidence for.

The number of remaining clusters is the number of alleles.  Reads with a much
higher divergence to their partners than the pool (failed / very noisy reads)
are excluded: a few of them share wrong bases often enough to pose as a
haplotype.

Validated on HG002 20x ONT (2004 containers; trio-derived read haplotypes as
truth), see ``experiments/ava_phasing``.
"""

from __future__ import annotations

import logging
import shlex
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import tool_timeouts

log = logging.getLogger(__name__)

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class PhasingParams:
    # minimap2: the ava-ont preset without -X (all pairs in both directions)
    minimap_params: str = "-k15 -w5 -e0 -m100 -r2k"
    # align every read pair once (--dual=no) and derive the other direction
    # by inverting the alignment: half the alignment work
    align_pairs_once: bool = True
    # SNV sites
    min_alt: int = 3
    min_ref: int = 2
    min_alt_frac: float = 0.2
    homopolymer_len: int = 4
    indel_window: int = 2
    recurrence: int = 1
    # SV sites
    use_sv_sites: bool = True
    min_sv: int = 20
    sv_merge_dist: int = 100
    sv_pad: int = 50
    min_clip: int = 100
    # reads
    # phase at most this many reads (the longest); the all-vs-all cost grows
    # with the square of the read count. The others stay unassigned.
    max_reads: int = 50
    lowq_factor: float = 3.0
    lowq_min_excess: float = 0.05
    # clustering
    error_rate: float = 0.1
    min_group: int = 3
    min_discriminating: int = 2
    max_alleles: int = 4


@dataclass
class PhasingResult:
    status: str  # "phased", "single", "no_information", "failed"
    groups: dict[str, int] = field(default_factory=dict)  # readname -> allele
    n_alleles: int = 0
    unassigned: list[str] = field(default_factory=list)
    low_quality: list[str] = field(default_factory=list)
    n_snv_sites: int = 0
    n_sv_sites: int = 0

    def clusters(self) -> dict[int, list[str]]:
        out: dict[int, list[str]] = defaultdict(list)
        for r, g in sorted(self.groups.items()):
            out[g].append(r)
        return dict(sorted(out.items()))


# --------------------------------------------------------------------------- #
# all-vs-all alignments
# --------------------------------------------------------------------------- #
@dataclass
class PairAln:
    """Read q aligned to read t, in t's coordinates (numpy arrays)."""

    q: str
    t: str
    t_start: int
    t_end: int
    mism_pos: np.ndarray  # t positions of X columns
    mism_base: np.ndarray  # q's base there (uint8 ASCII)
    del_start: np.ndarray  # deletions in q: runs of t positions
    del_len: np.ndarray
    ins_after: np.ndarray  # t positions followed by an insertion in q
    indel_pos: np.ndarray  # sorted t positions of indels ...
    indel_len: np.ndarray  # ... with +insertion / -deletion length
    clip_left: int
    clip_right: int
    divergence: float  # (X + indel bases of indels <= 10 bp) per aligned base

    # dict / set views (tests, debugging)
    @property
    def mism(self) -> dict[int, str]:
        return {
            int(p): chr(b) for p, b in zip(self.mism_pos, self.mism_base, strict=True)
        }

    @property
    def dels(self) -> set[int]:
        return {
            int(p)
            for s, n in zip(self.del_start, self.del_len, strict=True)
            for p in range(s, s + n)
        }

    @property
    def indels(self) -> list[tuple[int, int]]:
        return [
            (int(p), int(n))
            for p, n in zip(self.indel_pos, self.indel_len, strict=True)
        ]

    def blocked_at(self, pos: np.ndarray, w: int) -> np.ndarray:
        """Whether each t position is within w of a deletion or next to an
        insertion (+-w): its column is unreliable in this alignment."""
        starts = np.concatenate((self.del_start - w, self.ins_after - w))
        ends = np.concatenate(
            (self.del_start + self.del_len - 1 + w, self.ins_after + w + 1)
        )
        return _in_intervals(pos, starts, ends)


def _in_intervals(x: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """x in the union of the closed intervals [starts, ends]."""
    if len(starts) == 0:
        return np.zeros(len(x), dtype=bool)
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], np.maximum.accumulate(ends[order])
    i = np.searchsorted(s, x, side="right") - 1
    return (i >= 0) & (e[np.maximum(i, 0)] >= x)


def _expand(starts: np.ndarray, lens: np.ndarray) -> np.ndarray:
    """Positions start .. start + len - 1 of every run, concatenated."""
    if len(lens) == 0:
        return np.zeros(0, dtype=np.int64)
    offs = np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
    return np.repeat(starts, lens) + offs


# query / reference consumed per CIGAR op (M I D N S H P = X)
_Q_ADV = np.array([1, 1, 0, 0, 1, 0, 0, 1, 1, 0], dtype=np.int64)
_R_ADV = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1, 0], dtype=np.int64)
_COMP = np.arange(256, dtype=np.uint8)
for _a, _b in zip(b"ACGTNacgtn", b"TGCANtgcan", strict=True):
    _COMP[_a] = _b


def _as_array(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode(), dtype=np.uint8)


def run_ava(
    reads: dict[str, SeqRecord],
    tmp_dir: Path,
    params: PhasingParams,
    threads: int,
    timeout: int,
) -> list[pysam.AlignedSegment]:
    """All-vs-all alignment with explicit mismatches, best alignment per pair."""
    fa = tmp_dir / "phasing_reads.fasta"
    with open(fa, "w") as f:
        SeqIO.write(list(reads.values()), f, "fasta")
    sam = tmp_dir / "phasing_ava.sam"
    # -D: skip the self diagonal; every other read is a secondary hit of the
    # query, so -N must allow one per read. --dual=no: each pair only once
    dual = " --dual=no" if params.align_pairs_once else ""
    cmd = (
        f"minimap2 -a {params.minimap_params} -D{dual} --eqx -t {threads} "
        f"--secondary=yes -N {len(reads)} -p 0.05 {fa} {fa}"
    )
    log.info(cmd)
    with open(sam, "w") as out:
        subprocess.run(
            shlex.split(cmd),
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=True,
            timeout=timeout,
        )
    best: dict[tuple[str, str], pysam.AlignedSegment] = {}
    with pysam.AlignmentFile(str(sam), "r", check_sq=False) as f:
        for a in f:
            if a.is_unmapped or a.is_supplementary or a.cigartuples is None:
                continue
            if a.query_name == a.reference_name:
                continue
            key = (a.query_name, a.reference_name)
            if params.align_pairs_once:
                key = min(key, key[::-1])
            if key not in best or a.get_tag("AS") > best[key].get_tag("AS"):
                best[key] = a
    return list(best.values())


class ReadArrays:
    """Read sequences as uint8 arrays, with cached reverse complements."""

    def __init__(self, seqs: dict[str, str]):
        self.fwd = {n: _as_array(s) for n, s in seqs.items()}
        self._rc: dict[str, np.ndarray] = {}

    def oriented(self, name: str, reverse: bool) -> np.ndarray:
        if not reverse:
            return self.fwd[name]
        if name not in self._rc:
            self._rc[name] = _COMP[self.fwd[name][::-1]]
        return self._rc[name]


_OPCODE = np.full(256, -1, dtype=np.int64)
for _k, _c in enumerate(b"MIDNSHP=XB"):
    _OPCODE[_c] = _k


def _cigar(aln: pysam.AlignedSegment):
    """ops, lengths and the oriented query / reference start of every op."""
    # parse the CIGAR string in numpy: much faster than cigartuples for the
    # hundreds of ops of a long read-to-read alignment
    cs = np.frombuffer(aln.cigarstring.encode(), dtype=np.uint8)
    is_op = cs >= 61  # letters and '='
    op_at = np.flatnonzero(is_op)
    digit_at = np.flatnonzero(~is_op)
    k = np.searchsorted(op_at, digit_at)
    lens = np.bincount(
        k,
        weights=(cs[digit_at] - 48) * 10.0 ** (op_at[k] - digit_at - 1),
        minlength=len(op_at),
    ).astype(np.int64)
    ops = _OPCODE[cs[op_at]]
    if (ops == 0).any():
        raise ValueError("parse_pair needs minimap2 --eqx alignments")
    q_adv, r_adv = lens * _Q_ADV[ops], lens * _R_ADV[ops]
    # a leading hard clip: oriented query coordinates start after it
    q0 = int(lens[0]) if ops[0] == 5 else 0
    qst = q0 + np.cumsum(q_adv) - q_adv
    rst = aln.reference_start + np.cumsum(r_adv) - r_adv
    return ops, lens, qst, rst


def _indels(ins_pos, ins_len, del_pos, del_len):
    pos = np.concatenate((ins_pos, del_pos))
    ln = np.concatenate((ins_len, -del_len))
    order = np.argsort(pos, kind="stable")
    return pos[order], ln[order]


def parse_pair(
    aln: pysam.AlignedSegment, seqs: "dict[str, str] | ReadArrays", _cig=None
) -> PairAln:
    """The alignment of q (aln's query) on t (aln's reference)."""
    arrs = seqs if isinstance(seqs, ReadArrays) else ReadArrays(seqs)
    # secondary records carry SEQ '*': take the read from `seqs`
    qa = arrs.oriented(aln.query_name, aln.is_reverse)
    ops, lens, qst, rst = _cig if _cig is not None else _cigar(aln)
    x, i, d = ops == 8, ops == 1, ops == 2
    small = (i | d) & (lens <= 10)
    n_x = int(lens[x].sum())
    n_eq = int(lens[ops == 7].sum())
    n_err = n_x + int(lens[small].sum())
    indel_pos, indel_len = _indels(rst[i], lens[i], rst[d], lens[d])
    return PairAln(
        q=aln.query_name,
        t=aln.reference_name,
        t_start=aln.reference_start,
        t_end=aln.reference_end,
        mism_pos=_expand(rst[x], lens[x]),
        mism_base=qa[_expand(qst[x], lens[x])],
        del_start=rst[d],
        del_len=lens[d],
        ins_after=rst[i] - 1,
        indel_pos=indel_pos,
        indel_len=indel_len,
        clip_left=int(lens[0]) if ops[0] in (4, 5) else 0,
        clip_right=int(lens[-1]) if ops[-1] in (4, 5) else 0,
        divergence=n_err / max(1, n_eq + n_err),
    )


def parse_pair_both(
    aln: pysam.AlignedSegment, seqs: "dict[str, str] | ReadArrays"
) -> tuple[PairAln, PairAln]:
    """The alignment of q on t, and the same alignment inverted: t on q.

    With ``--dual=no`` minimap2 aligns every read pair once; the inverted view
    makes q the target, as if t had been aligned to it: insertions become
    deletions and vice versa, mismatches carry t's base, and for a reverse
    alignment the coordinates are mirrored onto q's forward strand.
    """
    arrs = seqs if isinstance(seqs, ReadArrays) else ReadArrays(seqs)
    cig = _cigar(aln)
    fwd = parse_pair(aln, arrs, _cig=cig)
    ta = arrs.fwd[aln.reference_name]
    lq, lt = len(arrs.fwd[aln.query_name]), len(ta)
    rev = aln.is_reverse
    ops, lens, qst, rst = cig
    x, i, d = ops == 8, ops == 1, ops == 2
    qx = _expand(qst[x], lens[x])
    base = ta[_expand(rst[x], lens[x])]
    if rev:
        # oriented query position k is forward position lq - 1 - k
        mism_pos, mism_base = lq - 1 - qx, _COMP[base]
        del_start = lq - qst[i] - lens[i]  # q bases missing in t
        ins_after = lq - 1 - qst[d]  # t bases missing in q
    else:
        mism_pos, mism_base = qx, base
        del_start = qst[i]
        ins_after = qst[d] - 1
    indel_pos, indel_len = _indels(ins_after + 1, lens[d], del_start, lens[i])
    clip_l, clip_r = fwd.clip_left, fwd.clip_right
    q_start, q_end = clip_l, lq - clip_r
    r_start, r_end = aln.reference_start, aln.reference_end
    inv = PairAln(
        q=aln.reference_name,
        t=aln.query_name,
        t_start=lq - q_end if rev else q_start,
        t_end=lq - q_start if rev else q_end,
        mism_pos=mism_pos,
        mism_base=mism_base,
        del_start=del_start,
        del_len=lens[i],
        ins_after=ins_after,
        indel_pos=indel_pos,
        indel_len=indel_len,
        clip_left=lt - r_end if rev else r_start,
        clip_right=r_start if rev else lt - r_end,
        divergence=fwd.divergence,
    )
    return fwd, inv


def low_quality_reads(pairs: list[PairAln], params: PhasingParams) -> set[str]:
    """Reads whose median divergence to their partners is far above the pool's."""
    d: dict[str, list[float]] = defaultdict(list)
    for p in pairs:
        d[p.q].append(p.divergence)
        d[p.t].append(p.divergence)
    if not d:
        return set()
    div = {r: float(np.median(v)) for r, v in d.items()}
    med = float(np.median(list(div.values())))
    cut = max(params.lowq_factor * med, med + params.lowq_min_excess)
    return {r for r, v in div.items() if v > cut}


# --------------------------------------------------------------------------- #
# sites
# --------------------------------------------------------------------------- #
@dataclass
class Site:
    t: str
    pos: int
    kind: str  # "snv" or "sv"
    agree: list[str]  # reads with t's allele (t included)
    disagree: list[str]  # reads with the alternative allele
    depth: int


def _homopolymer_mask(seq: str, min_len: int, pad: int = 1) -> np.ndarray:
    n = len(seq)
    m = np.zeros(n, dtype=bool)
    if n == 0:
        return m
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    bounds = np.concatenate(([0], np.flatnonzero(a[1:] != a[:-1]) + 1, [n]))
    starts, ends = bounds[:-1], bounds[1:]
    long = ends - starts >= min_len
    for i, j in zip(starts[long], ends[long], strict=True):
        m[max(0, i - pad) : min(n, j + pad)] = True
    return m


def snv_sites(
    by_t: dict[str, list[PairAln]], seqs: dict[str, str], params: PhasingParams
) -> list[Site]:
    """Columns of every target read t where >= min_alt reads share one
    alternative base and >= min_ref reads (t included) carry t's base."""
    sites: list[Site] = []
    w = params.indel_window
    for t, plist in by_t.items():
        tlen = len(seqs[t])
        # columns where >= min_alt reads mismatch at all; most mismatch
        # columns are single sequencing errors
        all_pos = np.concatenate([p.mism_pos for p in plist])
        if len(all_pos) == 0:
            continue
        cand = np.flatnonzero(np.bincount(all_pos, minlength=tlen) >= params.min_alt)
        if len(cand) == 0:
            continue
        cand = cand[~_homopolymer_mask(seqs[t], params.homopolymer_len)[cand]]
        if len(cand) == 0:
            continue
        # B[i, c]: read i's mismatching base at candidate column c (0: none)
        col = np.full(tlen, -1, dtype=np.int64)
        col[cand] = np.arange(len(cand))
        B = np.zeros((len(plist), len(cand)), dtype=np.uint8)
        for i, p in enumerate(plist):
            c = col[p.mism_pos]
            sel = c >= 0
            B[i, c[sel]] = p.mism_base[sel]
        # the alternative base: the most frequent one, ties to the one seen
        # first in pair order (as collections.Counter.most_common)
        bases = np.unique(B[B > 0])
        counts = np.stack([(B == b).sum(axis=0) for b in bases])
        first = np.stack(
            [
                np.where((B == b).any(axis=0), (B == b).argmax(axis=0), len(plist))
                for b in bases
            ]
        )
        n_alt = counts.max(axis=0)
        best = np.where(counts == n_alt, first, len(plist) + 1).argmin(axis=0)
        alt = bases[best]
        keep = n_alt >= params.min_alt
        cand, B, alt = cand[keep], B[:, keep], alt[keep]
        if len(cand) == 0:
            continue
        # reads that cover a column outside their indels
        usable = np.stack(
            [
                (p.t_start <= cand) & (cand < p.t_end) & ~p.blocked_at(cand, w)
                for p in plist
            ]
        )
        agree = usable & (B == 0)
        disagree = usable & (B == alt[None, :])
        n_agree = agree.sum(axis=0) + 1  # t itself
        n_dis = disagree.sum(axis=0)
        depth = n_agree + n_dis
        ok = (n_agree >= params.min_ref) & (n_dis >= params.min_alt)
        ok &= np.minimum(n_agree, n_dis) / depth >= params.min_alt_frac
        qs = [p.q for p in plist]
        for c in np.flatnonzero(ok):
            sites.append(
                Site(
                    t,
                    int(cand[c]),
                    "snv",
                    [t] + [qs[i] for i in np.flatnonzero(agree[:, c])],
                    [qs[i] for i in np.flatnonzero(disagree[:, c])],
                    int(depth[c]),
                )
            )
    return sites


def recurrent_sites(
    sites: list[Site], min_support: int, min_conc: float = 0.9, min_dist: int = 20
) -> list[Site]:
    """Sites whose read split recurs at >= min_support other positions of the
    same target read."""
    by_t: dict[str, list[Site]] = defaultdict(list)
    for s in sites:
        by_t[s.t].append(s)
    keep_ids: set[int] = set()
    for sl in by_t.values():
        if len(sl) < 2:
            continue
        # +1 / -1: read has t's / the alternative allele at the site
        idx = {
            r: k for k, r in enumerate({r for s in sl for r in s.agree + s.disagree})
        }
        pos_m = np.zeros((len(sl), len(idx)), dtype=np.int32)
        neg_m = np.zeros((len(sl), len(idx)), dtype=np.int32)
        for i, s in enumerate(sl):
            pos_m[i, [idx[r] for r in s.agree]] = 1
            neg_m[i, [idx[r] for r in s.disagree]] = 1
        seen = pos_m + neg_m
        shared = seen @ seen.T
        same = pos_m @ pos_m.T + neg_m @ neg_m.T
        p = np.array([s.pos for s in sl])
        far = np.abs(p[:, None] - p[None, :]) >= min_dist
        np.fill_diagonal(far, False)
        with np.errstate(divide="ignore", invalid="ignore"):
            conc = same / shared
        support = (far & (shared >= 4) & (conc >= min_conc)).sum(axis=1)
        keep_ids.update(
            id(s) for s, n in zip(sl, support, strict=True) if n >= min_support
        )
    return [s for s in sites if id(s) in keep_ids]


def sv_sites(
    by_t: dict[str, list[PairAln]], seqs: dict[str, str], params: PhasingParams
) -> list[Site]:
    """SV-scale sites: net indel >= min_sv over a window, or a clip inside it."""
    sites: list[Site] = []
    for t, plist in by_t.items():
        tlen = len(seqs[t])
        ev: list[int] = []
        for p in plist:
            ev.extend(p.indel_pos[np.abs(p.indel_len) >= params.min_sv].tolist())
            if p.clip_left >= params.min_clip and p.t_start > params.min_clip:
                ev.append(p.t_start)
            if p.clip_right >= params.min_clip and p.t_end < tlen - params.min_clip:
                ev.append(p.t_end)
        if not ev:
            continue
        ev.sort()
        wins = [[ev[0], ev[0]]]
        for x in ev[1:]:
            if x - wins[-1][1] <= params.sv_merge_dist:
                wins[-1][1] = x
            else:
                wins.append([x, x])
        # net indel length of each read over [lo, hi] from prefix sums
        csum = [np.concatenate(([0], np.cumsum(p.indel_len))) for p in plist]
        for a, b in wins:
            lo, hi = a - params.sv_pad, b + params.sv_pad
            if lo < 0 or hi > tlen:
                continue
            agree, disagree = [t], []
            for p, cs in zip(plist, csum, strict=True):
                if (p.clip_left >= params.min_clip and lo <= p.t_start <= hi) or (
                    p.clip_right >= params.min_clip and lo <= p.t_end <= hi
                ):
                    disagree.append(p.q)
                    continue
                if not (p.t_start <= lo and p.t_end >= hi):
                    continue
                net = int(
                    cs[np.searchsorted(p.indel_pos, hi, side="right")]
                    - cs[np.searchsorted(p.indel_pos, lo, side="left")]
                )
                if abs(net) >= params.min_sv:
                    disagree.append(p.q)
                elif abs(net) < params.min_sv / 2:
                    agree.append(p.q)
            if len(disagree) < 2 or len(agree) < 2:
                continue
            sites.append(
                Site(t, (a + b) // 2, "sv", agree, disagree, len(agree) + len(disagree))
            )
    return sites


# --------------------------------------------------------------------------- #
# clustering
# --------------------------------------------------------------------------- #
def pair_weights(sites: list[Site], names: list[str], error_rate: float) -> np.ndarray:
    """Signed log-likelihood-ratio weights (same vs different haplotype)."""
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    A = np.zeros((n, n))
    D = np.zeros((n, n))
    for s in sites:
        ti = idx[s.t]
        # agreement is informative only at a balanced split: where 2 of 20
        # reads differ, every other read "agrees" with t whatever its haplotype
        w_agree = min(1.0, 2.0 * min(len(s.agree), len(s.disagree)) / max(1, s.depth))
        for q in s.agree:
            if q != s.t:
                A[ti, idx[q]] += w_agree
        for q in s.disagree:
            D[ti, idx[q]] += 1.0
    A = A + A.T
    D = D + D.T
    return np.log((1 - error_rate) / error_rate) * (A - D)


def correlation_cluster(W: np.ndarray) -> list[int]:
    """Greedy agglomeration: merge the two clusters with the largest summed
    inter-cluster weight while it is positive."""
    n = W.shape[0]
    clusters = [[i] for i in range(n)]
    C = W.copy()
    np.fill_diagonal(C, 0)
    active = list(range(n))
    while len(active) > 1:
        sub = C[np.ix_(active, active)]
        np.fill_diagonal(sub, -np.inf)
        k = int(np.argmax(sub))
        a, b = divmod(k, len(active))
        if sub[a, b] <= 0:
            break
        ia, ib = active[a], active[b]
        clusters[ia].extend(clusters[ib])
        clusters[ib] = []
        C[ia, :] += C[ib, :]
        C[:, ia] += C[:, ib]
        active.remove(ib)
    labels = [-1] * n
    for lab, members in enumerate(c for c in clusters if c):
        for i in members:
            labels[i] = lab
    return labels


def discriminating_positions(
    sites: list[Site],
    a_set: set[str],
    b_set: set[str],
    min_obs: int = 2,
    min_major: float = 0.8,
) -> int:
    """Distinct variant positions at which clusters a and b are each
    near-unanimous with different alleles (max over targets: every variant is
    seen once per target read)."""
    per_t: Counter = Counter()
    for s in sites:
        ag, di = set(s.agree), set(s.disagree)
        a0, a1 = len(ag & a_set), len(di & a_set)
        b0, b1 = len(ag & b_set), len(di & b_set)
        if a0 + a1 < min_obs or b0 + b1 < min_obs:
            continue
        fa, fb = a0 / (a0 + a1), b0 / (b0 + b1)
        if (fa >= min_major and fb <= 1 - min_major) or (
            fb >= min_major and fa <= 1 - min_major
        ):
            per_t[s.t] += 1
    return max(per_t.values(), default=0)


def refine_clusters(
    labels: list[int],
    names: list[str],
    sites: list[Site],
    W: np.ndarray,
    params: PhasingParams,
) -> list[int]:
    lab = np.array(labels)
    while True:
        cnt = Counter(lab[lab >= 0].tolist())
        big = [lb for lb, c in cnt.most_common() if c >= params.min_group]
        if len(big) < 2:
            break
        sets = {lb: {names[i] for i in np.where(lab == lb)[0]} for lb in big}
        best = None
        for i, a in enumerate(big):
            for b in big[i + 1 :]:
                d = discriminating_positions(sites, sets[a], sets[b])
                if best is None or d < best[0]:
                    best = (d, a, b)
        if best[0] >= params.min_discriminating:
            break
        lab[lab == best[2]] = best[1]
    cnt = Counter(lab[lab >= 0].tolist())
    big = [lb for lb, c in cnt.most_common() if c >= params.min_group]
    out = lab.copy()
    for i, lb in enumerate(lab):
        if lb in big:
            continue
        out[i] = -1
        if not big:
            continue
        scores = np.array([W[i, lab == b].sum() for b in big])
        order = np.argsort(scores)[::-1]
        top = scores[order[0]]
        second = scores[order[1]] if len(big) > 1 else 0.0
        if top > 0 and top > second:
            out[i] = big[order[0]]
    cnt = Counter(out[out >= 0].tolist())
    remap = {lb: k for k, (lb, _) in enumerate(cnt.most_common())}
    return [remap.get(lb, -1) for lb in out.tolist()]


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #
def phase_reads(
    reads: dict[str, SeqRecord],
    params: PhasingParams | None = None,
    threads: int = 1,
    timeout: int = 120,
    tmp_dir_path: Path | str | None = None,
) -> PhasingResult:
    """Phase `reads` (cut to the candidate region +- a phasing flank).

    status:
      "phased"         >= 2 alleles
      "single"         the reads form one allele
      "no_information" no read pair could be compared (no informative sites)
      "failed"         the all-vs-all alignment failed or timed out
    """
    params = params or PhasingParams()
    if len(reads) < 2 * params.min_group:
        return PhasingResult(status="no_information", unassigned=sorted(reads))
    all_reads = reads
    if params.max_reads and len(reads) > params.max_reads:
        keep = sorted(reads, key=lambda n: (-len(reads[n].seq), n))[: params.max_reads]
        reads = {n: reads[n] for n in sorted(keep)}
    seqs = {n: str(r.seq) for n, r in reads.items()}
    try:
        with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tmp:
            alns = run_ava(reads, Path(tmp), params, threads=threads, timeout=timeout)
            if params.align_pairs_once:
                arrs = ReadArrays(seqs)
                pairs = [p for a in alns for p in parse_pair_both(a, arrs)]
            else:
                arrs = ReadArrays(seqs)
                pairs = [parse_pair(a, arrs) for a in alns]
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
        log.warning(f"read phasing: all-vs-all alignment failed: {e}")
        if isinstance(e, subprocess.TimeoutExpired):
            tool_timeouts.record("phasing all-vs-all")
        return PhasingResult(status="failed", unassigned=sorted(all_reads))

    lowq = low_quality_reads(pairs, params)
    pairs = [p for p in pairs if p.q not in lowq and p.t not in lowq]
    by_t: dict[str, list[PairAln]] = defaultdict(list)
    for p in pairs:
        by_t[p.t].append(p)
    s_snv = snv_sites(by_t, seqs, params)
    if params.recurrence > 0:
        s_snv = recurrent_sites(s_snv, min_support=params.recurrence)
    s_sv = sv_sites(by_t, seqs, params) if params.use_sv_sites else []
    sites = s_snv + s_sv

    names = sorted(r for r in reads if r not in lowq)
    W = pair_weights(sites, names, params.error_rate)
    labels = refine_clusters(correlation_cluster(W), names, sites, W, params)

    groups = {n: lb for n, lb in zip(names, labels, strict=True) if lb >= 0}
    n_alleles = len(set(groups.values()))
    if n_alleles > params.max_alleles:
        # keep the largest; the reads of the others are left to be re-added
        keep = {
            lb for lb, _ in Counter(groups.values()).most_common(params.max_alleles)
        }
        groups = {n: lb for n, lb in groups.items() if lb in keep}
        n_alleles = params.max_alleles
    status = (
        "phased" if n_alleles >= 2 else "single" if n_alleles == 1 else "no_information"
    )
    res = PhasingResult(
        status=status,
        groups=groups,
        n_alleles=n_alleles,
        unassigned=sorted(n for n in all_reads if n not in groups and n not in lowq),
        low_quality=sorted(lowq),
        n_snv_sites=len(s_snv),
        n_sv_sites=len(s_sv),
    )
    log.info(
        f"read phasing: status={status} alleles={n_alleles} "
        f"sizes={sorted(Counter(groups.values()).values(), reverse=True)} "
        f"unassigned={len(res.unassigned)} low_quality={len(lowq)} "
        f"snv_sites={len(s_snv)} sv_sites={len(s_sv)}"
    )
    return res
