"""Prototype: phase reads from the all-vs-all alignments of one container.

Informative sites are found per target read t from the alignments q->t
(minimap2 --eqx): a column of t where >= min_alt other reads share one
alternative base and >= min_ref reads (t included) carry t's base.  Such a
column is a variant between haplotypes (sequencing errors do not recur on the
same column with the same base), so every read q aligned over it either agrees
with t or disagrees with t.  Summed over all targets this gives per read pair
(i, j): agree[i, j], disagree[i, j].
"""

from __future__ import annotations

import json
import sqlite3
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus as cmod
from svirlpool.localassembly import read_cache as read_cache_mod

WORK = Path(__file__).parent / "data" / "flank200"
BAM = Path("/home/mayv_c/biodata/local/alignments/HG002.minimap2.20x.softclipped.bam")


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def load_container(crID: int, db: Path = WORK / "crs_containers.db") -> dict:
    conts = cmod.load_crs_containers_from_db(db, [crID])
    cont = conts[crID]
    return {cr.crID: cr for cr in cont["crs"]}


def get_reads(
    crs_dict: dict, flank: int, buffer_clipped_length: int = 5000
) -> tuple[dict[str, SeqRecord], dict[str, SeqRecord], dict]:
    """Cut reads the way process_consensus_container does, with CR_CUT_FLANK=flank.

    Returns (cutreads, full read records, alignments per crID).
    """
    old = cmod.CR_CUT_FLANK
    cmod.CR_CUT_FLANK = flank
    try:
        cache = read_cache_mod.ReadSequenceCache(path_alignments=BAM)
        alns, recs = {}, {}
        for cr in crs_dict.values():
            a, s = cache.fetch_for_cr(cr)
            alns[cr.crID] = a
            recs.update(s)
        ivs = cmod.get_read_alignment_intervals_in_cr(
            crs=list(crs_dict.values()),
            dict_alignments=alns,
            buffer_clipped_length=buffer_clipped_length,
        )
        mx = cmod.get_max_extents_of_read_alignments_on_cr(ivs)
        cut = cmod.trim_reads(dict_alignments=alns, intervals=mx, read_records=recs)
    finally:
        cmod.CR_CUT_FLANK = old
    return cut, recs, alns


# --------------------------------------------------------------------------- #
# AVA
# --------------------------------------------------------------------------- #
def run_ava(
    reads: dict[str, SeqRecord], tmp: Path, params: str = "-k15 -w5 -e0 -m100 -r2k", extra: str = "",
    timeout: int = 120,
) -> list[pysam.AlignedSegment]:
    fa = tmp / "reads.fa"
    with open(fa, "w") as f:
        SeqIO.write(list(reads.values()), f, "fasta")
    sam = tmp / "ava.sam"
    # the ava-ont preset without -X: every pair in both directions, so each
    # target read gets the complete pileup of the others
    cmd = (f"minimap2 -a {params} -D --eqx -t 2 --secondary=yes -N 1000 -p 0.05 "
           f"{extra} {fa} {fa}")
    with open(sam, "w") as out:
        subprocess.run(cmd.split(), stdout=out, stderr=subprocess.DEVNULL,
                       check=True, timeout=timeout)
    # every other read is a secondary hit of the query: keep the best
    # alignment per (query, target) pair
    best: dict[tuple[str, str], pysam.AlignedSegment] = {}
    with pysam.AlignmentFile(str(sam), "r", check_sq=False) as f:
        for a in f:
            if a.is_unmapped or a.is_supplementary or a.query_name == a.reference_name:
                continue
            key = (a.query_name, a.reference_name)
            if key not in best or a.get_tag("AS") > best[key].get_tag("AS"):
                best[key] = a
    return list(best.values())


# --------------------------------------------------------------------------- #
# informative sites
# --------------------------------------------------------------------------- #
@dataclass
class PairAln:
    q: str
    t: str
    t_start: int
    t_end: int
    # t position -> q base, for X columns; t positions deleted in q
    mism: dict[int, str]
    dels: set[int]
    ins_after: dict[int, int]  # t position -> inserted length after it (>0)
    indels: list[tuple[int, int]] = field(default_factory=list)  # (t pos, +ins/-del)
    clip_left: int = 0  # unaligned query bases before / after the alignment
    clip_right: int = 0
    divergence: float = 0.0  # mismatches + small indel bases per aligned base


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_pair(aln: pysam.AlignedSegment, seqs: dict[str, str]) -> PairAln:
    # secondary records carry SEQ '*': take the read from `seqs`, in the
    # orientation of the alignment
    qs = seqs[aln.query_name]
    if aln.is_reverse:
        qs = qs.translate(_RC)[::-1]
    if aln.cigartuples[0][0] == 5:  # hard clip: SEQ excludes it
        qs = qs[aln.cigartuples[0][1]:]
    rpos, qpos = aln.reference_start, 0
    mism, dels, ins, indels = {}, set(), {}, []
    for op, ln in aln.cigartuples:
        if op == 7:  # =
            rpos += ln
            qpos += ln
        elif op == 8:  # X
            for k in range(ln):
                mism[rpos + k] = qs[qpos + k]
            rpos += ln
            qpos += ln
        elif op == 0:
            raise ValueError("run minimap2 with --eqx")
        elif op == 1:
            ins[rpos - 1] = ln
            indels.append((rpos, ln))
            qpos += ln
        elif op == 2:
            dels.update(range(rpos, rpos + ln))
            indels.append((rpos, -ln))
            rpos += ln
        elif op == 4:
            qpos += ln
    ct = aln.cigartuples
    n_eq = sum(ln for op, ln in ct if op == 7)
    n_err = len(mism) + sum(ln for op, ln in ct if op in (1, 2) and ln <= 10)
    cl = ct[0][1] if ct[0][0] in (4, 5) else 0
    cr = ct[-1][1] if ct[-1][0] in (4, 5) else 0
    return PairAln(aln.query_name, aln.reference_name, aln.reference_start,
                   aln.reference_end, mism, dels, ins, indels, cl, cr,
                   n_err / max(1, n_eq + n_err))


def homopolymer_mask(seq: str, min_len: int = 4, pad: int = 1) -> np.ndarray:
    """True at positions inside or within `pad` of a homopolymer >= min_len."""
    n = len(seq)
    m = np.zeros(n, dtype=bool)
    i = 0
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            m[max(0, i - pad):min(n, j + pad)] = True
        i = j
    return m


@dataclass
class Site:
    t: str
    pos: int
    ref_base: str
    alt_base: str
    n_ref: int  # incl. t
    n_alt: int
    depth: int
    agree: list[str] = field(default_factory=list)
    disagree: list[str] = field(default_factory=list)


def read_divergence(pairs: list[PairAln]) -> dict[str, float]:
    """Median divergence of each read to its AVA partners."""
    d: dict[str, list[float]] = defaultdict(list)
    for p in pairs:
        d[p.q].append(p.divergence)
        d[p.t].append(p.divergence)
    return {r: float(np.median(v)) for r, v in d.items()}


def low_quality_reads(pairs: list[PairAln], factor: float = 2.5, min_excess: float = 0.02) -> set[str]:
    """Reads whose divergence to their partners is far above the pool's.

    Their errors are frequent enough that a few of them share a wrong base at
    many columns, which the recurrence filter mistakes for a haplotype.
    """
    div = read_divergence(pairs)
    if not div:
        return set()
    med = float(np.median(list(div.values())))
    return {r for r, v in div.items() if v > max(factor * med, med + min_excess)}


def informative_sites(
    pairs: list[PairAln],
    seqs: dict[str, str],
    min_alt: int = 2,
    min_ref: int = 2,
    min_alt_frac: float = 0.15,
    hp_len: int = 4,
    indel_window: int = 2,
) -> list[Site]:
    by_t: dict[str, list[PairAln]] = defaultdict(list)
    for p in pairs:
        by_t[p.t].append(p)
    sites: list[Site] = []
    for t, plist in by_t.items():
        tseq = seqs[t]
        hp = homopolymer_mask(tseq, hp_len)
        # candidate columns: any X
        cand: dict[int, Counter] = defaultdict(Counter)
        for p in plist:
            for pos, b in p.mism.items():
                cand[pos][b] += 1
        for pos, cnt in cand.items():
            if hp[pos]:
                continue
            alt, n_alt = cnt.most_common(1)[0]
            if n_alt < min_alt:
                continue
            agree, disagree, depth = [t], [], 1
            for p in plist:
                if not (p.t_start <= pos < p.t_end):
                    continue
                # skip reads with an indel right at/near the column
                if any((pos + d) in p.dels for d in range(-indel_window, indel_window + 1)):
                    continue
                if any((pos + d) in p.ins_after for d in range(-indel_window - 1, indel_window + 1)):
                    continue
                depth += 1
                b = p.mism.get(pos)
                if b is None:
                    agree.append(p.q)
                elif b == alt:
                    disagree.append(p.q)
            n_ref = len(agree)
            if n_ref < min_ref or len(disagree) < min_alt:
                continue
            if len(disagree) / depth < min_alt_frac or n_ref / depth < min_alt_frac:
                continue
            sites.append(Site(t, pos, tseq[pos], alt, n_ref, len(disagree), depth,
                              agree, disagree))
    return sites


def indel_sites(
    pairs: list[PairAln],
    seqs: dict[str, str],
    min_sv: int = 20,
    merge_dist: int = 100,
    pad: int = 50,
    min_alt: int = 2,
    min_ref: int = 2,
    min_clip: int = 100,
) -> list[Site]:
    """SV-scale variant sites on each target read.

    A window of t where some query has an indel >= min_sv (or an interior
    clip >= min_clip) is a site.  Queries spanning the window (+-pad) with a
    net indel >= min_sv disagree with t, those with |net| < min_sv/2 agree;
    queries whose alignment is clipped (>= min_clip) inside the window
    disagree (large SV / breakend).  Net indel over the window makes the
    call robust against fragmented indels in repeats.
    """
    by_t: dict[str, list[PairAln]] = defaultdict(list)
    for p in pairs:
        by_t[p.t].append(p)
    sites: list[Site] = []
    for t, plist in by_t.items():
        tlen = len(seqs[t])
        ev: list[int] = []
        for p in plist:
            ev.extend(pos for pos, sz in p.indels if abs(sz) >= min_sv)
            if p.clip_left >= min_clip and p.t_start > min_clip:
                ev.append(p.t_start)
            if p.clip_right >= min_clip and p.t_end < tlen - min_clip:
                ev.append(p.t_end)
        if not ev:
            continue
        ev.sort()
        wins = [[ev[0], ev[0]]]
        for x in ev[1:]:
            if x - wins[-1][1] <= merge_dist:
                wins[-1][1] = x
            else:
                wins.append([x, x])
        for a, b in wins:
            lo, hi = a - pad, b + pad
            if lo < 0 or hi > tlen:
                continue
            agree, disagree, depth = [t], [], 1
            for p in plist:
                clipped_in = (
                    (p.clip_left >= min_clip and lo <= p.t_start <= hi)
                    or (p.clip_right >= min_clip and lo <= p.t_end <= hi)
                )
                if clipped_in:
                    disagree.append(p.q)
                    depth += 1
                    continue
                if not (p.t_start <= lo and p.t_end >= hi):
                    continue
                net = sum(sz for pos, sz in p.indels if lo <= pos <= hi)
                depth += 1
                if abs(net) >= min_sv:
                    disagree.append(p.q)
                elif abs(net) < min_sv / 2:
                    agree.append(p.q)
            if len(disagree) < min_alt or len(agree) < min_ref:
                continue
            sites.append(Site(t, (a + b) // 2, "sv", "sv", len(agree),
                              len(disagree), depth, agree, disagree))
    return sites


def site_concordance(a: Site, b: Site, min_shared: int = 4) -> float | None:
    """Fraction of shared reads on the same side in two sites of one target."""
    sa = {x: 1 for x in a.agree} | {x: -1 for x in a.disagree}
    sb = {x: 1 for x in b.agree} | {x: -1 for x in b.disagree}
    shared = sa.keys() & sb.keys()
    if len(shared) < min_shared:
        return None
    return sum(sa[x] == sb[x] for x in shared) / len(shared)


def recurrent_sites(
    sites: list[Site], min_support: int = 1, min_conc: float = 0.9, min_dist: int = 20
) -> list[Site]:
    """Keep sites whose read split recurs at >= min_support other positions of
    the same target read (a haplotype split is shared by all its het SNVs; the
    read subsets of coinciding sequencing errors are random)."""
    by_t: dict[str, list[Site]] = defaultdict(list)
    for s in sites:
        by_t[s.t].append(s)
    keep = []
    for t, sl in by_t.items():
        for i, a in enumerate(sl):
            sup = 0
            for j, b in enumerate(sl):
                if i == j or abs(a.pos - b.pos) < min_dist:
                    continue
                c = site_concordance(a, b)
                if c is not None and c >= min_conc:
                    sup += 1
                    if sup >= min_support:
                        break
            if sup >= min_support:
                keep.append(a)
    return keep


def pair_counts(sites: list[Site], names: list[str], mode: str = "target", weight_agree: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """agree / disagree counts per read pair.

    mode="target": only pairs (q, t) at t's sites (hifiasm-like).
    mode="all": every pair of reads observed at a site (co-occurrence).
    """
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    A = np.zeros((n, n))
    D = np.zeros((n, n))
    for s in sites:
        ti = idx[s.t]
        # agreement is only informative at a balanced split: at a site where
        # 2 of 20 reads differ, every other read "agrees" with t whatever its
        # haplotype
        w_agree = min(1.0, 2.0 * min(s.n_alt, s.n_ref) / max(1, s.depth)) if weight_agree else 1.0
        if mode == "target":
            for q in s.agree:
                if q != s.t:
                    A[ti, idx[q]] += w_agree
            for q in s.disagree:
                D[ti, idx[q]] += 1
        else:
            ag = [idx[x] for x in s.agree]
            di = [idx[x] for x in s.disagree]
            for grp in (ag, di):
                for a in grp:
                    for b in grp:
                        if a != b:
                            A[a, b] += 0.5
            for a in ag:
                for b in di:
                    D[a, b] += 0.5
                    D[b, a] += 0.5
    A = A + A.T
    D = D + D.T
    return A, D


# --------------------------------------------------------------------------- #
# clustering
# --------------------------------------------------------------------------- #
def phase_weights(A: np.ndarray, D: np.ndarray, err: float = 0.1) -> np.ndarray:
    """Log-likelihood ratio same vs different haplotype per pair.

    Per site observation: same hap -> agree with prob 1-err; different hap ->
    disagree with prob 1-err.
    """
    l = np.log((1 - err) / err)
    return l * (A - D)


def correlation_cluster(W: np.ndarray, min_merge: float = 0.0) -> list[int]:
    """Greedy average-linkage agglomeration on signed weights.

    Merge the two clusters with the largest summed inter-cluster weight while
    it is > min_merge.  Returns labels.
    """
    n = W.shape[0]
    clusters = [[i] for i in range(n)]
    S = W.copy()
    np.fill_diagonal(S, 0)
    # inter-cluster sums
    C = S.copy()
    active = list(range(n))
    while len(active) > 1:
        sub = C[np.ix_(active, active)]
        np.fill_diagonal(sub, -np.inf)
        k = np.argmax(sub)
        a, b = divmod(k, len(active))
        if sub[a, b] <= min_merge:
            break
        ia, ib = active[a], active[b]
        clusters[ia].extend(clusters[ib])
        clusters[ib] = []
        C[ia, :] += C[ib, :]
        C[:, ia] += C[:, ib]
        active.remove(ib)
    labels = [-1] * n
    for lab, ci in enumerate([c for c in clusters if c]):
        for i in ci:
            labels[i] = lab
    return labels


# --------------------------------------------------------------------------- #
# convenience
# --------------------------------------------------------------------------- #
def analyse(crID: int, flank: int, tmp: Path, **kw) -> dict:
    crs = load_container(crID)
    cut, recs, alns = get_reads(crs, flank)
    names = sorted(cut)
    seqs = {n: str(cut[n].seq) for n in names}
    tmp.mkdir(parents=True, exist_ok=True)
    ava = run_ava(cut, tmp)
    pairs = [parse_pair(a, seqs) for a in ava]
    lowq = low_quality_reads(pairs) if kw.pop("filter_lowq", True) else set()
    pairs = [p for p in pairs if p.q not in lowq and p.t not in lowq]
    recur = kw.pop("recur", 1)
    use_sv = kw.pop("use_sv", True)
    use_snv = kw.pop("use_snv", True)
    sv_kw = kw.pop("sv_kw", {})
    raw_sites = informative_sites(pairs, seqs, **kw)
    sites = recurrent_sites(raw_sites, min_support=recur) if recur > 0 else raw_sites
    if not use_snv:
        sites = []
    sv_sites = indel_sites(pairs, seqs, **sv_kw) if use_sv else []
    sites = sites + sv_sites
    A, D = pair_counts(sites, names, weight_agree=kw.pop("weight_agree", True))
    return dict(lowq=lowq, raw_sites=raw_sites, sv_sites=sv_sites, crs=crs, names=names, seqs=seqs, pairs=pairs, sites=sites, A=A, D=D,
                lengths={n: len(s) for n, s in seqs.items()})


# --------------------------------------------------------------------------- #
# cluster-level refinement: allele count from discriminating sites
# --------------------------------------------------------------------------- #
def discriminating_sites(
    sites: list[Site], A_set: set[str], B_set: set[str],
    min_obs: int = 2, min_major: float = 0.8,
) -> int:
    """Number of distinct variant positions at which clusters A and B are each
    near-unanimous and carry different alleles.

    The same variant is seen once per target read; the count is the maximum
    over targets of the discriminating positions on that target.
    """
    per_t: Counter = Counter()
    for s in sites:
        ag, di = set(s.agree), set(s.disagree)
        a0, a1 = len(ag & A_set), len(di & A_set)
        b0, b1 = len(ag & B_set), len(di & B_set)
        if a0 + a1 < min_obs or b0 + b1 < min_obs:
            continue
        fa, fb = a0 / (a0 + a1), b0 / (b0 + b1)
        if (fa >= min_major and fb <= 1 - min_major) or (fb >= min_major and fa <= 1 - min_major):
            per_t[s.t] += 1
    return max(per_t.values(), default=0)


def refine_clusters(
    labels: list[int], names: list[str], sites: list[Site], W: np.ndarray,
    min_group: int = 3, min_disc: int = 2, assign_margin: float = 0.0,
) -> tuple[list[int], dict]:
    """Merge clusters that are not separated by >= min_disc discriminating
    variant positions; then assign reads of small clusters to the big cluster
    they have positive net evidence for (else -1)."""
    lab = np.array(labels)
    meta: dict = {"merges": 0, "assigned": 0}
    while True:
        cnt = Counter(lab[lab >= 0])
        big = [l for l, c in cnt.most_common() if c >= min_group]
        if len(big) < 2:
            break
        sets = {l: {names[i] for i in np.where(lab == l)[0]} for l in big}
        best = None
        for i, a in enumerate(big):
            for b in big[i + 1:]:
                d = discriminating_sites(sites, sets[a], sets[b])
                if best is None or d < best[0]:
                    best = (d, a, b)
        if best[0] >= min_disc:
            break
        _, a, b = best
        lab[lab == b] = a
        meta["merges"] += 1
    cnt = Counter(lab[lab >= 0])
    big = [l for l, c in cnt.most_common() if c >= min_group]
    out = lab.copy()
    for i, l in enumerate(lab):
        if l in big:
            continue
        if not big:
            out[i] = -1
            continue
        scores = [W[i, lab == b].sum() for b in big]
        order = np.argsort(scores)[::-1]
        top = scores[order[0]]
        second = scores[order[1]] if len(big) > 1 else 0.0
        if top > 0 and top - second > assign_margin:
            out[i] = big[order[0]]
            meta["assigned"] += 1
        else:
            out[i] = -1
    # relabel 0..k-1 by size
    cnt = Counter(out[out >= 0])
    remap = {l: k for k, (l, _) in enumerate(cnt.most_common())}
    meta["n_alleles"] = len(remap)
    return [remap.get(l, -1) for l in out], meta
