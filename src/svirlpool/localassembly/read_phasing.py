"""Experimental: phase the reads of a candidate-region container from their
all-vs-all alignments, and derive the number of alleles from the phasing.

The legacy clustering in ``consensus`` compares reads by their SV-scale
signals only and needs the local copy number as the number of clusters.  Reads
of one haplotype, however, also share every small variant -- above all SNVs --
over their whole length, and the reads are much longer than the candidate
region.  This module aligns reads cut to the candidate region +- a phasing
flank all-vs-all (minimap2 ``--eqx``, both directions, so every read is the
target of all others) and then

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

log = logging.getLogger(__name__)

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class PhasingParams:
    # minimap2: the ava-ont preset without -X (all pairs in both directions)
    minimap_params: str = "-k15 -w5 -e0 -m100 -r2k"
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
    lowq_factor: float = 2.5
    lowq_min_excess: float = 0.02
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
    q: str
    t: str
    t_start: int
    t_end: int
    mism: dict[int, str]  # t position -> q base at X columns
    dels: set[int]  # t positions deleted in q
    ins_after: set[int]  # t positions followed by an insertion in q
    indels: list[tuple[int, int]]  # (t position, +ins / -del length)
    clip_left: int
    clip_right: int
    divergence: float  # (X + indel bases of indels <= 10 bp) per aligned base


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
    # query, so -N must allow one per read
    cmd = (
        f"minimap2 -a {params.minimap_params} -D --eqx -t {threads} "
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
            if key not in best or a.get_tag("AS") > best[key].get_tag("AS"):
                best[key] = a
    return list(best.values())


def parse_pair(aln: pysam.AlignedSegment, seqs: dict[str, str]) -> PairAln:
    # secondary records carry SEQ '*': take the read from `seqs`
    qs = seqs[aln.query_name]
    if aln.is_reverse:
        qs = qs.translate(_RC)[::-1]
    ct = aln.cigartuples
    if ct[0][0] == 5:
        qs = qs[ct[0][1] :]
    rpos, qpos = aln.reference_start, 0
    mism: dict[int, str] = {}
    dels: set[int] = set()
    ins_after: set[int] = set()
    indels: list[tuple[int, int]] = []
    n_eq = n_small = 0
    for op, ln in ct:
        if op == 7:  # =
            rpos += ln
            qpos += ln
            n_eq += ln
        elif op == 8:  # X
            for k in range(ln):
                mism[rpos + k] = qs[qpos + k]
            rpos += ln
            qpos += ln
        elif op == 1:
            ins_after.add(rpos - 1)
            indels.append((rpos, ln))
            qpos += ln
            n_small += ln if ln <= 10 else 0
        elif op == 2:
            dels.update(range(rpos, rpos + ln))
            indels.append((rpos, -ln))
            rpos += ln
            n_small += ln if ln <= 10 else 0
        elif op == 4:
            qpos += ln
        elif op == 0:
            raise ValueError("parse_pair needs minimap2 --eqx alignments")
    n_err = len(mism) + n_small
    return PairAln(
        q=aln.query_name,
        t=aln.reference_name,
        t_start=aln.reference_start,
        t_end=aln.reference_end,
        mism=mism,
        dels=dels,
        ins_after=ins_after,
        indels=indels,
        clip_left=ct[0][1] if ct[0][0] in (4, 5) else 0,
        clip_right=ct[-1][1] if ct[-1][0] in (4, 5) else 0,
        divergence=n_err / max(1, n_eq + n_err),
    )


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
    i = 0
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            m[max(0, i - pad) : min(n, j + pad)] = True
        i = j
    return m


def snv_sites(
    by_t: dict[str, list[PairAln]], seqs: dict[str, str], params: PhasingParams
) -> list[Site]:
    sites: list[Site] = []
    w = params.indel_window
    for t, plist in by_t.items():
        hp = _homopolymer_mask(seqs[t], params.homopolymer_len)
        cand: dict[int, Counter] = defaultdict(Counter)
        for p in plist:
            for pos, b in p.mism.items():
                cand[pos][b] += 1
        for pos in sorted(cand):
            if hp[pos]:
                continue
            alt, n_alt = cand[pos].most_common(1)[0]
            if n_alt < params.min_alt:
                continue
            agree, disagree = [t], []
            for p in plist:
                if not (p.t_start <= pos < p.t_end):
                    continue
                if any((pos + d) in p.dels for d in range(-w, w + 1)):
                    continue
                if any((pos + d) in p.ins_after for d in range(-w - 1, w + 1)):
                    continue
                b = p.mism.get(pos)
                if b is None:
                    agree.append(p.q)
                elif b == alt:
                    disagree.append(p.q)
            depth = len(agree) + len(disagree)
            if len(agree) < params.min_ref or len(disagree) < params.min_alt:
                continue
            if min(len(agree), len(disagree)) / depth < params.min_alt_frac:
                continue
            sites.append(Site(t, pos, "snv", agree, disagree, depth))
    return sites


def _concordance(a: Site, b: Site, min_shared: int = 4) -> float | None:
    sa = dict.fromkeys(a.agree, 1) | dict.fromkeys(a.disagree, -1)
    sb = dict.fromkeys(b.agree, 1) | dict.fromkeys(b.disagree, -1)
    shared = sa.keys() & sb.keys()
    if len(shared) < min_shared:
        return None
    return sum(sa[x] == sb[x] for x in shared) / len(shared)


def recurrent_sites(
    sites: list[Site], min_support: int, min_conc: float = 0.9, min_dist: int = 20
) -> list[Site]:
    """Sites whose read split recurs at >= min_support other positions of the
    same target read."""
    by_t: dict[str, list[Site]] = defaultdict(list)
    for s in sites:
        by_t[s.t].append(s)
    keep: list[Site] = []
    for sl in by_t.values():
        for i, a in enumerate(sl):
            sup = 0
            for j, b in enumerate(sl):
                if i == j or abs(a.pos - b.pos) < min_dist:
                    continue
                c = _concordance(a, b)
                if c is not None and c >= min_conc:
                    sup += 1
                    if sup >= min_support:
                        break
            if sup >= min_support:
                keep.append(a)
    return keep


def sv_sites(
    by_t: dict[str, list[PairAln]], seqs: dict[str, str], params: PhasingParams
) -> list[Site]:
    """SV-scale sites: net indel >= min_sv over a window, or a clip inside it."""
    sites: list[Site] = []
    for t, plist in by_t.items():
        tlen = len(seqs[t])
        ev: list[int] = []
        for p in plist:
            ev.extend(pos for pos, sz in p.indels if abs(sz) >= params.min_sv)
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
        for a, b in wins:
            lo, hi = a - params.sv_pad, b + params.sv_pad
            if lo < 0 or hi > tlen:
                continue
            agree, disagree = [t], []
            for p in plist:
                if (p.clip_left >= params.min_clip and lo <= p.t_start <= hi) or (
                    p.clip_right >= params.min_clip and lo <= p.t_end <= hi
                ):
                    disagree.append(p.q)
                    continue
                if not (p.t_start <= lo and p.t_end >= hi):
                    continue
                net = sum(sz for pos, sz in p.indels if lo <= pos <= hi)
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
    seqs = {n: str(r.seq) for n, r in reads.items()}
    try:
        with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tmp:
            alns = run_ava(reads, Path(tmp), params, threads=threads, timeout=timeout)
            pairs = [parse_pair(a, seqs) for a in alns]
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
        log.warning(f"read phasing: all-vs-all alignment failed: {e}")
        return PhasingResult(status="failed", unassigned=sorted(reads))

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
        unassigned=sorted(n for n in names if n not in groups),
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
