"""Truth helpers shared by container_truth_light.py, copied verbatim from
container_truth.py (which is a script and runs on import)."""

import bisect
import collections

TRUTH_VCF = {
    "T2T": "/home/mayv_c/development/svp_improvements/resources/truthsets/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz",
    "V5": "/home/mayv_c/biodata/local/benchmarks/HG002_GRCh38_v5.0q_stvar.filtered.vcf.gz",
}
TRUTH_BENCH = {
    "T2T": "/home/mayv_c/development/svp_improvements/resources/truthsets/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed",
    "V5": "/home/mayv_c/development/svp_improvements/resources/truthsets/HG002_GRCh38_v5.0q_stvar.benchmark.bed",
}
TRF = "/home/mayv_c/biodata/local/references/GRCh38/GRCh38.trf.bed"
MIN_SV = 20


def load_bed(path, gap=0):
    d = collections.defaultdict(list)
    for l in open(path):
        if l.startswith(("#", "track")):
            continue
        c, s, e = l.split()[:3]
        d[c].append((int(s), int(e)))
    merged = {}
    for c, L in d.items():
        m = []
        for s, e in sorted(L):
            if m and s <= m[-1][1] + gap:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[c] = m
    return merged


class Intervals:
    def __init__(self, merged):
        self.m = merged
        self.ends = {c: [e for _, e in L] for c, L in merged.items()}

    def overlapping(self, c, s, e):
        L = self.m.get(c, [])
        i = bisect.bisect_right(self.ends.get(c, []), s)  # first with end > s
        out = []
        while i < len(L) and L[i][0] < e:
            out.append(tuple(L[i]))
            i += 1
        return out

    def covered(self, c, s, e):
        return sum(min(e, b) - max(s, a) for a, b in self.overlapping(c, s, e))


def locus_span_fn(trf):
    def locus_span(c, s, e):
        lo, hi = s - 200, e + 200
        for ts, te in trf.overlapping(c, lo - 50, hi + 50):
            lo, hi = min(lo, ts - 50), max(hi, te + 50)
        return max(0, lo), hi

    return locus_span


def same_allele(a, b, r=0.9):
    return a != 0 and (a > 0) == (b > 0) and min(abs(a), abs(b)) / max(abs(a), abs(b)) >= r


def truth_at(vf, c, lo, hi):
    haps = [[], []]
    net = [0, 0]
    net_all = [0, 0]
    missing = False
    n_sym = 0
    for x in vf.fetch(c, lo, hi):
        if not x.alts:
            continue
        alt = x.alts[0]
        if alt.startswith("<"):
            n_sym += 1
            continue
        dlen = len(alt) - len(x.ref)
        gt = x.samples[0]["GT"]
        for h in (0, 1):
            if h >= len(gt) or gt[h] is None:
                if abs(dlen) >= MIN_SV:
                    missing = True
                continue
            if gt[h] > 0:
                net_all[h] += dlen
                if abs(dlen) >= MIN_SV:
                    haps[h].append(f"{dlen:+d}@{x.pos}")
                    net[h] += dlen
    carry = [bool(haps[0]), bool(haps[1])]
    if not any(carry):
        cat = "ref"
    elif carry[0] != carry[1]:
        cat = "het"
    elif same_allele(net[0], net[1]) or (net[0] == net[1] == 0):
        cat = "hom"
    else:
        cat = "cpx_het"
    n_alleles = 1 if cat in ("ref", "hom") else 2
    sub = "-"
    if cat == "cpx_het":
        a, b = net
        if abs(a) < MIN_SV or abs(b) < MIN_SV:
            sub = "one_net_small"
        elif (a > 0) != (b > 0):
            sub = "opposite_sign"
        elif abs(a - b) <= 10:
            sub = "near_hom_le10bp"
        else:
            sub = "size_diff"
    return dict(sub=sub, haps=haps, net=net, net_all=net_all, cat=cat, n_alleles=n_alleles,
                missing=missing, n_sym=n_sym)


