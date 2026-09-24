"""Per-container truth / diagnostic table for the allele-separation experiment.

For every crs container of an svirlpool run: container span, signals, TRF,
copy number (as consensus uses it), the phased truth per haplotype over the
locus span (T2TQ100 primary, V5 secondary), the consensus clusters svirlpool
made (reads per cluster, net indel of each consensus' alignment over the
locus), the calls, and truvari refine status.

Locus span = container span +-200 bp, extended by overlapping merged TRF
(TRF merged with 500 bp gaps) +-50 bp  (hap_labels.locus_span approach).

usage: container_truth.py <run_dir_snapshot> <out.tsv>
  run_dir_snapshot holds crs_containers.db, consensus_containers.txt,
  consensus.bam, copy_number_tracks.bed.gz, logs/consensus.batch_*.log,
  variants.vcf.gz, truvari/{V5,T2TQ100}/all/refine.{base,comp}.vcf.gz,
  {V5,T2TQ100}.all.bed
"""

import bisect
import collections
import glob
import itertools
import json
import re
import sqlite3
import sys
from pathlib import Path

import pysam

RUN = Path(sys.argv[1])
OUT = Path(sys.argv[2])

TRUTH = {
    "T2T": "/home/mayv_c/development/svp_improvements/resources/truthsets/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz",
    "V5": "/home/mayv_c/biodata/local/benchmarks/HG002_GRCh38_v5.0q_stvar.filtered.vcf.gz",
}
TRUVARI = {"T2T": "T2TQ100", "V5": "V5"}
TRF = "/home/mayv_c/biodata/local/references/GRCh38/GRCh38.trf.bed"
MIN_SV = 20
MIN_CONS_INDEL = 5  # indels on the consensus alignment counted into its net

# ---------------------------------------------------------------- intervals
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


trf = Intervals(load_bed(TRF, gap=500))
bench = {k: Intervals(load_bed(RUN / f"{TRUVARI[k]}.all.bed")) for k in TRUTH}


def locus_span(c, s, e):
    lo, hi = s - 200, e + 200
    for ts, te in trf.overlapping(c, lo - 50, hi + 50):
        lo, hi = min(lo, ts - 50), max(hi, te + 50)
    return max(0, lo), hi


def same_allele(a, b, r=0.9):
    return a != 0 and (a > 0) == (b > 0) and min(abs(a), abs(b)) / max(abs(a), abs(b)) >= r


def allele_match(a, b, r=0.7, ref_tol=MIN_SV):
    """consensus net a vs truth hap net b: both ~ref, or same sign and ratio >= r."""
    if abs(a) < ref_tol and abs(b) < ref_tol:
        return True
    if a == 0 or b == 0 or (a > 0) != (b > 0):
        return False
    return min(abs(a), abs(b)) / max(abs(a), abs(b)) >= r


# ---------------------------------------------------------------- containers
con = sqlite3.connect(f"file:{RUN / 'crs_containers.db'}?mode=ro", uri=True)
containers = {k: json.loads(v) for k, v in con.execute("select crID, data from containers")}

# ---------------------------------------------------------------- consensus logs
log_info = {}
re_prog = re.compile(r"Processing container \(representative crID (\d+)\)")
for f in glob.glob(str(RUN / "logs" / "consensus.batch_*.log")):
    cur = None
    for line in open(f):
        m = re_prog.search(line)
        if m:
            cur = int(m.group(1))
            log_info[cur] = dict(
                cn_log=None, n_reads=None, n_size_outliers=0, lam_kmeans=set(),
                lam_other=0, n_asm_fail=0, n_final_fail=0, rescue=0, isolated=0,
                ava_timeout=0, no_consensus=False, cn_skipped=False,
            )
            continue
        if cur is None:
            continue
        d = log_info[cur]
        if "Maximum copy number for this container:" in line:
            d["cn_log"] = int(line.rsplit(":", 1)[1])
        elif " - number of reads:" in line:
            d["n_reads"] = int(line.rsplit(":", 1)[1])
        elif "KMeans: filtered" in line:
            d["n_size_outliers"] = int(re.search(r"filtered (\d+)", line).group(1))
        elif line.startswith("lamassemble "):
            fa = line.split()[-1]
            m2 = re.search(r"reads\.kmeans\.(\d+)\.", fa)
            if m2:
                d["lam_kmeans"].add(int(m2.group(1)))
            else:
                d["lam_other"] += 1
        elif "consensus assembly failed" in line:
            d["n_asm_fail"] += 1
        elif "final consensus generation failed" in line:
            d["n_final_fail"] += 1
        elif "Trying to rescue" in line and "KMeans" in line:
            d["rescue"] += 1
        elif "Trying to rescue" in line and "isolated" in line:
            d["isolated"] += 1
        elif "AVA alignment failed" in line:
            d["ava_timeout"] += 1
        elif "No consensus objects could be built" in line:
            d["no_consensus"] = True
        elif "exceeds threshold" in line:
            d["cn_skipped"] = True

# ---------------------------------------------------------------- consensus objects
cons = collections.defaultdict(dict)  # rep crID -> consID -> info
for line in open(RUN / "consensus_containers.txt"):
    d = json.loads(line)
    for cid, c in d["consensus_dicts"].items():
        rep = int(cid.split(".")[0])
        reads = {x[2] for x in c["intervals_cutread_alignments"]}
        cons[rep][cid] = dict(n_reads=len(reads), reads=reads, len=len(c["consensus_sequence"]))

# ---------------------------------------------------------------- consensus alignments
bam_segs = collections.defaultdict(list)
with pysam.AlignmentFile(str(RUN / "consensus.bam")) as bam:
    for a in bam:
        if a.is_unmapped or a.is_secondary:
            continue
        bam_segs[a.query_name].append(a)


def cons_net(cid, chrom, lo, hi):
    """net indel (>= MIN_CONS_INDEL) of a consensus' alignments over [lo,hi]."""
    segs = [a for a in bam_segs.get(cid, []) if a.reference_name == chrom
            and a.reference_start < hi and a.reference_end > lo]
    if not segs:
        return None, 0, False
    net = 0
    for a in segs:
        rpos = a.reference_start
        for op, n in a.cigartuples:
            if op in (0, 7, 8):
                rpos += n
            elif op == 2:
                if n >= MIN_CONS_INDEL:
                    ov = min(hi, rpos + n) - max(lo, rpos)
                    if ov > 0:
                        net -= ov
                rpos += n
            elif op == 1:
                if n >= MIN_CONS_INDEL and lo <= rpos <= hi:
                    net += n
    # split alignments: junction between consecutive same-strand segments
    def qint(a):
        ct = a.cigartuples
        lead = ct[0][1] if ct[0][0] in (4, 5) else 0
        trail = ct[-1][1] if ct[-1][0] in (4, 5) else 0
        ql = a.query_alignment_length
        return (trail, trail + ql) if a.is_reverse else (lead, lead + ql)

    segs.sort(key=lambda a: a.reference_start)
    for p, q in zip(segs, segs[1:]):
        if p.is_reverse != q.is_reverse:
            continue
        junc = p.reference_end
        if not (lo <= junc <= hi):
            continue
        rgap = q.reference_start - p.reference_end
        pq, qq = qint(p), qint(q)
        qgap = (pq[0] - qq[1]) if p.is_reverse else (qq[0] - pq[1])
        if abs(qgap - rgap) >= MIN_CONS_INDEL:
            net += qgap - rgap
    covers = min(a.reference_start for a in segs) <= lo and max(a.reference_end for a in segs) >= hi
    return net, len(segs), covers


# ---------------------------------------------------------------- truth
truth_vf = {k: pysam.VariantFile(v) for k, v in TRUTH.items()}


def truth_at(k, c, lo, hi):
    haps = [[], []]
    net = [0, 0]
    net_all = [0, 0]
    missing = False
    n_sym = 0
    for x in truth_vf[k].fetch(c, lo, hi):
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


# ---------------------------------------------------------------- calls + truvari
calls_vf = pysam.VariantFile(str(RUN / "variants.vcf.gz"))
refine_base = {k: pysam.VariantFile(str(RUN / "truvari" / TRUVARI[k] / "all" / "refine.base.vcf.gz")) for k in TRUTH}
refine_comp = {k: pysam.VariantFile(str(RUN / "truvari" / TRUVARI[k] / "all" / "refine.comp.vcf.gz")) for k in TRUTH}
comp_status = {}
for k in TRUTH:
    st = {}
    for x in pysam.VariantFile(str(RUN / "truvari" / TRUVARI[k] / "all" / "refine.comp.vcf.gz")):
        st[x.id] = x.samples[0]["BD"]
    comp_status[k] = st
calls_by_cons = collections.defaultdict(list)
for x in calls_vf:
    for cid in x.info.get("CONSENSUSIDs", ()):
        calls_by_cons[int(cid.split(":")[-1].split(".")[0])].append(x.id)

ALN = pysam.AlignmentFile(json.load(open(RUN / "config.json"))["alignments"])


def mapq_profile(c, s, e):
    """primary alignments over the container span: count and fraction MAPQ < 10."""
    n = low = 0
    for a in ALN.fetch(c, s, e):
        if a.is_secondary or a.is_supplementary or a.is_unmapped:
            continue
        n += 1
        low += a.mapping_quality < 10
    return n, round(low / n, 3) if n else 0.0


cn_tbx = pysam.TabixFile(str(RUN / "copy_number_tracks.bed.gz"))


def cn_raw(regions):
    m = 0
    for c, s, e in regions:
        try:
            for row in cn_tbx.fetch(c, s, e):
                m = max(m, int(row.split("\t")[3]))
        except ValueError:
            pass
    return m


def fmt_gt(gt):
    return "/".join("." if a is None else str(a) for a in gt)


# ---------------------------------------------------------------- main loop
rows = []
for rep in sorted(containers):
    crs = containers[rep]["crs"]
    c = crs[0]["chr"]
    s = min(cr["referenceStart"] for cr in crs)
    e = max(cr["referenceEnd"] for cr in crs)
    lo, hi = locus_span(c, s, e)
    sigs = [x for cr in crs for x in cr["sv_signals"]]
    regions = [(cr["chr"], cr["referenceStart"], cr["referenceEnd"]) for cr in crs]
    trf_iv = trf.overlapping(c, s, e)
    li = log_info.get(rep, {})
    raw = cn_raw(regions)
    n_aln, frac_lowmapq = mapq_profile(c, s, e)
    r = dict(
        crID=rep, chr=c, start=s, end=e, span_len=e - s, locus_start=lo, locus_end=hi,
        n_crs=len(crs), n_sv_signals=len(sigs),
        n_signal_reads=len({x["readname"] for x in sigs}),
        n_reads=li.get("n_reads"),
        trf=bool(trf_iv), trf_merged_bp=sum(b - a for a, b in trf_iv),
        trf_frac_span=round(trf.covered(c, s, e) / max(1, e - s), 3),
        cn_raw=raw, cn=max(2, raw), cn_log=li.get("cn_log"),
        n_aln_primary=n_aln, frac_mapq_lt10=frac_lowmapq,
        # no segdup annotation available: CN>2, >=20% MAPQ<10, or >=2x median depth (20x)
        segdup_proxy=max(2, raw) > 2 or frac_lowmapq >= 0.2 or n_aln >= 40,
    )
    for k in ("T2T", "V5"):
        t = truth_at(k, c, lo, hi)
        r[f"{k}_cat"] = t["cat"]
        r[f"{k}_n_alleles"] = t["n_alleles"]
        r[f"{k}_cpx_sub"] = t["sub"]
        r[f"{k}_max_abs_net"] = max(abs(t["net"][0]), abs(t["net"][1]))
        r[f"{k}_hap1"] = ",".join(t["haps"][0]) or "-"
        r[f"{k}_hap2"] = ",".join(t["haps"][1]) or "-"
        r[f"{k}_net"] = f"{t['net'][0]},{t['net'][1]}"
        r[f"{k}_net_all"] = f"{t['net_all'][0]},{t['net_all'][1]}"
        r[f"{k}_hap_missing"] = t["missing"]
        r[f"{k}_bench_frac"] = round(bench[k].covered(c, lo, hi) / (hi - lo), 3)
        r[f"_{k}_t"] = t
    # svirlpool clusters
    cs = cons.get(rep, {})
    cids = sorted(cs, key=lambda x: int(x.split(".")[1]))
    r["n_consensus"] = len(cs)
    r["cons_ids"] = ",".join(cids) or "-"
    r["reads_per_cons"] = ",".join(str(cs[x]["n_reads"]) for x in cids) or "-"
    r["n_reads_in_cons"] = len(set().union(*[cs[x]["reads"] for x in cids])) if cids else 0
    nets = []
    ncov = 0
    split = 0
    for x in cids:
        n, nseg, cov = cons_net(x, c, lo, hi)
        nets.append(n)
        ncov += cov
        split += nseg > 1
    r["cons_net"] = ",".join("NA" if n is None else str(n) for n in nets) or "-"
    r["n_cons_unaligned_locus"] = sum(n is None for n in nets)
    r["n_cons_covering_locus"] = ncov
    r["n_cons_split_aln"] = split
    # distinct consensus alleles (greedy: same allele if allele_match r=0.9)
    reps_ = []
    for n in nets:
        if n is None:
            continue
        if not any(allele_match(n, m, r=0.9) for m in reps_):
            reps_.append(n)
    r["n_distinct_cons_alleles"] = len(reps_)
    km = li.get("lam_kmeans", set())
    method = ("kmeans" if km else "") + ("+" if km and li.get("lam_other") else "") + ("other" if li.get("lam_other") else "")
    r["cluster_method"] = method or "none"
    r["kmeans_k"] = len(km)
    r["n_size_outliers"] = li.get("n_size_outliers", 0)
    r["n_asm_failed"] = li.get("n_asm_fail", 0) + li.get("n_final_fail", 0)
    r["asm_rescue"] = li.get("rescue", 0) + li.get("isolated", 0)
    r["ava_timeout"] = li.get("ava_timeout", 0)
    r["no_consensus"] = len(cs) == 0
    r["cn_skipped"] = li.get("cn_skipped", False)
    # truth alleles recovered by some consensus
    for k in ("T2T", "V5"):
        t = r[f"_{k}_t"]
        tn = [t["net"][0]] if t["cat"] in ("ref", "hom") else t["net"]
        # one-to-one: each consensus can stand for one truth allele only
        cn_ok = [n for n in nets if n is not None]
        best = 0
        for perm in itertools.permutations(range(len(cn_ok)), min(len(cn_ok), len(tn))):
            best = max(best, sum(allele_match(cn_ok[j], tn[i]) for i, j in enumerate(perm)))
        r[f"{k}_alleles_recovered"] = best
    # calls
    calls = []
    n_pass = 0
    gt_wrong = 0
    for x in calls_vf.fetch(c, lo, hi):
        svlen = x.info.get("SVLEN")
        svlen = svlen[0] if isinstance(svlen, tuple) else svlen
        st = x.info.get("SVTYPE")
        gt = fmt_gt(x.samples[0]["GT"])
        filt = ",".join(x.filter.keys()) or "."
        signed = abs(svlen) * (1 if st == "INS" else -1) if st in ("INS", "DEL") else 0
        # "ref" = inside a truvari refine region (status only on the harmonised records)
        v5 = comp_status["V5"].get(x.id, "ref" if filt == "PASS" else "-")
        t2 = comp_status["T2T"].get(x.id, "ref" if filt == "PASS" else "-")
        exp = ""
        if filt == "PASS":
            n_pass += 1
            if signed:
                t = r["_T2T_t"]
                nm = sum(allele_match(signed, h, r=0.7, ref_tol=0) for h in t["net"])
                exp = {0: "0/0", 1: "0/1", 2: "1/1"}[nm]
                if exp != "0/0" and exp != gt.replace("|", "/"):
                    gt_wrong += 1
        calls.append(f"{st}{signed:+d}:{gt}:{filt}:DV{x.samples[0].get('DV')}/TC{x.samples[0].get('TC')}:V5{v5}:T2T{t2}" + (f":hapGT{exp}" if exp else ""))
    r["n_calls_pass"] = n_pass
    r["calls"] = ";".join(calls) or "-"
    r["n_calls_from_container"] = len(calls_by_cons.get(rep, []))
    r["n_calls_gt_wrong_T2T"] = gt_wrong
    for k in ("V5", "T2T"):
        tp = fn = 0
        lst = []
        for x in refine_base[k].fetch(c, lo, hi):
            bd = x.samples[0]["BD"]
            svlen = len(x.alts[0]) - len(x.ref) if x.alts and not x.alts[0].startswith("<") else 0
            lst.append(f"{svlen:+d}:{fmt_gt(x.samples[0]['GT'])}:{bd}")
            tp += bd == "TP"
            fn += bd == "FN"
        r[f"truv_{k}_TP"] = tp
        r[f"truv_{k}_FN"] = fn
        # comp side from refine.comp by position: inside refined regions the
        # calls are phab-harmonised records without the original ID
        r[f"truv_{k}_FP"] = sum(x.samples[0]["BD"] == "FP" for x in refine_comp[k].fetch(c, lo, hi))
        r[f"truv_{k}_base"] = ";".join(lst) or "-"
    for k in ("T2T", "V5"):
        nt, nc = r[f"{k}_n_alleles"], r["n_consensus"]
        r[f"cluster_vs_truth_{k}"] = "under" if nc < nt else ("equal" if nc == nt else "over")
    r["cluster_vs_truth"] = r["cluster_vs_truth_T2T"]
    for k in ("T2T", "V5"):
        del r[f"_{k}_t"]
    rows.append(r)

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    cols = list(rows[0])
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"wrote {len(rows)} rows to {OUT}")
