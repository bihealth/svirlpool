"""Summaries of container_truth.tsv and a stratified test set of containers.

usage: select_test_set.py data/container_truth.tsv data/test_set.tsv
"""

import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 200)
d = pd.read_csv(sys.argv[1], sep="\t")
rng = np.random.default_rng(7)

# ------------------------------------------------------------------ summary
print("== T2T category x n_consensus")
print(pd.crosstab(d.T2T_cat, d.n_consensus, margins=True))
print("\n== V5 category x n_consensus")
print(pd.crosstab(d.V5_cat, d.n_consensus, margins=True))
d["trf_s"] = np.where(d.trf, "TRF", "nonTRF")
d["cn_s"] = np.where(d.cn > 2, "CN>2", np.where(d.segdup_proxy, "CN2_segdup_proxy", "CN2"))
print("\n== cluster_vs_truth (T2T) by category x TRF")
print(pd.crosstab([d.T2T_cat, d.trf_s], d.cluster_vs_truth, margins=True))
print("\n== cluster_vs_truth (T2T) by category x CN / segdup proxy")
print(pd.crosstab([d.T2T_cat, d.cn_s], d.cluster_vs_truth, margins=True))
print("\n== mismatch rate by TRF x CN")
print(d.assign(mis=d.cluster_vs_truth != "equal").groupby(["trf_s", "cn_s"]).mis.agg(["sum", "count", "mean"]))
print("\n== cluster method x cluster_vs_truth")
print(pd.crosstab([d.cluster_method, d.T2T_cat], d.cluster_vs_truth))
print("\n== truth alleles recovered by a consensus (T2T) by category x n_consensus")
print(pd.crosstab([d.T2T_cat, d.n_consensus], d.T2T_alleles_recovered))
g = d.groupby(["T2T_cat", "cluster_vs_truth"])
print("\n== truvari refine (V5) + GT errors by category x cluster_vs_truth")
print(g[["truv_V5_TP", "truv_V5_FN", "truv_V5_FP", "n_calls_gt_wrong_T2T"]].sum().assign(n=g.size()))
print("\n== T2T vs V5 category")
print(pd.crosstab(d.T2T_cat, d.V5_cat))

for k in ("T2T", "V5"):
    d[f"{k}_max_abs_net_all"] = d[f"{k}_net_all"].map(lambda v: max(abs(int(x)) for x in v.split(",")))
print("\n== truth 'ref' containers whose small (<20 bp) indels sum to >= 20 bp on a haplotype (T2T):",
      int(((d.T2T_cat == "ref") & (d.T2T_max_abs_net_all >= 20)).sum()), "of", int((d.T2T_cat == "ref").sum()))

# ------------------------------------------------------------------ test set
clean = (
    (d.n_reads <= 100)
    & (d.T2T_bench_frac >= 0.99)
    & (d.V5_bench_frac >= 0.99)
    & ~d.T2T_hap_missing
    & (d.n_consensus > 0)
)
agree = d.T2T_cat == d.V5_cat


def wrong_reasons(r):
    out = []
    if r.cluster_vs_truth != "equal":
        out.append(f"clusters_{r.cluster_vs_truth}({r.n_consensus} vs {r.T2T_n_alleles} truth alleles)")
    if r.T2T_alleles_recovered < r.T2T_n_alleles:
        out.append(f"allele_missing({r.T2T_alleles_recovered}/{r.T2T_n_alleles} recovered)")
    if r.truv_V5_FN > 0:
        out.append(f"V5_FN={r.truv_V5_FN}")
    if r.truv_V5_FP > 0:
        out.append(f"V5_FP={r.truv_V5_FP}")
    if r.n_calls_gt_wrong_T2T > 0:
        out.append(f"GT_wrong={r.n_calls_gt_wrong_T2T}")
    if r.T2T_cat == "ref" and r.n_calls_pass > 0:
        out.append(f"PASS_calls_at_ref={r.n_calls_pass}")
    if r.n_asm_failed > 0:
        out.append(f"asm_failed={r.n_asm_failed}")
    return out


d["wrong_reasons"] = d.apply(lambda r: ";".join(wrong_reasons(r)), axis=1)
d["is_wrong"] = d.wrong_reasons != ""

strata = [
    ("het_nonTRF", 10, clean & agree & (d.T2T_cat == "het") & ~d.trf & ~d.segdup_proxy & (d.T2T_max_abs_net >= 30)),
    ("hom_nonTRF", 10, clean & agree & (d.T2T_cat == "hom") & ~d.trf & ~d.segdup_proxy & (d.T2T_max_abs_net >= 30)),
    ("het_TRF", 10, clean & agree & (d.T2T_cat == "het") & d.trf & ~d.segdup_proxy & (d.T2T_max_abs_net >= 30)),
    ("hom_TRF", 8, clean & agree & (d.T2T_cat == "hom") & d.trf & ~d.segdup_proxy & (d.T2T_max_abs_net >= 30)),
    ("cpx_het", 8, clean & agree & (d.T2T_cat == "cpx_het") & ~d.segdup_proxy
     & d.T2T_cpx_sub.isin(["size_diff", "opposite_sign"]) & (d.T2T_max_abs_net >= 50)),
    ("ref", 6, clean & agree & (d.T2T_cat == "ref") & ~d.segdup_proxy
     & (d.T2T_max_abs_net_all < 10) & (d.V5_max_abs_net_all < 10)),
    # all CN>2 containers lie outside the benchmark regions -> 5 benchmarked
    # segdup-proxy containers here, plus one CN=4 container appended below
    ("cn_gt2_or_segdup", 5, (d.n_reads <= 100) & (d.n_consensus > 0) & d.segdup_proxy
     & (d.cn == 2) & (d.T2T_bench_frac >= 0.99) & (d.V5_bench_frac >= 0.99)),
]


def primary_kind(reasons):
    return reasons.split(";")[0].split("(")[0].split("=")[0] if reasons else "right"


picked = []
used = set()
for name, n, mask in strata:
    pool = d[mask & ~d.crID.isin(used)]
    wrong = pool[pool.is_wrong]
    right = pool[~pool.is_wrong]
    n_wrong = min(len(wrong), (n + 1) // 2)
    # diversify the wrong ones over their primary failure kind, round robin
    kinds = {}
    for idx in rng.permutation(len(wrong)):
        r = wrong.iloc[idx]
        kinds.setdefault(primary_kind(r.wrong_reasons), []).append(r)
    sel_w = []
    while len(sel_w) < n_wrong and any(kinds.values()):
        for k in sorted(kinds, key=lambda k: -len(kinds[k])):
            if kinds[k] and len(sel_w) < n_wrong:
                sel_w.append(kinds[k].pop())
    n_right = min(len(right), n - len(sel_w))
    sel_r = [right.iloc[i] for i in rng.choice(len(right), n_right, replace=False)] if n_right else []
    for r, kind in [(x, "wrong") for x in sel_w] + [(x, "right") for x in sel_r]:
        used.add(r.crID)
        if kind == "wrong":
            reason = "svirlpool wrong: " + r.wrong_reasons
        else:
            reason = "svirlpool right: cluster count = truth alleles, all truth alleles recovered, no FN/FP/GT error"
        picked.append(dict(
            crID=r.crID, stratum=name, outcome=kind, reason=reason, n_reads=r.n_reads,
            chr=r.chr, start=r.start, end=r.end, trf=r.trf, cn=r.cn, segdup_proxy=r.segdup_proxy,
            T2T_cat=r.T2T_cat, T2T_net=r.T2T_net, V5_net=r.V5_net,
            n_consensus=r.n_consensus, reads_per_cons=r.reads_per_cons, cons_net=r.cons_net,
            cluster_method=r.cluster_method, calls=r.calls,
        ))
    print(f"{name}: pool {len(pool)} (wrong {len(wrong)}, right {len(right)}) -> picked {len(sel_w)} wrong + {len(sel_r)} right")

# one CN=4 container (chr18 p-arm start, outside both benchmarks: truth unreliable)
r = d[(d.cn > 2) & (d.n_reads <= 100) & (d.n_consensus >= 2) & (d.T2T_cat == "het")].sort_values("crID").iloc[0]
picked.append(dict(
    crID=r.crID, stratum="cn_gt2_or_segdup", outcome="wrong" if r.is_wrong else "right",
    reason=f"CN={r.cn} (tracks), locus outside the T2TQ100/V5 benchmark regions, truth labels unreliable; "
    + (r.wrong_reasons or "cluster count matches truth"),
    n_reads=r.n_reads, chr=r.chr, start=r.start, end=r.end, trf=r.trf, cn=r.cn, segdup_proxy=r.segdup_proxy,
    T2T_cat=r.T2T_cat, T2T_net=r.T2T_net, V5_net=r.V5_net, n_consensus=r.n_consensus,
    reads_per_cons=r.reads_per_cons, cons_net=r.cons_net, cluster_method=r.cluster_method, calls=r.calls,
))
ts = pd.DataFrame(picked)
ts.to_csv(sys.argv[2], sep="\t", index=False)
print(f"\nwrote {len(ts)} containers to {sys.argv[2]}")
print(ts.groupby(["stratum", "outcome"]).size().unstack(fill_value=0))
print("n_reads: ", ts.n_reads.describe().to_dict())
