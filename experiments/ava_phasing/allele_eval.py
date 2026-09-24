"""Allele-level evaluation of read partitions against trio labels + phased truth.

usage: allele_eval.py <batch_eval.tsv> [truth_col_prefix=T2T]
"""
import json, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd
sys.path.insert(0, ".")
import batch_eval as be

res = pd.read_csv(sys.argv[1], sep="\t")
res = res[res["_labels"].notna()]
truth = pd.read_csv("data/container_truth.tsv", sep="\t").set_index("crID")
P = sys.argv[2] if len(sys.argv) > 2 else "T2T"
TRIO = be.TRIO
CUR = be.CUR


def allele_of_haps(cat):
    # pat=hap1, mat=hap2 ; same allele for hom/ref
    return {"pat": "A1", "mat": "A1" if cat in ("hom", "ref") else "A2"}


def evaluate(labels: dict, trio: dict, cat: str, min_cl=3):
    amap = allele_of_haps(cat)
    rows = [(labels[r], amap[trio[r]]) for r in labels if trio.get(r) in amap and labels[r] >= 0]
    n_lab = sum(1 for r in labels if trio.get(r) in amap)
    if len(rows) < 4:
        return None
    by_cl = defaultdict(Counter)
    for cl, al in rows:
        by_cl[cl][al] += 1
    pure = sum(c.most_common(1)[0][1] for c in by_cl.values()) / len(rows)
    sizes = Counter(cl for cl, _ in rows)
    maj = {cl: c.most_common(1)[0][0] for cl, c in by_cl.items() if sizes[cl] >= min_cl}
    alleles = set(amap.values())
    present = {al for _, al in rows}
    recovered = present <= set(maj.values())
    n_big = sum(1 for cl in sizes if sizes[cl] >= min_cl)
    return dict(purity=pure, recovered=recovered, n_big=n_big, assigned=len(rows) / max(n_lab, 1),
                n_alleles=len(present))


out = []
for _, r in res.iterrows():
    c = int(r.crID)
    if c not in truth.index:
        continue
    cat = truth.loc[c, f"{P}_cat"]
    if cat not in ("het", "hom", "cpx_het", "ref") or truth.loc[c, f"{P}_bench_frac"] < 1:
        continue
    trio = TRIO.get(c, {})
    new = {k: v for k, v in json.loads(r["_labels"]).items()}
    cur_raw = CUR.get(str(c), {})
    cur = {rn: k for k, (cid, rns) in enumerate(sorted(cur_raw.items())) for rn in rns}
    # evaluate both on the reads both know about
    common = set(new) & set(cur) if cur else set(new)
    en = evaluate({k: new[k] for k in common}, trio, cat)
    ec = evaluate({k: cur[k] for k in common}, trio, cat) if cur else dict(purity=np.nan, recovered=False, n_big=0, assigned=0, n_alleles=0)
    if en is None or ec is None:
        continue
    out.append(dict(crID=c, cat=cat, trf=bool(truth.loc[c, "trf"]), maxnet=truth.loc[c, f"{P}_max_abs_net"],
                    **{f"new_{k}": v for k, v in en.items()}, **{f"cur_{k}": v for k, v in ec.items()}))
d = pd.DataFrame(out)
d.to_csv(sys.argv[1].replace(".tsv", f".allele_{P}.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200)
g = d.groupby(["cat", "trf"]).agg(n=("crID", "size"),
    new_pur=("new_purity", "mean"), cur_pur=("cur_purity", "mean"),
    new_rec=("new_recovered", "mean"), cur_rec=("cur_recovered", "mean"),
    new_k=("new_n_big", "mean"), cur_k=("cur_n_big", "mean"), new_asg=("new_assigned", "mean"))
print(g.round(3))
tot = d.agg(dict(new_purity="mean", cur_purity="mean", new_recovered="mean", cur_recovered="mean"))
print(tot.round(3).to_dict())
