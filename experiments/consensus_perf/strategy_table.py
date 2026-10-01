"""--clustering-strategy accurate vs fast in svp_improvements (trio): truvari
refined F1 / precision / recall of HG002, Mendelian consistency of the family
call, svirlpool run and consensus-stage times per sample, and the clustering
route of every container.

usage: strategy_table.py [variant ...]   (default: cs_accurate cs_balanced cs_fast)
"""

import glob
import json
import sys
from collections import Counter

import pandas as pd

R = "/home/mayv_c/development/svp_improvements/results"
B = "/home/mayv_c/development/svp_improvements/benchmarks/svirlpool"
SAMPLES = ("HG002", "HG003", "HG004")
variants = sys.argv[1:] or ["cs_accurate", "cs_balanced", "cs_fast"]


def bench_s(path):
    lines = open(path).read().strip().splitlines()
    return float(dict(zip(lines[0].split("\t"), lines[-1].split("\t"), strict=True))["s"])


rows = []
for v in variants:
    row = {"variant": v}
    for ts in ("V5", "T2TQ100"):
        for rs in ("all", "non_trf"):
            s = json.load(open(f"{R}/{v}/20x/truvari/{ts}/{rs}/refine.variant_summary.json"))
            row[f"{ts} {rs} F1"] = s["f1"]
            row[f"{ts} {rs} P"] = s["precision"]
            row[f"{ts} {rs} R"] = s["recall"]
    row["GT conc V5"] = json.load(open(f"{R}/{v}/20x/truvari/V5/all/summary.json")).get(
        "gt_concordance", float("nan"))
    m = pd.read_csv(f"{R}/{v}/20x/mendel/mendel.tsv", sep="\t")
    a = m[(m.svtype == "all") & (m.size_bin == "all") & (m.denominator == "informative")]
    cons = a[a.status == "consistent"]
    inc = a[a.status == "inconsistent"]
    row["MC"] = float(cons.percentage.iloc[0])
    row["MC consistent"] = int(cons["count"].iloc[0])
    row["MC inconsistent"] = int(inc["count"].iloc[0])
    for smp in SAMPLES:
        row[f"run s {smp}"] = bench_s(f"{B}/run.{v}.20x.{smp}.tsv")
        secs = [bench_s(b) for b in glob.glob(
            f"{R}/{v}/20x/{smp}/work/benchmarks/consensus/consensus.batch_*.txt")]
        row[f"consensus s {smp}"] = sum(secs)
        routes = Counter()
        for f in glob.glob(f"{R}/{v}/20x/{smp}/work/consensus/*/consensus.batch_*.jsonl"):
            for line in open(f):
                cds = json.loads(line).get("consensus_dicts") or {}
                ms = {(c.get("clustering_meta_data") or {}).get("method", "none") for c in cds.values()}
                routes[",".join(sorted(ms)) if ms else "no consensus"] += 1
        row[f"routes {smp}"] = dict(routes.most_common())
    rows.append(row)

d = pd.DataFrame(rows).set_index("variant").T
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 120)
print(d.to_string())
