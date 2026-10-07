"""Can the local read depth recorded in the CR signals (ExtendedSVsignal.coverage,
the depth at the signal) flag the containers that are slow and resolve badly,
before the consensus stage?

Per HG002 container: depth = max over its CRs of the median signal coverage
(signalstrength_to_crs.cr_depth), relative to the median of that over all
containers (the sample's typical depth); joined with container_value.py
(seconds, last-level timeout, Q100 representation, truvari TP / FP calls).

usage: coverage_gate.py <svp root> <variant> [--truth T2TQ100] [--show 2399,3163]
"""

from __future__ import annotations

import argparse
import json
import sqlite3

import container_value as CV
import numpy as np
import pandas as pd


def container_depths(db: str) -> pd.DataFrame:
    rows = {}
    for cid, data in sqlite3.connect(db).execute("select crID, data from containers"):
        crs = json.loads(data)["crs"]
        meds, maxs, nsig, span = [], [], 0, 0
        for cr in crs:
            cov = [s["coverage"] for s in cr["sv_signals"]]
            meds.append(float(np.median(cov)) if cov else 0.0)
            maxs.append(max(cov, default=0))
            nsig += len(cov)
            span += cr["referenceEnd"] - cr["referenceStart"]
        rows[cid] = {"depth": max(meds), "depth_max": max(maxs), "n_signals": nsig,
                     "span": span, "chr": crs[0]["chr"], "start": crs[0]["referenceStart"]}
    return pd.DataFrame.from_dict(rows, orient="index")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("variant")
    ap.add_argument("--truth", default="T2TQ100")
    ap.add_argument("--show", default="")
    a = ap.parse_args()
    base = f"{a.root}/results/{a.variant}/20x"
    work = f"{base}/HG002/work"
    c = CV.containers(work).join(container_depths(work + "/crs_containers.db"), how="inner")
    q = f"{base}/HG002/consensus_q100.tsv"
    try:
        c["repr"] = CV.representation(q).reindex(c.index)
    except FileNotFoundError:
        c["repr"] = np.nan
    tv = f"{base}/truvari/{a.truth}/all"
    try:
        c["tp"] = CV.per_container(CV.call_containers(f"{tv}/tp-comp.vcf.gz")).reindex(c.index).fillna(0)
        c["fp"] = CV.per_container(CV.call_containers(f"{tv}/fp.vcf.gz")).reindex(c.index).fillna(0)
    except Exception:
        c["tp"] = c["fp"] = np.nan
    med = c.depth.median()
    c["rel"] = c.depth / med
    pd.set_option("display.width", 250)
    fmt = lambda x: f"{x:.3f}"  # noqa: E731
    print(f"== {a.variant}: median container depth {med:.1f}")
    if a.show:
        ids = [int(x) for x in a.show.split(",")]
        print(c.loc[c.index.intersection(ids)].to_string(float_format=fmt))
    c["bin"] = pd.cut(c.rel, [0, 1.5, 2, 2.5, 3, 4, 6, np.inf], right=False)
    rows = []
    for b, g in c.groupby("bin", observed=True):
        rows.append({"depth/median": str(b), **CV.summarize(g),
                     "median_s": g.seconds.median(), "gt_120s": int((g.seconds > 120).sum())})
    rows.append({"depth/median": "all", **CV.summarize(c), "median_s": c.seconds.median(),
                 "gt_120s": int((c.seconds > 120).sum())})
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))
    print("\n== a depth gate: containers with depth > k x median excluded")
    rows = []
    for k in (2, 2.5, 3, 4, 6):
        g = c[c.rel > k]
        rows.append({"k": k, **CV.summarize(g),
                     "share_of_hours": g.seconds.sum() / c.seconds.sum(),
                     "of_gt_120s": f"{int((g.seconds > 120).sum())}/{int((c.seconds > 120).sum())}"})
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))
    print("\n== the containers over 120 s: depth / median")
    s = c[c.seconds > 120].sort_values("seconds", ascending=False)
    print(s[["chr", "start", "span", "n_signals", "depth", "rel", "seconds", "last", "repr", "tp", "fp"]]
          .head(40).to_string(float_format=fmt))


if __name__ == "__main__":
    main()
