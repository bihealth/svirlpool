"""Candidate regions per container (one consensus job): a region subset
(svp_tiered15, 15% of the autosomes) vs a whole-genome run of the same sample.

Containers join CRs only through split reads (one significant BND in each,
>= max(3, depth / 3) reads per pair). A region subset produces no CR where a
read's other alignment falls outside the regions, so its links are cut.

usage: container_sizes.py <crs_containers.db> [<crs_containers.db> ...]
         [--regions regions15.bed]   (also: whole-genome containers touching the regions)
"""

from __future__ import annotations

import argparse
import ast
import json
import sqlite3

import numpy as np
import pandas as pd


def containers(path: str) -> pd.DataFrame:
    rows = []
    for cid, data in sqlite3.connect(path).execute("select crID, data from containers"):
        d = json.loads(data)
        crs = d["crs"]
        if isinstance(crs, str):  # serialized as a python repr
            crs = ast.literal_eval(crs)
        chrs = [c["chr"] for c in crs]
        rows.append({
            "cid": cid, "n_crs": len(crs), "n_chr": len(set(chrs)),
            "crs": [(c["chr"], c["referenceStart"], c["referenceEnd"]) for c in crs],
            "cr_bp": sum(c["referenceEnd"] - c["referenceStart"] for c in crs),
            "n_reads": len({s["readname"] for c in crs for s in c["sv_signals"]}),
        })
    return pd.DataFrame(rows)


def summary(name: str, c: pd.DataFrame) -> dict:
    n = c.n_crs.values
    return {
        "set": name, "containers": len(c), "crs": int(n.sum()),
        "crs_per_container": n.mean(),
        "containers_gt1_pct": 100 * (n > 1).mean(),
        "crs_in_multi_pct": 100 * n[n > 1].sum() / n.sum(),
        "containers_ge10": int((n >= 10).sum()),
        "crs_in_ge10_pct": 100 * n[n >= 10].sum() / n.sum(),
        "multi_chrom": int((c.n_chr > 1).sum()),
        "max_crs": int(n.max()),
        "max_reads": int(c.n_reads.max()),
        "p99_crs": float(np.quantile(n, .99)),
    }


def touching(c: pd.DataFrame, bed: str) -> pd.DataFrame:
    b = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
    by = {k: g[["s", "e"]].values for k, g in b.groupby("c")}

    def inside(cr):
        a = by.get(cr[0])
        return a is not None and bool(((a[:, 0] < cr[2]) & (a[:, 1] > cr[1])).any())

    return c[[any(inside(cr) for cr in crs) for crs in c.crs]]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dbs", nargs="+")
    ap.add_argument("--regions")
    a = ap.parse_args()
    rows, big = [], []
    for db in a.dbs:
        c = containers(db)
        name = db.split("/")[-5:-1]
        rows.append(summary("/".join(name), c))
        if a.regions:
            rows.append(summary("/".join(name) + " touching regions", touching(c, a.regions)))
        top = c.sort_values("n_crs", ascending=False).head(8)
        for r in top.itertuples():
            chrs = sorted({x[0] for x in r.crs})
            big.append({"set": "/".join(name), "cid": r.cid, "n_crs": r.n_crs, "n_reads": r.n_reads,
                        "chroms": ",".join(chrs[:6]) + ("..." if len(chrs) > 6 else ""),
                        "cr_kb": r.cr_bp / 1000})
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print()
    print(pd.DataFrame(big).to_string(index=False, float_format=lambda x: f"{x:.1f}"))
