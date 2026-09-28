"""Compare error-rate sweeps of several clustering modes side by side.

Reads the per-container ``eval.tsv`` of each sweep (``eval_rates.py``) and
prints, per added error rate and mode: containers with >= 2 consensuses,
trio pair accuracy, fraction of truth alleles recovered, containers with all
alleles recovered per T2TQ100 category and TRF stratum, consensus
differences and CPU seconds. Then, per category, the recovery of the modes
next to each other.

usage: compare.py <label>=<eval.tsv> [<label>=<eval.tsv> ...] [--out table.tsv]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for rate, g in df.groupby("rate"):
        e = g[g.evaluable]
        row = {
            "rate": rate,
            "multi": int((g.n_consensus >= 2).sum()),
            "cons/cont": round(g.n_consensus.mean(), 2),
            "pair_acc": round(g.pair_acc.mean(), 3),
            "alleles_rec": round(e.recovered.sum() / e.n_truth_alleles.sum(), 3),
        }
        for cat in ("het", "cpx_het", "hom", "ref"):
            row[cat] = round(e[e.cat == cat].all_recovered.mean(), 3)
        row["nonTRF"] = round(e[~e.trf].all_recovered.mean(), 3)
        row["TRF"] = round(e[e.trf].all_recovered.mean(), 3)
        row["cons_diff%"] = round(100 * g["diff"].sum() / max(1, g.aligned.sum()), 2)
        row["cpu_s"] = round(g.secs.sum())
        row["escalated"] = int((g.attempts > 1).sum())
        out.append(row)
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("evals", nargs="+", help="label=path/to/eval.tsv")
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    parts = []
    for spec in a.evals:
        label, path = spec.split("=", 1)
        df = pd.read_csv(path, sep="\t")
        s = summarize(df)
        s.insert(0, "mode", label)
        parts.append(s)
        print(f"== {label}  ({df[df.evaluable].crID.nunique()} evaluable)")
        print(s.drop(columns="mode").to_string(index=False))
        print("   clustering methods (all rates):", df.method.value_counts().to_dict())
        print()
    both = pd.concat(parts)
    for col in ("het", "cpx_het", "hom", "alleles_rec", "pair_acc", "cpu_s"):
        w = both.pivot(index="rate", columns="mode", values=col)
        print(f"-- {col}")
        print(w.to_string())
        print()
    if a.out:
        both.to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
