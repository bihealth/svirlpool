#!/usr/bin/env python
"""QC numbers for the trio read labels produced by trio_read_labels.py."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

D = Path(__file__).resolve().parent / "data"
lab = pd.read_csv(D / "trio_read_labels.tsv.gz", sep="\t")
summ = pd.read_csv(D / "trio_container_summary.tsv.gz", sep="\t")
sites = pd.read_csv(D / "trio_informative_sites.tsv.gz", sep="\t")
q = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]


def qs(x):
    return " ".join(f"{k:g}:{v:.3g}" for k, v in zip(q, np.quantile(x, q)))


print(f"sites used: {len(sites)}  (father_gt/mother_gt combos)")
print(sites.groupby(["father_gt", "mother_gt"]).size().to_string())
summ["sites_per_kb"] = summ.n_sites / summ.window_bp * 1000
print("\ninformative sites per container window (+-30kb):", qs(summ.n_sites))
print("sites per kb:", qs(summ.sites_per_kb))
print("containers with 0 sites:", (summ.n_sites == 0).sum(), " <10 sites:", (summ.n_sites < 10).sum())
print("by chr (0-site containers):", summ[summ.n_sites == 0].chr.value_counts().to_dict())

summ["frac_lab"] = (summ.n_pat + summ.n_mat) / summ.n_reads.clip(lower=1)
print("\nreads per container:", qs(summ.n_reads))
print("fraction labelled per container:", qs(summ.frac_lab))
print("overall: rows", len(lab), lab.label.value_counts().to_dict())
print("reads with 0 votes:", (lab.n_pat + lab.n_mat == 0).sum())

L = lab[lab.label != "unknown"]
imp = np.minimum(L.n_pat, L.n_mat) / (L.n_pat + L.n_mat)
print("\npurity min/(pat+mat) over labelled reads:", qs(imp), " mean", round(imp.mean(), 4))
print("  frac ==0:", round((imp == 0).mean(), 4), " frac >0.1:", round((imp > 0.1).mean(), 4))
V = lab[lab.n_pat + lab.n_mat > 0]
allimp = np.minimum(V.n_pat, V.n_mat) / (V.n_pat + V.n_mat)
print("all reads with >=1 vote: minority frac mean", round(allimp.mean(), 4),
      "; votes per read", qs(V.n_pat + V.n_mat))
print("  per-site error estimate (minority votes / all votes, labelled reads):",
      round(np.minimum(L.n_pat, L.n_mat).sum() / (L.n_pat + L.n_mat).sum(), 4))

# read-level consistency across containers
g = L.groupby("read_name").label.nunique()
print("\nreads labelled in >1 container:", (L.read_name.value_counts() > 1).sum(),
      " with conflicting labels:", (g > 1).sum())

s = summ[(summ.n_pat + summ.n_mat) >= 6].copy()
s["pat_frac"] = s.n_pat / (s.n_pat + s.n_mat)
s["p"] = [binomtest(int(a), int(a + b), 0.5).pvalue for a, b in zip(s.n_pat, s.n_mat)]
print(f"\ncontainers with >=6 labelled reads: {len(s)} / {len(summ)}")
print("pat fraction:", qs(s.pat_frac))
auto = s[~s.chr.isin(["chrX", "chrY"])]
print("autosomes pat fraction:", qs(auto.pat_frac), " n", len(auto))
print("binomial p<0.001:", (s.p < 0.001).sum(), " (autosomes", (auto.p < 0.001).sum(), ")")
sk = s[(s.p < 0.001) | (s.pat_frac <= 0.1) | (s.pat_frac >= 0.9)].sort_values("p")
pd.set_option("display.width", 200)
print("extremely skewed containers:")
print(sk[["crID", "chr", "start", "end", "n_sites", "n_reads", "n_pat", "n_mat", "n_unknown", "pat_frac", "p"]]
      .to_string(index=False, max_rows=80))
sx = summ[summ.chr.isin(["chrX", "chrY"])]
print("\nchrX/chrY containers:", len(sx), " pat/mat/unknown:", sx.n_pat.sum(), sx.n_mat.sum(), sx.n_unknown.sum())
if len(sys.argv) > 1:
    for cr in map(int, sys.argv[1:]):
        print(summ[summ.crID == cr].to_string(index=False))
