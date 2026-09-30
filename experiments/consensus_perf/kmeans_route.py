"""Which loci can take the KMeans route instead of the read phasing?

Reads kmeans_gate.py output (with per-read summed indels) and describes each
accepted KMeans partition by features known before any phasing:

  k          accepted clusters
  spread     max over clusters of the mean read distance to the centroid (bp)
  sep        min distance between centroids (bp; k >= 2)
  min_frac   smallest cluster's share of the reads
  net        k = 1: distance of the centroid from (0, 0), i.e. the SV size
  trf        share of the container's SV signals in tandem repeats

Then it compares the containers where KMeans loses an SV allele that the
phasing recovers with the rest, and scores candidate routing rules: share of
containers routed to KMeans, phasing time saved, alleles lost / gained vs
phasing on the containers with truth.

usage: kmeans_route.py <kmeans_gate.tsv>
"""

import json
import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 220)
d = pd.read_csv(sys.argv[1], sep="\t")
d = d[~d.skipped.astype(bool)].copy()


def features(r):
    if r.gate_k == 0:
        return pd.Series(dtype=float)
    groups = json.loads(r.kmeans_groups)
    indels = json.loads(r.indels)
    names = [rn for rn in groups if rn in indels]
    X = np.array([indels[rn] for rn in names], dtype=float)
    lab = np.array([groups[rn] for rn in names])
    cents, spreads, fracs = [], [], []
    for c in sorted(set(lab)):
        pts = X[lab == c]
        cen = pts.mean(axis=0)
        cents.append(cen)
        spreads.append(np.linalg.norm(pts - cen, axis=1).mean())
        fracs.append(len(pts) / len(X))
    sep = min(
        (np.linalg.norm(a - b) for i, a in enumerate(cents) for b in cents[i + 1 :]),
        default=np.nan,
    )
    return pd.Series({
        "spread": max(spreads), "sep": sep, "min_frac": min(fracs),
        "net": np.linalg.norm(cents[0]) if len(cents) == 1 else np.nan,
        "n_km_reads": len(X),
    })


d = d.join(d.apply(features, axis=1))
d["truth"] = d.bench.eq(True) & d.km_recovered.notna() & d.ph_recovered.notna()
t = d[d.truth].copy()
t["kr"] = t.km_recovered.astype(bool)
t["pr"] = t.ph_recovered.astype(bool)
t["lost"] = ~t.kr & t.pr
t["gained"] = t.kr & ~t.pr

def q(s):
    return " / ".join(f"{v:.2f}" for v in s.quantile([0.1, 0.5, 0.9]))


print("feature quantiles (q10 / q50 / q90), lost vs rest, per k")
for k in (1, 2):
    g = t[t.gate_k == k]
    cols = ["spread", "net"] if k == 1 else ["spread", "sep", "min_frac"]
    for c in cols + ["sig_trf_frac"]:
        print(f"  k={k} {c:12s} lost ({g.lost.sum():3d}): {q(g[g.lost][c]):24s} rest ({(~g.lost).sum()}): {q(g[~g.lost][c])}")
print("\nk = 1 containers by truth category:")
print(pd.crosstab(t[t.gate_k == 1].cat, t[t.gate_k == 1].lost).to_string())

total_phase = d.t_phase.sum()


def score(name, route):
    """route: boolean Series over d (True = KMeans instead of phasing)."""
    r = route.reindex(t.index).fillna(False).astype(bool)
    print(f"  {name:58s} routed {route.mean():6.1%}  saved {d.t_phase[route].sum() / total_phase:6.1%}"
          f"  lost {int((t.lost & r).sum()):4d}  gained {int((t.gained & r).sum()):3d}"
          f"  (truth routed {r.sum()})")


acc = d.gate_k > 0
k1, k2 = d.gate_k == 1, d.gate_k >= 2
nontr = d.sig_trf_frac < 0.5
print("\nrouting rules (phasing time is the 1-thread sum over all containers):")
score("gate accepts (k >= 1)", acc)
score("k >= 2", k2)
for s in (50, 100, 200):
    score(f"k >= 2, sep >= {s}", k2 & (d.sep >= s))
for s in (50, 100):
    for f in (0.2, 0.3):
        score(f"k >= 2, sep >= {s}, min_frac >= {f}", k2 & (d.sep >= s) & (d.min_frac >= f))
score("k >= 2, sep >= 4 * spread", k2 & (d.sep >= 4 * d.spread))
score("k >= 2, sep >= 8 * spread", k2 & (d.sep >= 8 * d.spread))
score("k >= 2, non-TR (< 50% signals in TRF)", k2 & nontr)
score("k = 1, net >= 50, spread <= 5", k1 & (d.net >= 50) & (d.spread <= 5))
score("k = 1, spread <= 2", k1 & (d.spread <= 2))
score("k = 1, non-TR", k1 & nontr)
score("k = 1, non-TR, net >= 50", k1 & nontr & (d.net >= 50))
score("non-TR, any accepted k", acc & nontr)

print("\nrules on the smallest-cluster share and the k = 1 SV size:")
for f in (0.15, 0.2, 0.25, 0.3):
    score(f"k >= 2, min_frac >= {f}", k2 & (d.min_frac >= f))
for f in (0.2, 0.25):
    score(f"k >= 2, min_frac >= {f}, sep >= 30", k2 & (d.min_frac >= f) & (d.sep >= 30))
for n in (50, 100, 200):
    score(f"k = 1, net >= {n}", k1 & (d.net >= n))
    score(f"k = 1, net >= {n}, spread <= 10", k1 & (d.net >= n) & (d.spread <= 10))
score("k = 1, non-TR, net >= 50", k1 & nontr & (d.net >= 50))
for f in (0.2, 0.25):
    base = k2 & (d.min_frac >= f)
    score(f"[k>=2, min_frac>={f}] or [k=1, non-TR, net>=50]", base | (k1 & nontr & (d.net >= 50)))
    score(f"[k>=2, min_frac>={f}] or [k=1, net>=100, spread<=10]", base | (k1 & (d.net >= 100) & (d.spread <= 10)))
