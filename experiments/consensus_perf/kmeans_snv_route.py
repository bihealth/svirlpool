"""Routing rules with the cheap reference-SNV split (``ref_snv_haplotypes``).

Reads kmeans_gate.py output (with snv_* columns) and scores, per rule, the
partition each container ends up with -- KMeans, the SNV split, or the read
phasing -- against always phasing: containers routed away from the phasing,
phasing time saved net of the SNV check, SV alleles lost / gained vs phasing,
and the mean allele recovery and purity over the containers with truth.

  hap_diff   distance between the two SNV haplotypes' median summed indels
             (bp); a het SV the KMeans missed shows up here

usage: kmeans_snv_route.py <kmeans_gate.tsv>
"""

import json
import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 220)
d = pd.read_csv(sys.argv[1], sep="\t")
d = d[~d.skipped.astype(bool)].copy()


def km_features(r):
    out = {"min_frac": np.nan, "net": np.nan, "hap_diff": np.nan}
    indels = json.loads(r.indels)
    if r.gate_k > 0:
        groups = json.loads(r.kmeans_groups)
        lab = np.array([groups[rn] for rn in groups if rn in indels])
        out["min_frac"] = np.bincount(lab).min() / len(lab) if len(set(lab)) > 1 else 1.0
        if r.gate_k == 1:
            X = np.array([indels[rn] for rn in groups if rn in indels], dtype=float)
            out["net"] = float(np.linalg.norm(X.mean(axis=0)))
    if r.snv_split:
        haps = json.loads(r.snv_groups)
        med = [
            np.median(np.array([indels[rn] for rn, h in haps.items() if h == k and rn in indels],
                               dtype=float).reshape(-1, 2), axis=0)
            for k in (0, 1)
        ]
        out["hap_diff"] = float(np.linalg.norm(med[0] - med[1]))
    return pd.Series(out)


d = d.join(d.apply(km_features, axis=1))
d["truth"] = d.bench.eq(True) & d.ph_recovered.notna()
for c in ("km", "ph", "snv"):
    d[f"{c}_rec"] = d[f"{c}_recovered"].map({True: True, False: False, "True": True, "False": False})
total_phase = d.t_phase.sum()
k1, k2 = d.gate_k == 1, d.gate_k >= 2
nontr = d.sig_trf_frac < 0.5
split = d.snv_split.astype(bool)

print(f"containers {len(d)}; SNV split found {split.mean():.1%}; SNV check {d.t_snv.sum():.0f} s "
      f"(mean {d.t_snv.mean() * 1000:.0f} ms) vs phasing {total_phase:.0f} s (mean {d.t_phase.mean() * 1000:.0f} ms)")
t = d[d.truth & d.snv_rec.notna()]
print(f"\nSNV split alone vs phasing, containers with truth and an SNV split ({len(t)}):")
print(f"  recovered snv {t.snv_rec.mean():.3f} / phasing {t.ph_rec.mean():.3f}; purity snv {t.snv_purity.mean():.3f} / "
      f"phasing {t.ph_purity.mean():.3f}; assigned snv {t.snv_assigned.mean():.3f} / phasing {t.ph_assigned.mean():.3f}")

g = d[k1 & d.truth & d.km_rec.notna()].copy()
g["km_lost"] = ~g.km_rec.astype(bool) & g.ph_rec.astype(bool)
print("\nk = 1 containers: SNV split and haplotype indel difference, KMeans lost vs rest")
print(pd.crosstab([g.km_lost, g.snv_split], g.cat).to_string())
q = lambda s: " / ".join(f"{v:.1f}" for v in s.quantile([0.1, 0.5, 0.9]))  # noqa: E731
print(f"  hap_diff q10/q50/q90  lost & split: {q(g[g.km_lost & g.snv_split].hap_diff)}   "
      f"rest & split: {q(g[~g.km_lost & g.snv_split].hap_diff)}")


def score(name, route):
    """route: Series of 'km' / 'snv' / 'ph' per container."""
    rec = pd.Series(np.nan, index=d.index, dtype=object)
    pur = pd.Series(np.nan, index=d.index)
    for c in ("km", "snv", "ph"):
        m = route == c
        rec[m] = d.loc[m, f"{c}_rec"]
        pur[m] = d.loc[m, f"{c}_purity"]
    tt = d.truth & rec.notna() & d.ph_rec.notna()
    r, p = rec[tt].astype(bool), d.ph_rec[tt].astype(bool)
    away = route != "ph"
    snv_cost = d.t_snv[d._snv_needed].sum() if "_snv_needed" in d else 0.0
    saved = (d.t_phase[away].sum() - snv_cost) / total_phase
    print(f"  {name:62s} away {away.mean():6.1%}  saved {saved:6.1%}  lost {int((~r & p).sum()):4d}  "
          f"gained {int((r & ~p).sum()):3d}  recovered {r.mean():.4f}  purity {pur[tt].mean():.4f}")


def rule(k2_km, k1_route):
    route = pd.Series("ph", index=d.index)
    route[k2 & k2_km] = "km"
    route[k1] = k1_route[k1]
    return route


print("\nrouting rules (baseline: always phase; saved = phasing time avoided - SNV checks run):")
d["_snv_needed"] = False
score("always phase", pd.Series("ph", index=d.index))
k2ok = d.min_frac >= 0.2
score("previous best: k>=2 minfrac>=0.2 | k=1 non-TR net>=50",
      rule(k2ok, pd.Series(np.where(nontr & (d.net >= 50), "km", "ph"), index=d.index)))
d["_snv_needed"] = k1
for D in (5, 10, 20, 40):
    score(f"k>=2 minfrac>=0.2 | k=1: km unless SNV split & hap_diff >= {D} (then phase)",
          rule(k2ok, pd.Series(np.where(split & (d.hap_diff >= D), "ph", "km"), index=d.index)))
for D in (5, 10, 20):
    score(f"   same, but SNV partition instead of phasing (hap_diff >= {D})",
          rule(k2ok, pd.Series(np.where(split & (d.hap_diff >= D), "snv", "km"), index=d.index)))
score("k>=2 minfrac>=0.2 | k=1: km unless SNV split (then phase)",
      rule(k2ok, pd.Series(np.where(split, "ph", "km"), index=d.index)))
d["_snv_needed"] = True
score("SNV split -> SNV partition, else phase", pd.Series(np.where(split, "snv", "ph"), index=d.index))
for D in (5, 10):
    route = rule(k2ok, pd.Series(np.where(split & (d.hap_diff >= D), "snv", "km"), index=d.index))
    route[(d.gate_k == 0) | (k2 & ~k2ok)] = np.where(split[(d.gate_k == 0) | (k2 & ~k2ok)], "snv", "ph")
    score(f"k-means rule, k=1 hap_diff>={D} -> SNV, rejected/lopsided: SNV if split else phase", route)

print("\ncombinations: KMeans where reliable, the SNV partition where a split exists, phasing for the rest")
d["_snv_needed"] = True
for S in (2, 5, 10, 20):
    ok = split & (d.snv_sites >= S)
    score(f"SNV split with >= {S} sites -> SNV, else phase", pd.Series(np.where(ok, "snv", "ph"), index=d.index))
    route = pd.Series(np.where(ok, "snv", "ph"), index=d.index)
    route[k2 & k2ok] = "km"
    score(f"k>=2 minfrac>=0.2 -> KMeans; SNV (>= {S} sites) -> SNV; else phase", route)
    route[k1 & ~ok] = "km"
    score("   + k=1 without SNV split -> KMeans", route)
d["_snv_needed"] = ~(k2 & k2ok)
route = pd.Series(np.where(split, "snv", "ph"), index=d.index)
route[k2 & k2ok] = "km"
score("k>=2 minfrac>=0.2 -> KMeans (no SNV check); else SNV split -> SNV; else phase", route)

# where the SNV partition loses against the phasing
t = d[d.truth & d.snv_rec.notna()].copy()
t["snv_lost"] = ~t.snv_rec.astype(bool) & t.ph_rec.astype(bool)
t["snv_gained"] = t.snv_rec.astype(bool) & ~t.ph_rec.astype(bool)
t["stratum"] = t.cat + t.trf.map({True: " TRF", False: " non-TRF"})
print("\nSNV partition vs phasing by stratum (containers with an SNV split and truth):")
print(t.groupby("stratum").agg(n=("crID", "size"), snv_rec=("snv_rec", "mean"), ph_rec=("ph_rec", "mean"),
                               lost=("snv_lost", "sum"), gained=("snv_gained", "sum"),
                               sites_med=("snv_sites", "median")).round(3).to_string())
print("\nSNV lost containers: sites q10/q50/q90", q(t[t.snv_lost].snv_sites), " rest", q(t[~t.snv_lost].snv_sites))
