"""Summarise ref_vs_ava_eval.py: arm A (all-vs-all sites) vs arm B (reference
sites), paired per container.

Read level (raw phasing groups vs trio read labels): containers where both
arms assign >= 2 labelled reads. Allele level (the partition production would
assemble, vs T2TQ100): containers whose truth locus matches and that are fully
benchmarked. Strata: TR = >= 50% of the container's SV signals in a tandem
repeat; label coverage = trio-labelled reads / reads (the trio labels come
from reference alignments, so low coverage marks reads the reference handles
badly).

usage: ref_vs_ava_summary.py <ref_vs_ava_eval.tsv>
"""

import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 200)
d = pd.read_csv(sys.argv[1], sep="\t")
d = d[~d.skipped.astype(bool)].copy()
d["TR"] = d.sig_trf_frac >= 0.5
d["label_cov"] = d.trio_total / d.n_reads
d["cov_cls"] = pd.cut(d.label_cov, [-0.01, 0.6, 0.8, 1.01], labels=["<0.6", "0.6-0.8", ">=0.8"])

print(f"containers {len(d)}")
print(f"time  A {d.A_secs.sum():.0f} s  B {d.B_secs.sum():.0f} s  "
      f"(per container median A {d.A_secs.median():.2f} s, B {d.B_secs.median() * 1000:.0f} ms)")
for a in "AB":
    print(f"{a}: status {d[f'{a}_status'].value_counts().to_dict()}  "
          f"alleles {dict(sorted(d[f'{a}_n_alleles'].value_counts().items()))}  "
          f"median SNV sites {d[f'{a}_n_snv'].median():.0f}  SV sites {d[f'{a}_n_sv'].median():.0f}  "
          f"low-quality reads {d[f'{a}_n_lowq'].sum()}")
print("\nallele count, A (rows) vs B (columns):")
print(pd.crosstab(d.A_n_alleles, d.B_n_alleles).to_string())


def read_level(x, label):
    x = x[(x.A_trio_assigned >= 2) & (x.B_trio_assigned >= 2)]
    if len(x) == 0:
        return
    da = x.A_acc - x.B_acc
    print(f"  {label:28s} n {len(x):5d}  acc A {x.A_acc.mean():.4f}  B {x.B_acc.mean():.4f}  "
          f"perfect A {(x.A_acc == 1).mean():.3f}  B {(x.B_acc == 1).mean():.3f}  "
          f"assigned A {x.A_trio_assigned.sum() / x.trio_total.sum():.3f}  "
          f"B {x.B_trio_assigned.sum() / x.trio_total.sum():.3f}  "
          f"A better {(da > 1e-9).sum():4d}  B better {(da < -1e-9).sum():4d}")


print("\nread level (pair accuracy of assigned trio-labelled reads; paired containers):")
read_level(d, "all")
for tr in (False, True):
    read_level(d[d.TR == tr], "TR" if tr else "non-TR")
for c in d.cov_cls.cat.categories:
    read_level(d[d.cov_cls == c], f"label coverage {c}")

t = d[d.truth_locus_ok.eq(True) & d.bench.eq(True) & d.A_recovered.notna() & d.B_recovered.notna()].copy()
for a in "AB":
    t[f"{a}_rec"] = t[f"{a}_recovered"].astype(str).eq("True")
t["stratum"] = t.cat + np.where(t.trf.astype(str).eq("True"), " TRF", " non-TRF")


def allele_level(x, label):
    lost = (x.A_rec & ~x.B_rec).sum()
    gained = (~x.A_rec & x.B_rec).sum()
    print(f"  {label:22s} n {len(x):5d}  recovered A {x.A_rec.mean():.3f}  B {x.B_rec.mean():.3f}  "
          f"purity A {x.A_purity.mean():.4f}  B {x.B_purity.mean():.4f}  "
          f"A only {lost:3d}  B only {gained:3d}")


print("\nallele level (production partition vs T2TQ100; A only / B only = allele set recovered by one arm only):")
allele_level(t, "all")
for s, x in t.groupby("stratum"):
    allele_level(x, s)
