"""Summarise kmeans_gate.py: how often the KMeans gate accepts, the phasing time
it would save, and the SV-allele accuracy of KMeans vs phasing on the
containers it accepts.

usage: kmeans_gate_summary.py <kmeans_gate.tsv>
"""

import sys

import pandas as pd

pd.set_option("display.width", 200)
d = pd.read_csv(sys.argv[1], sep="\t")
run = d[~d.skipped.astype(bool)].copy()
print(f"containers {len(d)}, skipped (CN > 4) {int(d.skipped.sum())}, processed {len(run)}")
acc = run.gate_k > 0
print(f"gate accepts {acc.sum()} ({acc.mean():.1%}); k: {run.gate_k.value_counts().sort_index().to_dict()}")
print(
    f"time (s, 1 thread): summed indels {run.t_indels.sum():.0f}, gate {run.t_gate.sum():.0f}, "
    f"phasing {run.t_phase.sum():.0f} (accepted {run.t_phase[acc].sum():.0f}, "
    f"rejected {run.t_phase[~acc].sum():.0f})"
)
print("phasing status on accepted:", run[acc].phase_status.value_counts().to_dict())

ev = run[acc & run.bench.eq(True) & run.km_recovered.notna() & run.ph_recovered.notna()].copy()
for c in ("km_recovered", "ph_recovered"):
    ev[c] = ev[c].astype(bool)
print(f"\naccepted containers with truth (T2TQ100 bench, >= 4 trio-labelled reads): {len(ev)}")


def table(g):
    return pd.Series({
        "n": len(g),
        "km_recovered": g.km_recovered.mean(), "ph_recovered": g.ph_recovered.mean(),
        "km_purity": g.km_purity.mean(), "ph_purity": g.ph_purity.mean(),
        "km_assigned": g.km_assigned.mean(), "ph_assigned": g.ph_assigned.mean(),
        "km_worse": int((~g.km_recovered & g.ph_recovered).sum()),
        "km_better": int((g.km_recovered & ~g.ph_recovered).sum()),
    })


ev["stratum"] = ev.cat + ev.trf.map({True: " TRF", False: " non-TRF"})
print(table(ev).round(3).to_string())
print()
print(ev.groupby("stratum").apply(table, include_groups=False).round(3).to_string())
print()
print(ev.groupby("gate_k").apply(table, include_groups=False).round(3).to_string())
worse = ev[~ev.km_recovered & ev.ph_recovered]
print(f"\nKMeans loses an allele the phasing recovers ({len(worse)}):")
print(worse[["crID", "chr", "start", "cat", "trf", "n_reads", "gate_k", "phase_n_alleles",
             "km_purity", "ph_purity"]].to_string(index=False))
