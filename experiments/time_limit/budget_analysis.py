"""Per-container comparison of consensus budget variants (the cluster tuning
experiments budget10 / budget10b / budget15; README section 7).

Reads the tables container_cost.py writes, <costdir>/<variant>.<sample>.tsv
(HG002 with truth), and prints per variant: summed container wall and
thread-seconds, escalations and their tools, containers dropped at the time
limit / at --max-assembly-bp, containers that retried a larger LAST -m; whether
two gate variants drop the same containers; the value (HG002 TP / FP calls,
Q100 representation) and cost of the gated containers priced in an uncapped
variant; a sweep of gate thresholds on the largest assembly of each
container; and the -m 10,50 saving by assembly size.

usage: budget_analysis.py <costdir> --variants a,b,... --uncapped V --m10 V
          [--gates V1,V2] [--default V --tmp V]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

SAMPLES = ["HG002", "HG003", "HG004"]
FILL = ("time_dropped", "size_dropped", "m_retries", "max_call_bp", "escalations")


def load(costdir: str, variants: list[str]) -> dict:
    t = {}
    for v in variants:
        for s in SAMPLES:
            d = pd.read_csv(f"{costdir}/{v}.{s}.tsv", sep="\t").set_index("crID")
            for c in FILL:
                d[c] = d[c].fillna(0) if c in d else 0
            t[v, s] = d
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("costdir")
    ap.add_argument("--variants", required=True)
    ap.add_argument("--uncapped", required=True, help="variant without any container limit, -m 50")
    ap.add_argument("--m10", required=True, help="variant with -m 10,50 and no limit")
    ap.add_argument("--gates", default="", help="two variants with the same gate")
    ap.add_argument("--default", default="")
    ap.add_argument("--tmp", default="", help="the default with node-local temporary files")
    a = ap.parse_args()
    variants = a.variants.split(",")
    T = load(a.costdir, variants)
    pd.set_option("display.width", 250)
    fmt = lambda x: f"{x:.2f}"  # noqa: E731

    print("== per variant (trio sums)")
    rows = []
    for v in variants:
        x = pd.concat([T[v, s] for s in SAMPLES])
        rows.append({"variant": v, "wall_h": x.seconds.sum() / 3600, "thread_h": x.thread_seconds.sum() / 3600,
                     "median_s": x.seconds.median(), "p99_s": x.seconds.quantile(.99), "max_s": x.seconds.max(),
                     "escalated": int((x.escalations > 0).sum()), "time_drop": int(x.time_dropped.sum()),
                     "size_drop": int(x.size_dropped.sum()), "m_retry": int((x.m_retries > 0).sum()),
                     "max_call_kb": x.max_call_bp.max() / 1e3})
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))

    print("\n== escalation tools (trio)")
    for v in variants:
        x = pd.concat([T[v, s] for s in SAMPLES])
        print(v, x.esc_tools.dropna().str.split(";").explode().value_counts().to_dict())
    if a.default and a.tmp:
        for s in SAMPLES:
            j = T[a.default, s][["escalations", "sec_level0"]].join(
                T[a.tmp, s][["escalations"]], lsuffix="_def", rsuffix="_tmp", how="inner")
            only = j[(j.escalations_tmp > 0) & (j.escalations_def == 0)]
            print(f"{s}: escalated only with {a.tmp}: {len(only)} (level-1 s in {a.default}: "
                  f"median {only.sec_level0.median():.1f}); only in {a.default}: "
                  f"{int(((j.escalations_def > 0) & (j.escalations_tmp == 0)).sum())}")

    if a.gates:
        g1, g2 = a.gates.split(",")
        print(f"\n== gate determinism: containers dropped at --max-assembly-bp, {g1} vs {g2}")
        for s in SAMPLES:
            x = set(T[g1, s].index[T[g1, s].size_dropped > 0])
            y = set(T[g2, s].index[T[g2, s].size_dropped > 0])
            print(f"{s}: {len(x)} / {len(y)}, common {len(x & y)}, only {g1} {sorted(x - y)}, only {g2} {sorted(y - x)}")
        ref, m10 = T[a.uncapped, "HG002"], T[a.m10, "HG002"]
        ids = T[g1, "HG002"].index[T[g1, "HG002"].size_dropped > 0]
        x = ref.loc[ref.index.intersection(ids)]
        print(f"\n== {g1} HG002 gated containers, priced in {a.uncapped}: {len(x)}: TP {x.tp.sum():.0f} "
              f"FP {x.fp.sum():.0f}, repr {x.repr.mean():.3f}; thread-h {x.thread_seconds.sum() / 3600:.2f} "
              f"of {ref.thread_seconds.sum() / 3600:.2f} ({a.uncapped}), "
              f"{m10.loc[m10.index.intersection(ids)].thread_seconds.sum() / 3600:.2f} of "
              f"{m10.thread_seconds.sum() / 3600:.2f} ({a.m10})")

    ref, m10 = T[a.uncapped, "HG002"], T[a.m10, "HG002"]
    print(f"\n== gate sweep on the largest assembly per container ({a.uncapped} inputs, HG002 truth)")
    rows = []
    for kb in (150, 200, 250, 300, 400, 500, 700):
        sel = ref.max_call_bp > kb * 1e3
        r = {"gate_kb": kb, "HG002": int(sel.sum()), "TP_lost": ref.tp[sel].sum(), "FP_removed": ref.fp[sel].sum(),
             "thread_h_saved_uncapped": ref.thread_seconds[sel].sum() / 3600,
             "thread_h_saved_m10": m10.loc[m10.index.intersection(ref.index[sel])].thread_seconds.sum() / 3600}
        for s in ("HG003", "HG004"):
            r[s] = int((T[a.uncapped, s].max_call_bp > kb * 1e3).sum())
        rows.append(r)
    print(pd.DataFrame(rows).to_string(index=False, float_format=fmt))

    print(f"\n== -m 10,50 by the largest assembly (HG002 thread-s, {a.uncapped} vs {a.m10})")
    j = ref[["thread_seconds", "max_call_bp"]].join(m10[["thread_seconds", "m_retries"]], rsuffix="_m10", how="inner")
    j["bin"] = pd.cut(j.max_call_bp, [0, 5e4, 1e5, 2e5, 3e5, 5e5, np.inf])
    print(j.groupby("bin", observed=True).agg(
        n=("thread_seconds", "size"), m50=("thread_seconds", "sum"), m10=("thread_seconds_m10", "sum"),
        retries=("m_retries", "sum")).assign(ratio=lambda g: g.m10 / g.m50).to_string(float_format=lambda x: f"{x:.1f}"))


if __name__ == "__main__":
    main()
