"""Tiered phasing: when can a container keep the cheap reference-site phasing
(arm B) instead of the all-vs-all phasing (arm A)?

A rule routes a container to B from B's own result and production-available
features (ref_route_features.py); the hybrid is scored against always-A:

  routed     share of containers kept on B
  saved      share of arm-A phasing time avoided
  bad/better containers where the hybrid's trio pair accuracy is > 0.05
             below/above always-A (read level, paired containers)
  lost/gained containers whose T2TQ100 allele set is recovered by A only /
             by B only, among the routed (allele level)

Rules: B finds n alleles in AL, B partition discordance <= DISC, mean allele
imbalance at kept SNV columns <= AFD, smallest B group share >= MINF, median
read divergence to the reference <= DIV. Selection: the most time saved with
lost <= gained and (bad - better) <= BUDGET x paired containers.

Checks: (1) rules fitted on one set, scored on another (disjoint regions);
(2) leave-one-chromosome-out within a set (fit on 21 autosomes, apply to the
22nd).

usage: ref_route_rules.py NAME=eval.tsv[.gz],features.tsv ... [--fit NAME] [--test NAME]
"""

import argparse
import itertools

import numpy as np
import pandas as pd

GRID = {
    "al": [(2,), (1, 2)],
    "disc": [0.02, 0.03, 0.04, 0.05, 0.07, 1.0],
    "afd": [0.15, 0.2, 0.25, 1.0],
    "minf": [0.0, 0.2, 0.3],
    "div": [0.02, 0.03, 1.0],
}
COMBOS = list(itertools.product(*GRID.values()))
SIMPLE = ((2,), 0.02, 1.0, 0.3, 1.0)


def load(spec):
    name, files = spec.split("=", 1)
    ev, ft = files.split(",")
    d = pd.read_csv(ev, sep="\t")
    d = d[~d.skipped.astype(bool)]
    f = pd.read_csv(ft, sep="\t")
    m = d.merge(f, on="crID", suffixes=("", "_f"))
    m["A_rec_b"] = m.A_recovered.astype(str).eq("True")
    m["B_rec_b"] = m.B_recovered.astype(str).eq("True")
    m["has_allele"] = (m.truth_locus_ok.eq(True) & m.bench.eq(True)
                       & m.A_recovered.notna() & m.B_recovered.notna())
    m["paired"] = (m.A_trio_assigned >= 2) & (m.B_trio_assigned >= 2)
    m["set"] = name
    return m.reset_index(drop=True)


def rule(x, c):
    al, disc, afd, minf, div = c
    return (x.B_n_alleles_f.isin(al) & (x.discord.fillna(0) <= disc)
            & (x.af_dev.fillna(0) <= afd) & (x.B_min_frac >= minf)
            & (x.div_med <= div)).to_numpy()


def evaluate(x, route):
    p = x.paired.to_numpy()
    xa, xb = x.A_acc.to_numpy(), x.B_acc.to_numpy()
    bad = ((xa - xb) > 0.05) & route & p
    better = ((xb - xa) > 0.05) & route & p
    acc = np.where(route, xb, xa)[p]
    h = x.has_allele.to_numpy()
    ra, rb = x.A_rec_b.to_numpy(), x.B_rec_b.to_numpy()
    return {
        "n": len(x), "paired": int(p.sum()), "routed": route.mean(),
        "saved": x.A_secs.to_numpy()[route].sum() / max(x.A_secs.sum(), 1e-9),
        "acc": acc.mean(), "accA": xa[p].mean(), "bad": int(bad.sum()), "better": int(better.sum()),
        "lost": int((ra & ~rb & route & h).sum()), "gained": int((~ra & rb & route & h).sum()),
    }


def select(x, rules, budget):
    best = None
    for c, r in rules.items():
        e = evaluate(x, r)
        if e["lost"] <= e["gained"] and e["bad"] - e["better"] <= budget * e["paired"]:
            if best is None or e["saved"] > best[0]:
                best = (e["saved"], c)
    return None if best is None else best[1]


def fmt(e):
    return (f"n {e['n']:5d}  routed {e['routed']:.3f}  saved {e['saved']:.3f}  "
            f"acc {e['acc']:.4f} (A {e['accA']:.4f})  bad {e['bad']:3d} better {e['better']:3d}  "
            f"alleles lost {e['lost']:2d} gained {e['gained']:2d}")


def loco(x, budget):
    rules = {c: rule(x, c) for c in COMBOS}
    route = np.zeros(len(x), bool)
    chosen = []
    chrs = x.chr.to_numpy()
    for ch in sorted(set(chrs)):
        tr, te = chrs != ch, chrs == ch
        c = select(x[tr], {k: v[tr] for k, v in rules.items()}, budget)
        chosen.append(c)
        if c is not None:
            route[te] = rules[c][te]
    return evaluate(x, route), pd.Series([str(c) for c in chosen]).value_counts()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("sets", nargs="+")
    p.add_argument("--fit")
    p.add_argument("--test")
    p.add_argument("--budgets", default="0,0.002,0.005")
    a = p.parse_args()
    sets = {s.split("=")[0]: load(s) for s in a.sets}
    if len(sets) > 1:
        sets["all"] = pd.concat(sets.values(), ignore_index=True)
    budgets = [float(b) for b in a.budgets.split(",")]
    for name, x in sets.items():
        print(f"\n== {name}: always B  {fmt(evaluate(x, np.ones(len(x), bool)))}")
        print(f"   simple rule {SIMPLE}: {fmt(evaluate(x, rule(x, SIMPLE)))}")
    if a.fit and a.test:
        f, t = sets[a.fit], sets[a.test]
        print(f"\n== fitted on {a.fit}, tested on {a.test}")
        for b in budgets:
            c = select(f, {k: rule(f, k) for k in COMBOS}, b)
            if c is None:
                print(f"  budget {b}: no rule")
                continue
            print(f"  budget {b}: rule {c}\n    fit  {fmt(evaluate(f, rule(f, c)))}\n    test {fmt(evaluate(t, rule(t, c)))}")
    for name, x in sets.items():
        print(f"\n== leave-one-chromosome-out on {name}")
        for b in budgets:
            e, chosen = loco(x, b)
            print(f"  budget {b}: {fmt(e)}")
            print(f"    rules chosen: {chosen.head(4).to_dict()}")


if __name__ == "__main__":
    main()
