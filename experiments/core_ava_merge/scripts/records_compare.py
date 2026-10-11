"""Compare old vs core-alignment merge records from harness.py (records.pkl).

usage: records_compare.py DIR
Only records with size >= 50 (the call's min size) count.
"""

import pickle
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


def groups(rec):
    g = defaultdict(int)
    for s, c, t, rs, re_, size, nreads in rec[1]:
        g[(s, c)] += size
    return g


def cross_pairs(recs):
    out = set()
    for r in recs:
        pats = [p[:5] for p in r[1]]
        for a, b in combinations(sorted(pats), 2):
            if a[0] != b[0]:
                out.add((a, b))
    return out


def describe(recs, name):
    recs = [r for r in recs if r[2] >= 50]
    c = Counter()
    for r in recs:
        g = groups(r)
        samples = {s for s, _ in g}
        c["records"] += 1
        c[f"n_samples={len(samples)}"] += 1
        if len(g) > 1:
            sizes = list(g.values())
            if min(sizes) / max(sizes) < 0.7:
                c["size_discordant_groups(>30%)"] += 1
        child = [v for (s, _), v in g.items() if s == "HG002"]
        if len(child) > 1 and min(child) / max(child) < 0.9:
            c["child_two_alleles(>10%)"] += 1
        same_cons_multi = Counter((p[0], p[1]) for p in r[1])
        if any(v > 1 for v in same_cons_multi.values()):
            c["records_with_pieces_of_one_consensus"] += 1
    print(name, dict(sorted(c.items())))
    return recs


def main():
    d = Path(sys.argv[1])
    old, new = pickle.load(open(d / "records.pkl", "rb"))
    o = describe(old, "old")
    n = describe(new, "new")
    po, pn = cross_pairs(o), cross_pairs(n)
    print("cross-sample pattern pairs co-merged: old", len(po), "new", len(pn), "added", len(pn - po), "removed", len(po - pn))
    # record-level: patterns whose record changed
    def rec_of(recs):
        m = {}
        for i, r in enumerate(recs):
            for p in r[1]:
                m[p[:5]] = i
        return m
    ro, rn = rec_of(old), rec_of(new)
    changed = sum(1 for k in ro if k in rn and frozenset(p[:5] for p in old[ro[k]][1]) != frozenset(p[:5] for p in new[rn[k]][1]))
    print("patterns whose record composition changed:", changed, "of", len(ro))


if __name__ == "__main__":
    main()
