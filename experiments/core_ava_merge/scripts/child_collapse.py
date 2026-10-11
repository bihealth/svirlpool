"""Records (harness, new merge) with two HG002 consensuses whose summed sizes differ > 10%.

For each: the HG002 consensus groups, and the pair decisions between composites of
the two child consensuses (from the pair log) plus whether old records had them
together.
usage: child_collapse.py DIR [n_examples]
"""

import json
import pickle
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

d = Path(sys.argv[1])
n_ex = int(sys.argv[2]) if len(sys.argv) > 2 else 12
old, new = pickle.load(open(d / "records.pkl", "rb"))


def child_groups(rec):
    g = defaultdict(int)
    for s, c, t, rs, re_, size, nreads in rec[1]:
        if s == "HG002":
            g[c] += size
    return g


def flagged(recs):
    out = []
    for r in recs:
        if r[2] < 50:
            continue
        g = child_groups(r)
        if len(g) > 1 and min(g.values()) / max(g.values()) < 0.9:
            out.append(r)
    return out


fo, fn = flagged(old), flagged(new)
key = lambda r: frozenset(p[:5] for p in r[1])  # noqa: E731
so = {key(r) for r in fo}
only_new = [r for r in fn if key(r) not in so]
print("flagged old", len(fo), "new", len(fn), "new-only", len(only_new))
# classify new-only: are the two child consensuses siblings (same crID) or different CRs?
c = Counter()
for r in only_new:
    g = child_groups(r)
    crs = {cid.split(".")[0] for cid in g}
    c["siblings" if len(crs) == 1 else "different_crIDs"] += 1
print(c)
import json
json.dump([[r[0], r[2], sorted({(s, c) for s, c, *_ in r[1]}), [list(p) for p in r[1]]] for r in only_new], open(d / "child_collapse.json", "w"))
random.seed(2)
for r in random.sample(only_new, min(n_ex, len(only_new))):
    g = child_groups(r)
    print(r[0], r[2], dict(g), sorted({(s, c) for s, c, *_ in r[1]})[:8])
