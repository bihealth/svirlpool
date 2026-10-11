"""Classify an arm's new Mendel violations (arm inconsistent, base counterpart consistent)."""

import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from compare_mendel import find, index, status  # noqa: E402
from eval_labels import annotate, read_vcf  # noqa: E402

arm = sys.argv[1]
base = read_vcf(annotate("base_c8804b1"))
new = read_vcf(annotate(arm))
bi = index(base)
bk = defaultdict(list)
for r in base:
    bk[(r["chrom"], r["pos"], r["svtype"], r["svlen"])].append(r)
c = Counter()
ex = defaultdict(list)
for r in new:
    if status(r) != "inconsistent":
        continue
    o = find(r, bk, bi)
    if not o or status(o) != "consistent":
        continue
    sa = {x.split(":")[0] for x in r["cons"]}
    sb = {x.split(":")[0] for x in o["cons"]}
    ca, cb = set(r["cons"]), set(o["cons"])
    if ca < cb:
        k = "split (arm record lost consensuses of the base record)"
    elif ca > cb:
        k = "joined (arm record gained consensuses)"
    elif ca == cb:
        k = "same consensuses (representative/size/GT changed)"
    else:
        k = "regrouped"
    k2 = "child_alt" if r["gts"]["HG002"] not in ("0/0", "./.", "0") else "child_ref"
    c[(k, k2)] += 1
    ex[k].append((r, o))
for k, v in sorted(c.items()):
    print(*k, v, sep="\t")
random.seed(1)
for k, rows in ex.items():
    print("==", k)
    for r, o in random.sample(rows, min(5, len(rows))):
        print("  arm ", r["chrom"], r["pos"], r["svtype"], r["svlen"], [r["gts"][s] for s in ("HG002", "HG003", "HG004")], r["cons"][:6])
        print("  base", o["chrom"], o["pos"], o["svtype"], o["svlen"], [o["gts"][s] for s in ("HG002", "HG003", "HG004")], o["cons"][:6])
