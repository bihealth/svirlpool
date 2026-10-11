"""A1 records whose child composite links to the transmitting parent at the pair
level (reachable.tsv best=link) but whose arm record is not consistent."""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from eval_labels import annotate, read_vcf  # noqa: E402

ROOT = Path("/home/mayv_c/development/svp_mendel")
arm = sys.argv[1]
reach = {}
for line in open(ROOT / "proto/reachable.tsv"):
    t = line.rstrip("\n").split("\t")
    reach[t[1]] = t[7]
labels = {r["idx"]: r for r in csv.DictReader(open(ROOT / "tr_substring/out/record_classes.tsv"), delimiter="\t")}
base = read_vcf(ROOT / "tr_substring/data/family.mndl.vcf.gz")
bk = {(r["chrom"], r["pos"], r["svtype"], r["svlen"]): r for r in base}
recs = read_vcf(annotate(arm))
by_key = defaultdict(list)
by_cons = defaultdict(list)
for r in recs:
    by_key[(r["chrom"], r["pos"], r["svtype"], r["svlen"])].append(r)
    for c in r["cons"]:
        by_cons[c].append(r)
c = Counter()
examples = []
for idx, best in reach.items():
    lab = labels[idx]
    if lab["klass"] != "A1" or lab["cls"] != "a":
        continue
    key = (lab["chrom"], int(lab["pos"]), lab["svtype"], int(lab["svlen"]))
    b = bk[key]
    child = [x for x in b["cons"] if x.startswith("HG002:")]
    cand = [r for r in by_key.get(key, []) if set(r["cons"]) & set(child)]
    if not cand:
        cand = sorted((r for x in child for r in by_cons.get(x, []) if r["svtype"] == lab["svtype"]
                       and r["chrom"] == lab["chrom"] and abs(r["pos"] - key[1]) <= 2000), key=lambda r: abs(r["pos"] - key[1]))
    if not cand:
        c[(best, "lost")] += 1
        continue
    r = cand[0]
    tp = lab["t_parent"]
    joined = any(x.startswith(tp + ":") for x in r["cons"])
    ok = r["mndl"] == "consistent" and r["filter"] == "PASS"
    # any arm record carrying the child consensus that is consistent and joined?
    other_ok = any(o["mndl"] == "consistent" and any(x.startswith(tp + ":") for x in o["cons"])
                   for x in child for o in by_cons.get(x, []) if abs(o["pos"] - key[1]) <= 2000)
    state = ("consistent" if ok else "joined_inconsistent" if joined else "other_piece_ok" if other_ok else "not_joined")
    c[(best, state)] += 1
    if best == "link" and state in ("joined_inconsistent", "not_joined"):
        examples.append((state, lab["chrom"], lab["pos"], lab["svtype"], lab["svlen"], lab["gt"], tp, r["pos"], r["svlen"],
                         [r["gts"][s] for s in ("HG002", "HG003", "HG004")], r["cons"][:5]))
for k, v in sorted(c.items()):
    print(*k, v, sep="\t")
for e in examples[:25]:
    print(*e)
