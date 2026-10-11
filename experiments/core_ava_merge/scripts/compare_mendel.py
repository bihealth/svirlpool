"""Record-level Mendel transitions between the base arm and an arm.

usage: compare_mendel.py ARM

Each PASS >= 50 bp informative record of one call is linked to the records of the
other call that share a consensus ID, have the same SVTYPE and lie within 2 kb.
- fixed: base inconsistent record whose arm counterparts (via its HG002
  consensus) are all consistent;
- new violation: arm inconsistent record whose base counterparts were all
  consistent (or that has no base counterpart).
Also: size spread of the arm's merged records (CONSENSUSIDs from >= 2 samples).
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from eval_labels import annotate, read_vcf  # noqa: E402

ROOT = Path("/home/mayv_c/development/svp_mendel")


def index(recs):
    by_cons = defaultdict(list)
    for r in recs:
        for c in r["cons"]:
            by_cons[c].append(r)
    return by_cons


def counterparts(r, by_cons, only_child=False):
    out = {}
    for c in r["cons"]:
        if only_child and not c.startswith("HG002:"):
            continue
        for o in by_cons.get(c, []):
            if o["svtype"] == r["svtype"] and o["chrom"] == r["chrom"] and abs(o["pos"] - r["pos"]) <= 2000:
                out[id(o)] = o
    return list(out.values())


def find(r, by_key, by_cons):
    """The other call's record: same CHROM/POS/SVTYPE/SVLEN, else the nearest
    record sharing a consensus ID (same type, <= 2 kb)."""
    exact = by_key.get((r["chrom"], r["pos"], r["svtype"], r["svlen"]))
    if exact:
        return exact[0]
    cps = counterparts(r, by_cons)
    return min(cps, key=lambda o: abs(o["pos"] - r["pos"])) if cps else None


def status(r):
    if r["filter"] != "PASS" or r["svlen"] < 50:
        return "filtered"
    return r["mndl"] or "."


def main():
    arm = sys.argv[1]
    base = read_vcf(annotate("base_c8804b1"))
    new = read_vcf(annotate(arm))
    bi, ni = index(base), index(new)
    bk, nk = defaultdict(list), defaultdict(list)
    for r in base:
        bk[(r["chrom"], r["pos"], r["svtype"], r["svlen"])].append(r)
    for r in new:
        nk[(r["chrom"], r["pos"], r["svtype"], r["svlen"])].append(r)
    c = Counter()
    for r in base:
        if status(r) != "inconsistent":
            continue
        o = find(r, nk, ni)
        c["base_inconsistent -> arm " + (status(o) if o else "none")] += 1
    for r in new:
        if status(r) != "inconsistent":
            continue
        o = find(r, bk, bi)
        c["arm_inconsistent <- base " + (status(o) if o else "none")] += 1
    # merged records: size spread between the samples' groups is not in the VCF;
    # count records joining samples
    for name, recs in (("base", base), ("arm", new)):
        n = Counter()
        for r in recs:
            if status(r) == "filtered":
                continue
            s = {x.split(":")[0] for x in r["cons"]}
            n[len(s)] += 1
            child = [x for x in r["cons"] if x.startswith("HG002:")]
            if len({x.split(".")[0] for x in child}) < len(child):
                n["child_sibling_consensuses"] += 1
        print(name, "PASS>=50 records by #samples in CONSENSUSIDs:", dict(sorted(n.items(), key=str)))
    for k, v in sorted(c.items()):
        print(k, v, sep="\t")


if __name__ == "__main__":
    main()
