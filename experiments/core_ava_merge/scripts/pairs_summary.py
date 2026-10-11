"""Summarize the judged pairs (pairs.<chr>.jsonl) of a harness or sv-calling run.

usage: pairs_summary.py DIR
"""

import json
import sys
from collections import Counter
from pathlib import Path


def sample_set(s):
    inside = s[s.index("[") + 1 : -1]
    return {x.split(":")[0] for x in inside.split(",")}


def main():
    d = Path(sys.argv[1])
    rows = [json.loads(line) for f in sorted(d.glob("pairs.*.jsonl")) for line in open(f)]
    print(len(rows), "pairs")
    c = Counter()
    ratio_added = Counter()
    ratio_removed = Counter()
    for r in rows:
        sa, sb = sample_set(r["a"]), sample_set(r["b"])
        rel = "same_sample" if sa == sb else "cross_sample"
        ev = r["evaluable"]
        old = r["old"]
        la = r["link_aln"] if ev else None
        c[(rel, "evaluable" if ev else "not_evaluable", f"old={old}", f"aln={la}")] += 1
        big = max(r["size_a"], r["size_b"], 1)
        rat = round(min(r["size_a"], r["size_b"]) / big, 1)
        if ev and la and not old:
            ratio_added[(rel, rat)] += 1
        if ev and old and not la:
            why = "distinct" if r.get("distinct") else ("x_same_no_overlap" if (r.get("x_same") or r.get("y_same")) else "no_overlap")
            ratio_removed[(rel, why)] += 1
    for k, v in sorted(c.items()):
        print(*k, v, sep="\t")
    print("added links by size ratio (min/max):")
    for k, v in sorted(ratio_added.items()):
        print(*k, v, sep="\t")
    print("removed links by reason:")
    for k, v in sorted(ratio_removed.items()):
        print(*k, v, sep="\t")


if __name__ == "__main__":
    main()
