"""Per container: consensuses and phasing metadata of two svp_improvements runs.

usage: compare_consensus.py <variant_a> <variant_b>
"""

import glob
import json
import sys
from collections import Counter

R = "/home/mayv_c/development/svp_improvements/results"


def load(v):
    out = {}
    for f in glob.glob(f"{R}/{v}/20x/HG002/work/consensus/*/consensus.batch_*.jsonl"):
        for line in open(f):
            d = json.loads(line)
            cons = d["consensus_dicts"]
            if not cons:
                continue
            crIDs = next(iter(cons.values()))["crIDs"]
            key = min(crIDs)
            metas = [c.get("clustering_meta_data") or {} for c in cons.values()]
            out[key] = {
                "n": len(cons),
                "method": metas[0].get("method"),
                "status": metas[0].get("phasing_status"),
                "n_reads": sorted(
                    (len(c["intervals_cutread_alignments"]) for c in cons.values()),
                    reverse=True,
                ),
                "lens": sorted(len(c["consensus_sequence"]) for c in cons.values()),
            }
    return out


a, b = load(sys.argv[1]), load(sys.argv[2])
keys = sorted(a.keys() | b.keys())
print("containers", len(a), len(b))
print("n consensus a", Counter(x["n"] for x in a.values()))
print("n consensus b", Counter(x["n"] for x in b.values()))
diff_n = [k for k in keys if (a.get(k) or {}).get("n") != (b.get(k) or {}).get("n")]
print("containers with a different consensus count:", len(diff_n))
print(
    "transitions",
    Counter(((a.get(k) or {}).get("n"), (b.get(k) or {}).get("n")) for k in diff_n),
)
print(
    "method a->b",
    Counter(
        ((a.get(k) or {}).get("method"), (b.get(k) or {}).get("method")) for k in keys
    ).most_common(8),
)
same_n = [k for k in keys if k in a and k in b and a[k]["n"] == b[k]["n"]]
print(
    "same count, different read counts:",
    sum(a[k]["n_reads"] != b[k]["n_reads"] for k in same_n),
    "different consensus lengths:",
    sum(a[k]["lens"] != b[k]["lens"] for k in same_n),
)
