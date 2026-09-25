"""F1 (truvari refined) and consensus-stage time of svp_improvements variants.

usage: e2e_table.py <variant> [...]
"""

import glob
import json
import re
import sys

R = "/home/mayv_c/development/svp_improvements/results"
print(
    "variant\tV5 all\tV5 non-TRF\tT2TQ100 all\tT2TQ100 non-TRF\tGT conc. V5\t"
    "consensus s (sum over batches)\tlongest batch s"
)
for v in sys.argv[1:]:
    row = [v]
    for ts in ("V5", "T2TQ100"):
        for rs in ("all", "non_trf"):
            with open(
                f"{R}/{v}/20x/truvari/{ts}/{rs}/refine.variant_summary.json"
            ) as f:
                row.append(f"{json.load(f)['f1']:.4f}")
    with open(f"{R}/{v}/20x/truvari/V5/all/summary.json") as f:
        row.append(f"{json.load(f).get('gt_concordance', float('nan')):.3f}")
    secs = []
    for b in glob.glob(
        f"{R}/{v}/20x/HG002/work/benchmarks/consensus/consensus.batch_*.txt"
    ):
        if not re.search(r"batch_\d+", b):
            continue
        lines = open(b).read().strip().splitlines()
        secs.append(
            float(
                dict(zip(lines[0].split("\t"), lines[-1].split("\t"), strict=True))["s"]
            )
        )
    # variants with db_from reuse another variant's consensus
    row += [f"{sum(secs):.0f}", f"{max(secs):.0f}"] if secs else ["-", "-"]
    print("\t".join(row))
