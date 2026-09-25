"""TP / FP / FN (truvari refined) of svp_improvements variants per stratum.

usage: counts.py <variant> [...]
"""

import json
import sys

R = "/home/mayv_c/development/svp_improvements/results"
print("variant\tstratum\tTP-base\tTP-comp\tFP\tFN\tprecision\trecall\tF1")
for v in sys.argv[1:]:
    for ts in ("V5", "T2TQ100"):
        for rs in ("all", "non_trf"):
            with open(f"{R}/{v}/20x/truvari/{ts}/{rs}/refine.variant_summary.json") as f:
                s = json.load(f)
            print(
                f"{v}\t{ts}/{rs}\t{s['TP-base']}\t{s['TP-comp']}\t{s['FP']}\t{s['FN']}\t"
                f"{s['precision']:.4f}\t{s['recall']:.4f}\t{s['f1']:.4f}"
            )
