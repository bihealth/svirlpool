import json
import sys

for arm in sys.argv[1:]:
    for t in ("T2TQ100", "V5"):
        d = json.load(open(f"/home/mayv_c/development/svp_mendel/arms/{arm}/truvari_{t}/summary.json"))
        print(arm, t, {k: d[k] for k in ("TP-base", "TP-comp", "FP", "FN")})
