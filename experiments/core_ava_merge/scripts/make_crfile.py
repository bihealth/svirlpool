"""Candidate-regions file with the crIDs of every pattern on the given chromosomes."""

import pickle
import sys
from collections import defaultdict

comps = pickle.load(open("/home/mayv_c/development/svp_mendel/proto/all_svComposites.pkl", "rb"))
chroms = set(sys.argv[2:])
crs = defaultdict(set)
for s in comps:
    for p in s.svPatterns:
        if p.chr in chroms:
            crs[p.samplename].add(int(p.consensusID.split(".")[0]))
with open(sys.argv[1], "w") as f:
    for sample, ids in sorted(crs.items()):
        f.write(f"{sample}\t{','.join(map(str, sorted(ids)))}\n")
print({k: len(v) for k, v in crs.items()})
