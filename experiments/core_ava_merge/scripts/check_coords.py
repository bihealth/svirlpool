"""Check that pattern read_start/read_end are forward positions on the padded consensus."""

import pickle
from collections import Counter
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam

P = Path("/home/mayv_c/development/svp_mendel/proto")
comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
cores = pickle.load(open(P / "cores.pkl", "rb"))
c = Counter()
dels = Counter()
for svc in comps:
    for p in svc.svPatterns:
        key = (p.samplename, p.consensusID)
        if key not in cores:
            c["no_core"] += 1
            continue
        cs, core = cores[key]
        if isinstance(p, SVpatterns.SVpatternInsertion):
            seq = p.get_sequence()
            if seq is None:
                c["ins_noseq"] += 1
                continue
            s = p.read_start - cs
            e = p.read_end - cs
            sub = core[s:e]
            rev = p.SVprimitives[0].aln_is_reverse
            if rev:
                sub = cam.revcomp(sub)
            c[("ins", rev, sub == seq.upper(), len(p.SVprimitives))] += 1
        elif isinstance(p, SVpatterns.SVpatternDeletion):
            d = p.read_end - p.read_start
            dels[min(d, 5)] += 1
            c[("del_inside_core", 0 <= p.read_start - cs <= len(core))] += 1
print(c)
print("DEL read_end - read_start:", dels)
n_multi = Counter(len({(p.samplename, p.consensusID) for p in svc.svPatterns}) for svc in comps)
print("consensuses per composite", n_multi)
print("patterns per composite", Counter(min(len(svc.svPatterns), 5) for svc in comps))
strand = Counter(cam.composite_strand(svc) for svc in comps)
print("strand", strand)
