"""Are two svp_improvements variants' results the same?

Per sample: consensus batch outputs byte-identical (and which containers
differ), VCF records identical after dropping IDs and sorting; then the family
VCF the same way.

usage: compare_arms.py <variant A> <variant B>
"""

import glob
import gzip
import json
import os
import sys

R = "/home/mayv_c/development/svp_improvements/results"
SAMPLES = ("HG002", "HG003", "HG004")
a, b = sys.argv[1:3]


def vcf_records(path):
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            x[2] = "."
            out.append("\t".join(x))
    return sorted(out)


def containers(path):
    out = {}
    for line in open(path):
        cds = json.loads(line).get("consensus_dicts") or {}
        # consensus IDs are "<min crID>.<cluster>": the container is the prefix
        key = min((int(k.split(".")[0]) for k in cds), default=-1)
        out[str(key)] = out.get(str(key), "") + line
    return out


for smp in SAMPLES:
    fa = sorted(glob.glob(f"{R}/{a}/20x/{smp}/work/consensus/*/consensus.batch_*.jsonl"))
    same = differ = 0
    diff_containers = []
    for f in fa:
        g = f.replace(f"/{a}/", f"/{b}/")
        if not os.path.exists(g):
            differ += 1
            continue
        if open(f, "rb").read() == open(g, "rb").read():
            same += 1
            continue
        differ += 1
        ca, cb = containers(f), containers(g)
        diff_containers += [k for k in ca.keys() | cb.keys() if ca.get(k) != cb.get(k)]
    va = vcf_records(f"{R}/{a}/20x/{smp}/variants.vcf.gz")
    vb = vcf_records(f"{R}/{b}/20x/{smp}/variants.vcf.gz")
    sa, sb = set(va), set(vb)
    print(f"{smp}: consensus batches identical {same}/{same + differ}; "
          f"differing containers {len(diff_containers)} {sorted(diff_containers)[:10]}")
    print(f"       VCF records {len(va)} vs {len(vb)}: identical={va == vb}, "
          f"only {a} {len(sa - sb)}, only {b} {len(sb - sa)}")
fa = vcf_records(f"{R}/{a}/20x/family.vcf.gz")
fb = vcf_records(f"{R}/{b}/20x/family.vcf.gz")
print(f"family VCF records {len(fa)} vs {len(fb)}: identical={fa == fb}, "
      f"only {a} {len(set(fa) - set(fb))}, only {b} {len(set(fb) - set(fa))}")
