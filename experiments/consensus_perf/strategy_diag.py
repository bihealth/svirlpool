"""Where does --clustering-strategy fast lose? HG002 truvari bench FPs / FNs of
accurate and fast, by the route fast took in the container (fast_kmeans,
fast_snv, phased, phased_single). FPs are assigned to containers by their
CONSENSUSIDs, FNs by position (container regions +- 1 kb).

usage: strategy_diag.py [truthset=V5]
"""

import glob
import gzip
import json
import sys
from collections import Counter, defaultdict

R = "/home/mayv_c/development/svp_improvements/results"
TS = sys.argv[1] if len(sys.argv) > 1 else "V5"


def container_info(v):
    route, regions = {}, {}
    for f in glob.glob(f"{R}/{v}/20x/HG002/work/consensus/*/consensus.batch_*.jsonl"):
        for line in open(f):
            cds = json.loads(line).get("consensus_dicts") or {}
            for cid, c in cds.items():
                crid = int(cid.split(".")[0])
                route.setdefault(crid, (c.get("clustering_meta_data") or {}).get("method", "none"))
                regions[crid] = c["original_regions"]
    return route, regions


def records(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] != "#":
                x = line.split("\t", 8)
                yield x[0], int(x[1]), x[7]


route, regions = container_info("cs_fast")
by_chr = defaultdict(list)
for crid, regs in regions.items():
    for ch, s, e in regs:
        by_chr[ch].append((s - 1000, e + 1000, crid))


def locate(ch, pos):
    for s, e, crid in by_chr.get(ch, ()):
        if s <= pos <= e:
            return crid
    return None


print(f"containers by fast route: {dict(Counter(route.values()))}")
for rs in ("all", "non_trf"):
    print(f"\n{TS} {rs} (truvari bench, before refine):")
    tab = {}
    for v in ("cs_accurate", "cs_fast"):
        base = f"{R}/{v}/20x/truvari/{TS}/{rs}"
        fp, fn = Counter(), Counter()
        for _ch, _pos, info in records(f"{base}/fp.vcf.gz"):
            ids = [t for t in info.split(";") if t.startswith("CONSENSUSIDs=")]
            crids = {int(x.split(":")[1].split(".")[0]) for x in ids[0][13:].split(",")} if ids else set()
            fp[route.get(min(crids), "?") if crids else "?"] += 1
        for ch, pos, _info in records(f"{base}/fn.vcf.gz"):
            crid = locate(ch, pos)
            fn[route.get(crid, "outside containers") if crid is not None else "outside containers"] += 1
        tab[v] = (fp, fn)
    keys = sorted(set().union(*[set(fp) | set(fn) for fp, fn in tab.values()]))
    print(f"  {'fast route':22s} {'FP acc':>7s} {'FP fast':>8s} {'FN acc':>7s} {'FN fast':>8s}")
    for k in keys:
        a, f = tab["cs_accurate"], tab["cs_fast"]
        print(f"  {k:22s} {a[0][k]:7d} {f[0][k]:8d} {a[1][k]:7d} {f[1][k]:8d}")
