"""Which containers the FPs gained / lost between two variants come from, and
how the phasing ended there in each run.

usage: fp_containers.py <variant_a> <variant_b> [truthset/regionset]
"""

import glob
import gzip
import json
import re
import sys
from collections import defaultdict

R = "/home/mayv_c/development/svp_improvements/results"
va, vb = sys.argv[1], sys.argv[2]
stratum = sys.argv[3] if len(sys.argv) > 3 else "V5/non_trf"


def fps(v):
    out = set()
    with gzip.open(f"{R}/{v}/20x/truvari/{stratum}/fp.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t")
            out.add((x[0], int(x[1]), len(x[4]) - len(x[3])))
    return out


def regions(v):
    """container (min crID) -> list of (chr, start, end)"""
    out = {}
    for fn in glob.glob(f"{R}/{v}/20x/HG002/work/consensus/*/consensus.batch_*.jsonl"):
        for line in open(fn):
            cons = json.loads(line)["consensus_dicts"]
            if cons:
                c = next(iter(cons.values()))
                out[min(c["crIDs"])] = c["original_regions"]
    return out


def statuses(v):
    """container (representative crID) -> phasing log line"""
    out = {}
    for fn in glob.glob(f"{R}/{v}/20x/HG002/work/consensus/*/consensus.batch_*.log"):
        cur = None
        for line in open(fn):
            m = re.search(r"representative crID (\d+)\)", line)
            if m:
                cur = int(m.group(1))
            elif cur is not None and (
                "read phasing: status" in line or "all-vs-all alignment failed" in line
            ):
                s = re.search(r"status=\w+ alleles=\d+ sizes=\[[^\]]*\]", line)
                out[cur] = s.group(0) if s else "FAILED (timeout)"
    return out


reg = regions(vb)
sa, sb = statuses(va), statuses(vb)
fa, fb = fps(va), fps(vb)


def container_of(chrom, pos):
    best = None
    for cid, rs in reg.items():
        for c, s, e in rs:
            if c == chrom:
                d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
                if best is None or d < best[0]:
                    best = (d, cid)
    return best


for label, fset in ((f"FP only in {vb}", fb - fa), (f"FP only in {va}", fa - fb)):
    by_c = defaultdict(list)
    for chrom, pos, sv in sorted(fset):
        d, cid = container_of(chrom, pos)
        by_c[cid].append((chrom, pos, sv, d))
    print(f"== {label}: {len(fset)}")
    for cid, calls in sorted(by_c.items(), key=lambda kv: -len(kv[1])):
        print(
            f"  container {cid}: {len(calls)} FP ({calls[0][0]}:{calls[0][1]}, dist {calls[0][3]})"
            f"  {va}: {sa.get(cid, '-')}  |  {vb}: {sb.get(cid, '-')}"
        )
