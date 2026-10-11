"""Can the core-alignment merge reach the labelled records at all?

For each tr_substring seed record (base_c8804b1): are the child's composite(s)
(its HG002 consensus IDs, same type, <= 1 kb) and a composite of the
transmitting parent (t_parent) in one candidate group, and did the alignment
judge/link them? Uses the dumped composites and cores; aligns per group.

usage: reachable.py [klass ...]
"""

import csv
import gzip
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

ROOT = Path("/home/mayv_c/development/svp_mendel")
P = ROOT / "proto"


def base_records():
    out = {}
    with gzip.open(ROOT / "tr_substring/data/family.mndl.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t")
            info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
            out[(t[0], int(t[1]), info["SVTYPE"], abs(int(info["SVLEN"])))] = info["CONSENSUSIDs"].split(",")
    return out


def main():
    want = set(sys.argv[1:]) or {"A1", "B", "A2", "C", "D_cons", "D_gt"}
    labels = [r for r in csv.DictReader(open(ROOT / "tr_substring/out/record_classes.tsv"), delimiter="\t")
              if r["cls"] == "a" and r["klass"] in want]
    base = base_records()
    comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
    cores = pickle.load(open(P / "cores.pkl", "rb"))
    by_chr = defaultdict(list)
    for s in comps:
        if not issubclass(s.sv_type, SVpatterns.SVpatternAdjacency):
            by_chr[s.ref_start[0]].append(s)
    params = cam.CoreAlnParams()
    tally = defaultdict(Counter)
    rows = []
    cache = {}
    for lab in labels:
        key = (lab["chrom"], int(lab["pos"]), lab["svtype"], int(lab["svlen"]))
        cons = base.get(key)
        k = lab["klass"]
        tally[k]["n"] += 1
        if cons is None:
            tally[k]["no_base"] += 1
            continue
        child = {tuple(c.split(":")) for c in cons if c.startswith("HG002:")}
        tp = lab["t_parent"]
        chrom = lab["chrom"]
        if chrom not in cache:
            cs = by_chr[chrom]
            groups = M._candidate_groups_for_chromosome(chrom, cs, 150)
            gid = {}
            for gi, cc in enumerate(groups):
                for i in cc:
                    gid[id(cs[i])] = gi
            cache[chrom] = (cs, groups, gid)
        cs, groups, gid = cache[chrom]
        cx = [s for s in cs if cam.composite_consensus_key(s) in child and s.sv_type.get_sv_type() == lab["svtype"]
              and abs(s.ref_start[1] - key[1]) <= 1000]
        if not cx:
            tally[k]["child_composite_missing"] += 1
            continue
        g = {gid[id(s)] for s in cx}
        members = [cs[i] for gi in g for i in groups[gi]]
        parent = [s for s in members if s.svPatterns[0].samplename == tp and s.sv_type.get_sv_type() == lab["svtype"]]
        if not parent:
            tally[k]["no_tparent_composite_in_group"] += 1
            continue
        tally[k]["tparent_in_group"] += 1
        indels = [s for s in members if issubclass(s.sv_type, (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion))]
        placed, oriented = cam.place_composites(indels, cores)
        spans, drift = cam.core_reference_spans(indels, placed)
        keys = {p.key for p in placed if p is not None}
        st = cam.GroupStats()
        alns = cam.align_group_cores(oriented, keys, params, st, {kk: v[0] for kk, v in spans.items()}, drift)
        pos = {id(s): i for i, s in enumerate(indels)}
        verdicts = Counter()
        for x in cx:
            for y in parent:
                px, py = placed[pos[id(x)]], placed[pos[id(y)]]
                if px is None or py is None:
                    verdicts["unplaced"] += 1
                    continue
                d = cam.decide_pair(px, py, alns.get((px.key, py.key)), len(oriented[px.key]), len(oriented[py.key]),
                                    params, spans.get(px.key), spans.get(py.key), alns.get((py.key, px.key)))
                v = "link" if d.link else "distinct" if d.distinct else "twin" if d.both_same else (
                    "judged_unlinked" if d.evaluable else "not_evaluable:" + str(d.detail.get("reason")))
                verdicts[v] += 1
        best = next((v for v in ("link", "twin", "distinct", "judged_unlinked") if verdicts[v]), None) or (
            verdicts.most_common(1)[0][0] if verdicts else "none")
        tally[k]["best:" + best] += 1
        rows.append((k, lab["idx"], lab["chrom"], lab["pos"], lab["svtype"], lab["svlen"], tp, best, dict(verdicts)))
    for k in sorted(tally):
        print(k, dict(sorted(tally[k].items())))
    with open(P / "reachable.tsv", "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
