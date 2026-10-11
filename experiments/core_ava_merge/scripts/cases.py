"""Alignment verdict for the study's inconsistent child records and their parent partners.

For each row of partners.tsv with a partner record (has=True): the child record's
composites (its consensus IDs, same type, within 1 kb of POS) and the partner
record's composites are aligned core vs core and judged by decide_pair.
Writes cases_verdicts.tsv; prints a tally.

usage: cases.py [FRAC_TOL]
"""

import csv
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")
S = Path("/home/mayv_c/development/svp_mendel/study")


def main():
    frac = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
    params = cam.CoreAlnParams(frac_tol=frac)
    records = {r["id"]: r for r in csv.DictReader(open(S / "records.tsv"), delimiter="\t")}
    partners = [r for r in csv.DictReader(open(S / "partners.tsv"), delimiter="\t") if r["has"] == "True"]
    comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
    cores = pickle.load(open(P / "cores.pkl", "rb"))
    by_cons = defaultdict(list)
    for svc in comps:
        key = cam.composite_consensus_key(svc)
        if key is not None:
            by_cons[key].append(svc)

    def composites_of(rec):
        typ = rec["svtype"]
        pos = int(rec["pos"])
        out = []
        for cid in rec["consensus_ids"].split(","):
            s, c = cid.split(":")
            for svc in by_cons.get((s, c), []):
                if svc.sv_type.get_sv_type() != typ or svc.ref_start[0] != rec["chrom"]:
                    continue
                if abs(svc.ref_start[1] - pos) <= 1000:
                    out.append(svc)
        return out

    rows = []
    tally = Counter()
    for pr in partners:
        child = records[pr["id"]]
        partner = records.get(pr["pid"])
        if partner is None:
            continue
        a = composites_of(child)
        b = composites_of(partner)
        # keep only the child's own (HG002) composites and the parents' ones
        a = [x for x in a if x.svPatterns[0].samplename == "HG002"]
        b = [y for y in b if y.svPatterns[0].samplename != "HG002"]
        if not a or not b:
            tally["composites_not_found"] += 1
            rows.append({"id": pr["id"], "pid": pr["pid"], "cls": pr["cls"], "verdict": "composites_not_found"})
            continue
        group = a + b
        placed, oriented = cam.place_composites(group, cores)
        keys = {p.key for p in placed if p is not None}
        stats = cam.GroupStats()
        alns = cam.align_group_cores(oriented, keys, params, stats) if len(keys) > 1 else {}
        best = None
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                px, py = placed[i], placed[len(a) + j]
                if px is None or py is None:
                    verdict, det = "unplaced", {}
                else:
                    aln = alns.get((px.key, py.key))
                    if aln is None:
                        verdict, det = "no_alignment", {}
                    else:
                        d = cam.decide_pair(px, py, aln, len(oriented[px.key]), len(oriented[py.key]), params)
                        det = d.detail
                        if d.link:
                            verdict = "link"
                        elif d.distinct:
                            verdict = "distinct"
                        elif d.evaluable:
                            verdict = "evaluable_no_link"
                        else:
                            verdict = "not_evaluable"
                rank = ["link", "distinct", "evaluable_no_link", "not_evaluable", "no_alignment", "unplaced"].index(verdict)
                old = (M.can_merge_svComposites_insertions if child["svtype"] == "INS" else M.can_merge_svComposites_deletions)(
                    x, y, apriori_size_difference_fraction_tolerance=0.1, near=150, min_kmer_overlap=0.7
                )
                cand = (rank, verdict, x.get_size(), y.get_size(), det, old, px, py)
                if best is None or cand[0] < best[0]:
                    best = cand
        rank, verdict, sx, sy, det, old, px, py = best
        tally[verdict] += 1
        rows.append({
            "id": pr["id"], "pid": pr["pid"], "cls": pr["cls"], "chrom": child["chrom"], "pos": child["pos"],
            "svtype": child["svtype"], "child_len": child["svlen"], "partner_len": partner["svlen"],
            "n_child_comps": len(a), "n_partner_comps": len(b), "size_x": sx, "size_y": sy,
            "verdict": verdict, "old_rule": old,
            "x_indel": det.get("x_indel"), "y_indel": det.get("y_indel"),
            "x_eval": det.get("x_eval"), "y_eval": det.get("y_eval"),
            "x_on_y": det.get("x_on_y"), "y_on_x": det.get("y_on_x"),
            "x_key": px.key if px else None, "y_key": py.key if py else None,
        })
    out = P / f"cases_verdicts_{frac}.tsv"
    cols = sorted({k for r in rows for k in r}, key=lambda k: list(rows[-1].keys()).index(k) if k in rows[-1] else 99)
    with open(out, "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(out)
    for k, v in tally.most_common():
        print(k, v)
    by_cls = Counter((r["cls"], r["verdict"]) for r in rows)
    for k, v in sorted(by_cls.items()):
        print(*k, v, sep="\t")


if __name__ == "__main__":
    main()
