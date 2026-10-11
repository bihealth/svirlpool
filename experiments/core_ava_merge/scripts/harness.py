"""Run the old and the core-alignment vertical merge on the dumped composites.

usage: harness.py OUTDIR MODE FRAC_TOL [chr ...]
Writes OUTDIR/summary.json, OUTDIR/pairs.<chr>.jsonl, OUTDIR/records_{old,new}.pkl
(list of records, each a sorted list of pattern keys) for INS/DEL composites.
"""

import json
import pickle
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import attrs

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")
NEAR, KMER, TOL = 150, 0.7, 0.1


def pkey(p):
    return (p.samplename, p.consensusID, p.get_sv_type(), p.read_start, p.read_end, p.get_size(), len(p.get_supporting_reads()))


def job(args):
    chr_name, comps, cores, mode, params, outdir = args
    t0 = time.perf_counter()
    old = M.merge_svComposites_for_chromosome(chr_name, comps, TOL, NEAR, KMER)
    t_old = time.perf_counter() - t0
    t0 = time.perf_counter()
    new, stats = M.merge_svComposites_for_chromosome_core_aln(
        chr_name, comps, TOL, NEAR, KMER, cores, mode, params, outdir
    )
    t_new = time.perf_counter() - t0

    def recs(ms):
        return [
            (svc.sv_type.get_sv_type(), sorted(pkey(p) for p in svc.svPatterns), svc.get_size())
            for svc in ms
            if issubclass(svc.sv_type, (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion))
        ]

    return chr_name, t_old, t_new, attrs.asdict(stats), recs(old), recs(new)


def main():
    outdir = Path(sys.argv[1])
    mode = sys.argv[2]
    frac = float(sys.argv[3])
    chroms = [c for c in sys.argv[4:] if not c.startswith("--")]
    outdir.mkdir(parents=True, exist_ok=True)
    params = cam.CoreAlnParams(frac_tol=frac)
    write_pairs = "--no-pairs" not in sys.argv
    comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
    cores = pickle.load(open(P / "cores.pkl", "rb"))
    by_chr = defaultdict(list)
    for svc in comps:
        if issubclass(svc.sv_type, SVpatterns.SVpatternAdjacency):
            continue
        by_chr[svc.ref_start[0]].append(svc)
    if chroms:
        by_chr = {c: by_chr[c] for c in chroms}
    jobs = []
    for c, cs in sorted(by_chr.items(), key=lambda kv: -len(kv[1])):
        keys = {(p.samplename, p.consensusID) for s in cs for p in s.svPatterns}
        jobs.append((c, cs, {k: cores[k] for k in keys if k in cores}, mode, params, str(outdir) if write_pairs else None))
    summary = {"t_old": 0.0, "t_new": 0.0, "stats": cam.GroupStats()}
    rec_old, rec_new = [], []
    with ProcessPoolExecutor(8) as ex:
        for c, t_old, t_new, st, ro, rn in ex.map(job, jobs):
            summary["t_old"] += t_old
            summary["t_new"] += t_new
            summary["stats"].add(cam.GroupStats(**st))
            rec_old += ro
            rec_new += rn
            print(c, f"old {t_old:.1f}s new {t_new:.1f}s", st["minimap2_calls"], flush=True)
    summary["stats"] = attrs.asdict(summary["stats"])
    summary["n_old"] = len(rec_old)
    summary["n_new"] = len(rec_new)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=1))
    pickle.dump((rec_old, rec_new), open(outdir / "records.pkl", "wb"))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
