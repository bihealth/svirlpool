"""D5 probe: can the same core alignments genotype a parent that has no pattern in
a child-carrying record?

For every record of the core-alignment merge (align mode) that contains an HG002
INS/DEL composite (size >= 50) and lacks a pattern of parent P: look at P's
consensuses in the same candidate group (they are aligned anyway) and ask whether
the child's allele is present on them (alignment across the child SV's placement
range +- 150 bp has <= max(10, 0.1 * size) indel bases), absent (evaluable, more
differences), or not evaluable. A parent without any consensus in the group is
'no_consensus' (its consensuses there carry no INS/DEL pattern, or there are none).

usage: d5_probe.py OUTDIR [chr ...]
"""

import json
import pickle
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")
NEAR, KMER, TOL = 150, 0.7, 0.1
PARENTS = ("HG003", "HG004")


def presence(aln, x, len_x, params):
    px = cam.project(aln, x.amb_start, x.amb_end)
    if px is None:
        return "not_projected"
    if x.amb_start - aln.q0 < params.min_flank or aln.q1 - x.amb_end < params.min_flank:
        return "no_flank"
    w0 = max(aln.q0, x.amb_start - params.flank)
    w1 = min(aln.q1, x.amb_end + params.flank, len_x)
    indel, _ = aln.window_differences(w0, w1)
    return "present" if indel <= max(params.abs_tol, params.frac_tol * x.size) else "absent"


def job(args):
    chrom, cs, cores = args
    params = cam.CoreAlnParams()
    old_ins = lambda a, b: M.can_merge_svComposites_insertions(a, b, TOL, NEAR, KMER)  # noqa: E731
    old_del = lambda a, b: M.can_merge_svComposites_deletions(a, b, TOL, NEAR, KMER)  # noqa: E731
    tally = Counter()
    for cc in M._candidate_groups_for_chromosome(chrom, cs, NEAR):
        group = [cs[i] for i in cc]
        indels = [s for s in group if issubclass(s.sv_type, (SVpatterns.SVpatternInsertion, SVpatterns.SVpatternDeletion))]
        if not any(s.svPatterns[0].samplename == "HG002" for s in indels):
            continue
        placed, oriented = cam.place_composites(indels, cores)
        spans, drift = cam.core_reference_spans(indels, placed)
        keys = {p.key for p in placed if p is not None}
        st = cam.GroupStats()
        alns = cam.align_group_cores(oriented, keys, params, st, {k: v[0] for k, v in spans.items()}, drift) if len(keys) > 1 else {}
        ins = [s for s in indels if issubclass(s.sv_type, SVpatterns.SVpatternInsertion)]
        dels = [s for s in indels if issubclass(s.sv_type, SVpatterns.SVpatternDeletion)]
        pos = {id(s): i for i, s in enumerate(indels)}
        records = []
        for sub, rule in ((ins, old_ins), (dels, old_del)):
            if not sub:
                continue
            records += cam.merge_same_type_with_core_alignments(
                sub, [placed[pos[id(s)]] for s in sub], oriented, alns, rule, "align", params, st, None, spans
            )
        pat2comp = {id(p): s for s in indels for p in s.svPatterns}
        parent_keys = {pp: {k for k in keys if k[0] == pp} for pp in PARENTS}
        for rec in records:
            if rec.get_size() < 50:
                continue
            comps = {id(pat2comp[id(p)]): pat2comp[id(p)] for p in rec.svPatterns}.values()
            child = [s for s in comps if s.svPatterns[0].samplename == "HG002" and placed[pos[id(s)]] is not None]
            if not child:
                continue
            x = max(child, key=lambda s: s.get_size())
            px = placed[pos[id(x)]]
            samples = {p.samplename for p in rec.svPatterns}
            for pp in PARENTS:
                if pp in samples:
                    tally[(pp, "has_pattern")] += 1
                    continue
                verdicts = set()
                for k in parent_keys[pp]:
                    aln = alns.get((px.key, k))
                    verdicts.add(presence(aln, px, len(oriented[px.key]), params) if aln else "no_alignment")
                if not parent_keys[pp]:
                    v = "no_consensus_in_group"
                elif "present" in verdicts:
                    v = "present"
                elif "absent" in verdicts:
                    v = "absent"
                else:
                    v = "not_evaluable"
                tally[(pp, v)] += 1
    return dict(tally)


def main():
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    chroms = sys.argv[2:]
    comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
    cores = pickle.load(open(P / "cores.pkl", "rb"))
    by_chr = defaultdict(list)
    for s in comps:
        if not issubclass(s.sv_type, SVpatterns.SVpatternAdjacency):
            by_chr[s.ref_start[0]].append(s)
    if chroms:
        by_chr = {c: by_chr[c] for c in chroms}
    jobs = []
    for c, cs in sorted(by_chr.items(), key=lambda kv: -len(kv[1])):
        keys = {(p.samplename, p.consensusID) for s in cs for p in s.svPatterns}
        jobs.append((c, cs, {k: cores[k] for k in keys if k in cores}))
    total = Counter()
    with ProcessPoolExecutor(8) as ex:
        for t in ex.map(job, jobs):
            total.update(t)
    rows = {f"{k[0]}:{k[1]}": v for k, v in sorted(total.items())}
    (outdir / "d5_probe.json").write_text(json.dumps(rows, indent=1))
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
