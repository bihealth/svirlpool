"""Old vs core-alignment merge at given loci, with the pair decisions.

usage: loci.py chr:pos[:label] ...
"""

import pickle
import sys
from pathlib import Path

from svirlpool.localassembly import SVpatterns
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling import svcomposite_merging as M

P = Path("/home/mayv_c/development/svp_mendel/proto")


def short(svc):
    ids = ",".join(sorted({f"{p.samplename}:{p.consensusID}" for p in svc.svPatterns}))
    sizes = "+".join(str(p.get_size()) for p in svc.svPatterns)
    return f"{svc.sv_type.get_sv_type()} {svc.ref_start[1]} size={svc.get_size()} ({sizes}) [{ids}]"


def main():
    comps = pickle.load(open(P / "all_svComposites.pkl", "rb"))
    cores = pickle.load(open(P / "cores.pkl", "rb"))
    params = cam.CoreAlnParams()
    for arg in sys.argv[1:]:
        chrom, pos = arg.split(":")[:2]
        pos = int(pos)
        near = [s for s in comps if s.ref_start[0] == chrom and abs(s.ref_start[1] - pos) < 5000
                and not issubclass(s.sv_type, SVpatterns.SVpatternAdjacency)]
        print(f"\n######## {arg}: {len(near)} composites within 5 kb")
        for cc in M._candidate_groups_for_chromosome(chrom, near, 150):
            group = [near[i] for i in cc]
            if not any(abs(s.ref_start[1] - pos) < 1500 for s in group) or len(group) < 2:
                continue
            old = M.vertically_merged_svComposites_from_group(group, 150, 0.7, 0.1)
            log = []
            st = cam.GroupStats()
            new, rest = cam.merge_group_with_core_alignments(
                group, cores,
                lambda a, b: M.can_merge_svComposites_insertions(a, b, 0.1, 150, 0.7),
                lambda a, b: M.can_merge_svComposites_deletions(a, b, 0.1, 150, 0.7),
                "align", params, st, log,
            )
            print("-- old records:")
            for s in sorted(old, key=lambda s: s.ref_start):
                if s.get_size() >= 30 and abs(s.ref_start[1] - pos) < 1500:
                    print("   ", short(s))
            print("-- new records:")
            for s in sorted(new + rest, key=lambda s: s.ref_start):
                if s.get_size() >= 30 and abs(s.ref_start[1] - pos) < 1500:
                    print("   ", short(s))
            print("-- pair decisions (cross-sample, |pos - locus| < 1500, sizes >= 30):")
            for r in log:
                pa = int(r["a"].split(":")[1].split("[")[0])
                pb = int(r["b"].split(":")[1].split("[")[0])
                if min(r["size_a"], r["size_b"]) < 30 or abs(pa - pos) > 1500 or abs(pb - pos) > 1500:
                    continue
                print(f"   {r['a']} {r['size_a']} vs {r['b']} {r['size_b']}: old={r['old']} aln={r['link_aln']} "
                      f"distinct={r['distinct']} twin={r.get('twin')} final={r['final']} corr={r.get('corr')} indel={r.get('indel')} "
                      f"same={r.get('same')} win={r.get('window')} reason={r.get('reason')}")


if __name__ == "__main__":
    main()
