"""The truth columns the allele-level evaluation needs, for any container DB,
without a full run: per container the T2TQ100 / V5 category (ref / het / hom /
cpx_het over the locus span), the benchmark-region fraction of the locus and
the TRF flag. Locus span and categories as container_truth.py.

usage: container_truth_light.py <crs_containers.db> <out.tsv>
"""

import json
import sqlite3
import sys

import pysam
from container_truth_lib import TRF, TRUTH_BENCH, TRUTH_VCF, Intervals, load_bed, locus_span_fn, truth_at

DB, OUT = sys.argv[1], sys.argv[2]
trf = Intervals(load_bed(TRF, gap=500))
bench = {k: Intervals(load_bed(b)) for k, b in TRUTH_BENCH.items()}
vfs = {k: pysam.VariantFile(v) for k, v in TRUTH_VCF.items()}
locus_span = locus_span_fn(trf)
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
with open(OUT, "w") as f:
    cols = ["crID", "chr", "start", "end", "trf"]
    for k in TRUTH_VCF:
        cols += [f"{k}_cat", f"{k}_bench_frac"]
    f.write("\t".join(cols) + "\n")
    for crid, data in con.execute("SELECT crID, data FROM containers ORDER BY crID"):
        crs = json.loads(data)["crs"]
        c = crs[0]["chr"]
        s = min(cr["referenceStart"] for cr in crs)
        e = max(cr["referenceEnd"] for cr in crs)
        lo, hi = locus_span(c, s, e)
        row = [crid, c, s, e, bool(trf.overlapping(c, s, e))]
        for k in TRUTH_VCF:
            row += [truth_at(vfs[k], c, lo, hi)["cat"], round(bench[k].covered(c, lo, hi) / (hi - lo), 3)]
        f.write("\t".join(map(str, row)) + "\n")
print(f"-> {OUT}")
