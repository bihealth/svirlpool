#!/usr/bin/env python
"""Orthogonal check of the trio read labels against the phased GIAB v5.0q SV truth.

The labels themselves use no SV signal; here we only *evaluate* them: for isolated
heterozygous, phased DEL/INS (>=50 bp) in the GIAB Q100 VCF that lie within a container
span +-1 kb, every labelled read spanning the SV is typed as carrier / non-carrier from
its CIGAR, and we tabulate label x (carrier of hap1 / hap2). If hap1 == paternal and the
labels are clean, pat reads carry exactly the 1|0 SVs and mat reads the 0|1 SVs.
"""

from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import pysam

D = Path(__file__).resolve().parent / "data"
BAM = "/home/mayv_c/biodata/local/alignments/HG002.minimap2.20x.softclipped.bam"
VCF = "/home/mayv_c/biodata/local/benchmarks/HG002_GRCh38_v5.0q_stvar.filtered.vcf.gz"
PAD = 200
ISOLATE = 1000

summ = pd.read_csv(D / "trio_container_summary.tsv.gz", sep="\t")
lab = pd.read_csv(D / "trio_read_labels.tsv.gz", sep="\t")
lab = lab[lab.label != "unknown"]
labels = {cr: dict(zip(g.read_name, g.label)) for cr, g in lab.groupby("crID")}

vcf = pysam.VariantFile(VCF)
bam = pysam.AlignmentFile(BAM)
tab = Counter()
tabk = defaultdict(Counter)
per_sv = []
for cr, chrom, s, e in summ[["crID", "chr", "start", "end"]].itertuples(index=False):
    if cr not in labels:
        continue
    recs = list(vcf.fetch(chrom, max(0, s - 1000 - ISOLATE), e + 1000 + ISOLATE))
    for r in recs:
        if not (s - 1000 <= r.pos <= e + 1000):
            continue
        gt = r.samples[0]["GT"]
        if gt not in ((0, 1), (1, 0)) or not r.samples[0].phased:
            continue
        dl = len(r.ref) - len(r.alts[0])
        if abs(dl) < 50:
            continue
        kind = "DEL" if dl > 0 else "INS"
        L = abs(dl)
        rs, re_ = r.pos - 1, r.pos - 1 + len(r.ref)
        # isolated: no other truth variant >=10bp within ISOLATE bp
        if any(o is not r and abs(len(o.ref) - len(o.alts[0])) >= 10 and o.pos - 1 < re_ + ISOLATE and o.pos - 1 + len(o.ref) > rs - ISOLATE for o in recs):
            continue
        hap = 1 if gt == (1, 0) else 2
        c = Counter()
        seen = set()
        for a in bam.fetch(chrom, rs - PAD, re_ + PAD):
            if a.is_secondary or a.is_unmapped or a.mapping_quality < 1 or a.query_name in seen:
                continue
            lb = labels[cr].get(a.query_name)
            if lb is None or a.reference_start > rs - PAD or a.reference_end < re_ + PAD:
                continue
            seen.add(a.query_name)
            ref = a.reference_start
            ind = 0
            for op, l in a.cigartuples:
                if op in (0, 7, 8):
                    ref += l
                elif op == 2:
                    if kind == "DEL" and l >= 10 and rs - PAD <= ref <= re_ + PAD:
                        ind += l
                    ref += l
                elif op == 1 and kind == "INS" and l >= 10 and rs - PAD <= ref <= re_ + PAD:
                    ind += l
            call = "carrier" if ind >= 0.6 * L else ("noncarrier" if ind <= 0.3 * L else "ambig")
            if call == "ambig":
                continue
            carries_hap = hap if call == "carrier" else 3 - hap
            tab[(lb, f"hap{carries_hap}")] += 1
            tabk[kind][(lb, carries_hap)] += 1
            c[(lb, carries_hap)] += 1
        per_sv.append((cr, chrom, r.pos, kind, L, f"{gt[0]}|{gt[1]}", dict(c)))

print("label x haplotype implied by SV carrier status (read-SV pairs):")
for k in sorted(tab):
    print(" ", k, tab[k])
n_pat1, n_pat2 = tab[("pat", "hap1")], tab[("pat", "hap2")]
n_mat1, n_mat2 = tab[("mat", "hap1")], tab[("mat", "hap2")]
tot = n_pat1 + n_pat2 + n_mat1 + n_mat2
print(f"SVs tested: {len(per_sv)}; concordance assuming hap1=pat: {(n_pat1 + n_mat2) / max(tot, 1):.4f}")
for k, t in sorted(tabk.items()):
    ok = t[("pat", 1)] + t[("mat", 2)]
    print(f"  {k}: SVs {sum(p[3] == k for p in per_sv)}, read-SV pairs {sum(t.values())}, concordance {ok / max(sum(t.values()), 1):.4f}")
bad = [p for p in per_sv if (p[6].get(("pat", 2), 0) + p[6].get(("mat", 1), 0)) >= 2]
print(f"SVs with >=2 discordant reads (hap1=pat): {len(bad)}")
for p in bad[:30]:
    print("  ", p)
