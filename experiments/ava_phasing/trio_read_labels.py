#!/usr/bin/env python
"""Trio-SNV ground truth: assign HG002 ONT reads to the paternal / maternal haplotype.

Independent of any SV signal: only SNVs where the child (HG002) is heterozygous and
the parents (HG003 = father, HG004 = mother) make the parental origin of the two child
alleles unambiguous are used.

Pipeline
  1. windows: per container, span = min(referenceStart)..max(referenceEnd) of its CRs,
     window = span +- 30 kb; overlapping windows are merged per chromosome.
  2. pileup: `bcftools mpileup -X ont -q20 -Q10 -a AD,DP` on child+father+mother jointly
     (SNVs only), parallel over <=2 Mb chunks of the merged windows. Only sites where
     the child has >= 2 reads of a non-REF allele are kept (prefilter on raw AD).
     No `bcftools call`: parental genotype classes are defined on raw allele fractions,
     so AD is all we need and we avoid the (ONT-uncalibrated) genotype priors.
  3. sites: child het (AF 0.25-0.75, DP>=8), parents DP>=6 classified hom-ref
     (AF<=0.05) / het (0.25-0.75) / hom-alt (AF>=0.95); informative iff exactly one
     parent is homozygous or both are homozygous for different alleles; both-hom-same
     = Mendelian inconsistent (dropped). Further filters (counted): multi-allelic
     evidence (third allele > 10% of child reads), homopolymer >= 4 at/adjacent to the
     site (with either allele), and extreme depth (any sample DP > 2.5x its median).
  4. reads: every HG002 alignment (primary + supplementary, MAPQ>=1, no secondary /
     dup / qcfail) in a merged window; the base at each site (base quality >= 10) votes
     paternal or maternal; other bases are ignored. Votes of all alignments of a read
     name are pooled over the container's own +-30 kb window (one vote per site).
  5. labels: pat if n_pat>=2 and n_pat>=0.8*(n_pat+n_mat); mat symmetric; else unknown.
     Only read names with an alignment overlapping container span +- 1 kb are output.

Outputs (in data/):
  trio_read_labels.tsv.gz         crID read_name n_pat n_mat label
  trio_informative_sites.tsv.gz   chr pos(1-based) ref alt pat_allele mat_allele + AD
  trio_container_summary.tsv.gz   per-container QC numbers
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sqlite3
import shlex
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam

HERE = Path(__file__).resolve().parent
ALN = "/home/mayv_c/biodata/local/alignments"
DEFAULTS = dict(
    child=f"{ALN}/HG002.minimap2.20x.softclipped.bam",
    father=f"{ALN}/HG003.minimap2.20x.softclipped.bam",
    mother=f"{ALN}/HG004.minimap2.20x.softclipped.bam",
    ref="/home/mayv_c/biodata/local/references/GRCh38/GRCh38.fa",
    db="/home/mayv_c/development/svp_improvements/results/flank200/20x/HG002/work/crs_containers.db",
    outdir=str(HERE / "data"),
    tmp=os.path.join(os.environ.get("CLAUDE_JOB_DIR", "/home/mayv_c/.claude/jobs/bd7213cd"), "tmp", "trio"),
)

FLANK = 30_000
READ_FLANK = 1_000
CHUNK = 2_000_000


# ----------------------------------------------------------------------------- windows
def load_containers(db):
    con = sqlite3.connect(db)
    out = []
    for crid, data in con.execute("SELECT crID, data FROM containers"):
        crs = json.loads(data)["crs"]
        chrs = {c["chr"] for c in crs}
        assert len(chrs) == 1, (crid, chrs)
        out.append(
            (
                int(crid),
                crs[0]["chr"],
                min(c["referenceStart"] for c in crs),
                max(c["referenceEnd"] for c in crs),
            )
        )
    return sorted(out, key=lambda x: (x[1], x[2]))


def merge_windows(containers, chrom_len):
    by_chr = defaultdict(list)
    for _, chrom, s, e in containers:
        by_chr[chrom].append((max(0, s - FLANK), min(chrom_len[chrom], e + FLANK)))
    merged = []
    for chrom, iv in by_chr.items():
        iv.sort()
        cs, ce = iv[0]
        for s, e in iv[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                merged.append((chrom, cs, ce))
                cs, ce = s, e
        merged.append((chrom, cs, ce))
    return sorted(merged)


# ----------------------------------------------------------------------------- pileup
def run_pileup(job):
    chrom, s, e, args, out = job
    if os.path.exists(out):
        return out
    region = f"{chrom}:{s + 1}-{e}"
    mp = [
        "bcftools", "mpileup", "-X", "ont", "-I", "-q", "20", "-Q", "10", "-d", "1000",
        "-a", "FORMAT/AD,FORMAT/DP", "-f", args.ref, "-r", region, "-Ou",
        args.child, args.father, args.mother,
    ]
    # child (sample 0) must carry >=2 reads of the first ALT allele. NB: the filter has to be
    # applied with `view` (site level); `query -i` with FORMAT fields drops failing samples.
    cmd = (
        "set -o pipefail; " + " ".join(shlex.quote(x) for x in mp)
        + " 2>/dev/null | bcftools view -i 'FMT/AD[0:1]>=2' -Ou"
        + " | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'"
    )
    tmp = out + ".part"
    with open(tmp, "w") as fh:
        subprocess.run(["bash", "-c", cmd], stdout=fh, check=True)
    os.replace(tmp, out)
    return out


def hp_run(fa, chrom, pos0, allele):
    """Longest homopolymer touching pos0 when pos0 carries `allele` (0-based pos)."""
    s = max(0, pos0 - 12)
    seq = list(fa.fetch(chrom, s, pos0 + 13).upper())
    i = pos0 - s
    seq[i] = allele
    best = 0
    for j in (i - 1, i, i + 1):  # runs containing the site or its direct neighbours
        if j < 0 or j >= len(seq):
            continue
        b = seq[j]
        a = j
        while a - 1 >= 0 and seq[a - 1] == b:
            a -= 1
        z = j
        while z + 1 < len(seq) and seq[z + 1] == b:
            z += 1
        best = max(best, z - a + 1)
    return best


def classify_parent(af):
    if af <= 0.05:
        return "hr"
    if af >= 0.95:
        return "ha"
    if 0.25 <= af <= 0.75:
        return "het"
    return None


def select_sites(pileup_files, args, stats):
    fa = pysam.FastaFile(args.ref)
    rows = []
    for f in pileup_files:
        with open(f) as fh:
            for line in fh:
                chrom, pos, ref, alt, *ads = line.rstrip("\n").split("\t")
                alts = alt.split(",")
                ads = [[int(x) if x != "." else 0 for x in a.split(",")] for a in ads]
                # pick the child's most supported non-<*> alt allele
                cands = [k for k in range(1, len(alts) + 1) if alts[k - 1] != "<*>" and len(alts[k - 1]) == 1]
                if not cands or len(ref) != 1:
                    continue
                k = max(cands, key=lambda k: ads[0][k])
                rows.append((chrom, int(pos), ref.upper(), alts[k - 1].upper(), [(a[0], a[k], sum(a)) for a in ads]))
    stats["candidate_sites"] = len(rows)
    dp = np.array([[r[4][i][2] for i in range(3)] for r in rows])
    med = np.median(dp, axis=0)
    stats["median_dp_child_father_mother_at_candidates"] = med.tolist()
    max_dp = 2.5 * med

    sites = []
    c = defaultdict(int)
    for chrom, pos, ref, alt, ad in rows:
        (cr, ca, cd), (fr, fa_, fd), (mr, ma, md) = ad
        if cd < 8 or not (0.25 <= ca / cd <= 0.75):
            continue
        c["child_het"] += 1
        if fd < 6 or md < 6:
            c["parent_low_dp"] += 1
            continue
        if cd - cr - ca > 0.1 * cd:
            c["multiallelic"] += 1
            continue
        gf, gm = classify_parent(fa_ / fd), classify_parent(ma / md)
        if gf is None or gm is None:
            c["parent_ambiguous_af"] += 1
            continue
        if gf == gm and gf != "het":
            c["mendel_inconsistent"] += 1
            continue
        if gf == "het" and gm == "het":
            c["both_het_uninformative"] += 1
            continue
        c["informative"] += 1
        # parent homozygous for X transmits X; the child's other allele came from the other parent
        if gf == "hr" or gm == "ha":
            pat, mat = ref, alt
        else:  # gf == "ha" or gm == "hr"
            pat, mat = alt, ref
        if cd > max_dp[0] or fd > max_dp[1] or md > max_dp[2]:
            c["high_depth"] += 1
            continue
        hp = max(hp_run(fa, chrom, pos - 1, ref), hp_run(fa, chrom, pos - 1, alt))
        if hp >= 4:
            c["homopolymer_ge4"] += 1
            if not args.keep_homopolymer:
                continue
        sites.append((chrom, pos, ref, alt, pat, mat, gf, gm, ad, hp))
    stats["site_filter_counts"] = dict(c)
    stats["sites_used"] = len(sites)
    return sites


# ----------------------------------------------------------------------------- reads
def aln_query_pos(read, site_pos0):
    """Query positions for 0-based reference positions (sorted np array); -1 if not aligned (M/=/X)."""
    rs, qs, ln = [], [], []
    r, q = read.reference_start, 0
    for op, l in read.cigartuples:
        if op in (0, 7, 8):
            rs.append(r)
            qs.append(q)
            ln.append(l)
            r += l
            q += l
        elif op in (1, 4):
            q += l
        elif op in (2, 3):
            r += l
    rs, qs, ln = np.array(rs), np.array(qs), np.array(ln)
    b = np.searchsorted(rs, site_pos0, side="right") - 1
    ok = b >= 0
    bb = np.where(ok, b, 0)
    off = site_pos0 - rs[bb]
    ok &= off < ln[bb]
    return np.where(ok, qs[bb] + off, -1)


def vote_window(job):
    chrom, ws, we, site_pos0, pat_b, mat_b, bam_path, min_bq = job
    out = []  # (read_name, ref_start, ref_end, site_idx array, vote array (+1 pat / -1 mat))
    if len(site_pos0) == 0:
        site_pos0 = np.zeros(0, dtype=np.int64)
    bam = pysam.AlignmentFile(bam_path)
    for read in bam.fetch(chrom, ws, we):
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail or read.mapping_quality < 1:
            continue
        rs, re_ = read.reference_start, read.reference_end
        lo, hi = np.searchsorted(site_pos0, rs), np.searchsorted(site_pos0, re_)
        idx = np.arange(lo, hi)
        votes = np.zeros(0, dtype=np.int8)
        if len(idx):
            qp = aln_query_pos(read, site_pos0[idx])
            seq = read.query_sequence
            quals = read.query_qualities
            v = np.zeros(len(idx), dtype=np.int8)
            for j, p in enumerate(qp):
                if p < 0 or (quals is not None and quals[p] < min_bq):
                    continue
                base = seq[p]
                if base == pat_b[idx[j]]:
                    v[j] = 1
                elif base == mat_b[idx[j]]:
                    v[j] = -1
            keep = v != 0
            idx, votes = idx[keep], v[keep]
        out.append((read.query_name, rs, re_, idx.astype(np.int32), votes))
    return out


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k}", default=v)
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--min-bq", type=int, default=10)
    ap.add_argument("--keep-homopolymer", action="store_true", help="do not drop sites in/next to homopolymers >=4")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.tmp, exist_ok=True)
    stats = {}

    fa = pysam.FastaFile(args.ref)
    chrom_len = dict(zip(fa.references, fa.lengths))
    containers = load_containers(args.db)
    merged = merge_windows(containers, chrom_len)
    stats["containers"] = len(containers)
    stats["merged_windows"] = len(merged)
    stats["merged_window_bp"] = int(sum(e - s for _, s, e in merged))
    print(json.dumps(stats), file=sys.stderr)

    # --- pileup
    jobs = []
    for chrom, s, e in merged:
        for cs in range(s, e, CHUNK):
            ce = min(e, cs + CHUNK)
            jobs.append((chrom, cs, ce, args, os.path.join(args.tmp, f"pileup_{chrom}_{cs}_{ce}.tsv")))
    jobs.sort(key=lambda j: -(j[2] - j[1]))
    with ProcessPoolExecutor(args.threads) as ex:
        files = list(ex.map(run_pileup, jobs))
    print("pileup done", file=sys.stderr)

    # --- sites
    sites = select_sites(files, args, stats)
    print(json.dumps(stats), file=sys.stderr)
    with gzip.open(os.path.join(args.outdir, "trio_informative_sites.tsv.gz"), "wt") as fh:
        fh.write("chr\tpos\tref\talt\tpat_allele\tmat_allele\tfather_gt\tmother_gt\thp_run\t"
                 "child_ref\tchild_alt\tchild_dp\tfather_ref\tfather_alt\tfather_dp\tmother_ref\tmother_alt\tmother_dp\n")
        for chrom, pos, ref, alt, pat, mat, gf, gm, ad, hp in sites:
            fh.write("\t".join(map(str, [chrom, pos, ref, alt, pat, mat, gf, gm, hp, *[x for a in ad for x in a]])) + "\n")

    by_chr = defaultdict(list)
    for s in sites:
        by_chr[s[0]].append(s)
    vjobs, wsite = [], []
    for chrom, ws, we in merged:
        ss = [s for s in by_chr.get(chrom, []) if ws <= s[1] - 1 < we]
        pos0 = np.array([s[1] - 1 for s in ss], dtype=np.int64)
        vjobs.append((chrom, ws, we, pos0, [s[4] for s in ss], [s[5] for s in ss], args.child, args.min_bq))
        wsite.append(pos0)

    # --- votes
    with ProcessPoolExecutor(args.threads) as ex:
        wres = list(ex.map(vote_window, vjobs, chunksize=1))
    print("votes done", file=sys.stderr)

    # --- per-container labels
    win_of = {}
    for i, (chrom, ws, we) in enumerate(merged):
        win_of.setdefault(chrom, []).append((ws, we, i))
    lab_fh = gzip.open(os.path.join(args.outdir, "trio_read_labels.tsv.gz"), "wt")
    lab_fh.write("crID\tread_name\tn_pat\tn_mat\tlabel\n")
    sum_fh = gzip.open(os.path.join(args.outdir, "trio_container_summary.tsv.gz"), "wt")
    sum_fh.write("crID\tchr\tstart\tend\tn_sites\twindow_bp\tn_reads\tn_pat\tn_mat\tn_unknown\tn_novote\n")
    for crid, chrom, s, e in sorted(containers):
        wi = next(i for ws, we, i in win_of[chrom] if ws <= s and e <= we)
        pos0 = wsite[wi]
        cws, cwe = max(0, s - FLANK), e + FLANK
        slo, shi = np.searchsorted(pos0, cws), np.searchsorted(pos0, cwe)
        names = {r[0] for r in wres[wi] if r[1] < e + READ_FLANK and r[2] > s - READ_FLANK}
        per = defaultdict(dict)
        for name, rs, re_, idx, votes in wres[wi]:
            if name not in names:
                continue
            d = per[name]
            for i, v in zip(idx.tolist(), votes.tolist()):
                if slo <= i < shi and i not in d:
                    d[i] = v
        cnt = defaultdict(int)
        for name in sorted(names):
            vs = list(per[name].values())
            npat, nmat = vs.count(1), vs.count(-1)
            tot = npat + nmat
            if npat >= 2 and npat >= 0.8 * tot:
                lab = "pat"
            elif nmat >= 2 and nmat >= 0.8 * tot:
                lab = "mat"
            else:
                lab = "unknown"
            cnt[lab] += 1
            cnt["novote"] += tot == 0
            lab_fh.write(f"{crid}\t{name}\t{npat}\t{nmat}\t{lab}\n")
        sum_fh.write(f"{crid}\t{chrom}\t{s}\t{e}\t{shi - slo}\t{cwe - cws}\t{len(names)}\t"
                     f"{cnt['pat']}\t{cnt['mat']}\t{cnt['unknown']}\t{cnt['novote']}\n")
    lab_fh.close()
    sum_fh.close()
    with open(os.path.join(args.outdir, "trio_run_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    print(json.dumps(stats, indent=1), file=sys.stderr)


if __name__ == "__main__":
    main()
