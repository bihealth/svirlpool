"""Features for routing a container to the cheap reference-site phasing (arm B)
or the all-vs-all phasing (arm A), computable in production (no truth).

Per container, on the reads of ref_vs_ava_eval.py, arm B is rerun step by step
(as ref_read_phasing.phase_reads_reference) to record:

  mapq_mean, frac_mapq_lt20/_lt60   alignments of the phased reads in the window
  frac_multi_aln                    reads with > 1 alignment in the window
  div_med, div_iqr                  read divergence to the reference
  cand_cols, site_cols, kept_cols   SNV columns: candidates (top base >= 3 and
                                    second >= 2), with a site, after recurrence
  noise_ratio                       1 - kept_cols / site_cols
  third_base                        mean share of reads with a third base at
                                    kept columns
  af_dev                            mean |minor share - 0.5| at kept columns
  discord                           share of (read, kept column) observations
                                    that disagree with the read's B group
                                    majority (B's partition self-consistency)
  B_* (status, n_alleles, min group share, assigned share)

usage: ref_route_features.py <workdir> <out.tsv> [--procs 8] [--limit N]
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from kmeans_gate import BUFFER_CLIPPED, MAX_CN, PHASING_FLANK

from svirlpool.localassembly import consensus, read_phasing, ref_read_phasing as rrp
from svirlpool.localassembly import read_cache as read_cache_mod
from svirlpool.signalprocessing import copynumber_tracks

logging.disable(logging.WARNING)


def features(preads, alns, ref, chrom, start, end, params=None):
    p = params or read_phasing.PhasingParams()
    out = {}
    if len(preads) < 2 * p.min_group:
        out["B_status"] = "no_information"
        return out
    reads = read_phasing.select_reads(preads, p)
    w0 = max(0, start - PHASING_FLANK)
    w1 = min(ref.get_reference_length(chrom), end + PHASING_FLANK)
    refseq = ref.fetch(chrom, w0, w1).upper()
    ref_code = rrp._CODE[rrp._as_array(refseq)]
    by_read = defaultdict(list)
    mapq, seen = [], set()
    for a in alns:
        if (a.query_name not in reads or a.reference_name != chrom or a.is_secondary
                or a.is_unmapped or a.query_sequence is None or a.cigartuples is None
                or a.reference_end <= w0 or a.reference_start >= w1):
            continue
        key = (a.query_name, a.reference_start, a.reference_end)
        if key in seen:
            continue
        seen.add(key)
        mapq.append(a.mapping_quality)
        by_read[a.query_name].append(rrp._parse(a, ref_code, w0, w1))
    mapq = np.array(mapq) if mapq else np.zeros(1)
    out["mapq_mean"] = float(mapq.mean())
    out["frac_mapq_lt20"] = float((mapq < 20).mean())
    out["frac_mapq_lt60"] = float((mapq < 60).mean())
    out["frac_multi_aln"] = float(np.mean([len(v) > 1 for v in by_read.values()])) if by_read else 0.0
    div = []
    for al in by_read.values():
        e = sum(a.n_x + a.n_small for a in al)
        n = sum(a.n_eq for a in al) + e
        if n:
            div.append(e / n)
    div = np.array(div) if div else np.zeros(1)
    out["div_med"] = float(np.median(div))
    out["div_iqr"] = float(np.subtract(*np.percentile(div, [75, 25])))

    lowq = rrp.low_quality_reads_ref(by_read, p)
    names = sorted(r for r in reads if r not in lowq)
    # column matrix as in snv_sites_ref
    n_col = len(refseq)
    C = np.full((len(names), n_col), -1, dtype=np.int8)
    for i, r in enumerate(names):
        s = np.zeros(n_col, dtype=bool)
        for a in by_read.get(r, []):
            c = a.col_pos - w0
            ok = (a.col_code < 4) & ~rrp._blocked(a, a.col_pos, p.indel_window)
            tw = s[c]
            C[i, c[ok & ~tw]] = a.col_code[ok & ~tw]
            C[i, c[tw]] = -1
            s[c] = True
    counts = np.stack([(C == b).sum(axis=0) for b in range(4)], axis=1)
    srt = np.sort(counts, axis=1)
    hp = rrp._homopolymer_mask(refseq, p.homopolymer_len)
    cand = (srt[:, -1] >= p.min_alt) & (srt[:, -2] >= p.min_ref) & ~hp
    out["cand_cols"] = int(cand.sum())
    s_snv = rrp.snv_sites_ref(by_read, names, refseq, w0, p)
    site_cols = {s.pos for s in s_snv}
    kept = read_phasing.recurrent_sites(s_snv, min_support=p.recurrence)
    kept_cols = sorted({s.pos for s in kept})
    out["site_cols"] = len(site_cols)
    out["kept_cols"] = len(kept_cols)
    out["noise_ratio"] = 1 - len(kept_cols) / len(site_cols) if site_cols else np.nan
    s_sv = rrp.sv_sites_ref(by_read, names, w0, w1, p)
    res = read_phasing.phase_from_sites(preads, reads, kept, s_sv, lowq, p)
    out["B_status"] = res.status
    out["B_n_alleles"] = res.n_alleles
    out["B_n_sv"] = len(s_sv)
    g = res.groups
    cnt = np.bincount(list(g.values())) if g else np.zeros(1)
    out["B_min_frac"] = float(cnt.min() / cnt.sum()) if len(cnt) > 1 else 1.0
    out["B_assigned"] = len(g) / max(1, len(preads))
    if kept_cols:
        kc = np.array(kept_cols) - w0
        sub = counts[kc]
        o = np.sort(sub, axis=1)
        tot = sub.sum(axis=1)
        out["third_base"] = float(np.mean(o[:, -3] / np.maximum(tot, 1)))
        out["af_dev"] = float(np.mean(np.abs(o[:, -2] / np.maximum(o[:, -1] + o[:, -2], 1) - 0.5)))
        # B self-consistency
        idx = {r: i for i, r in enumerate(names)}
        dis = tot_obs = 0
        for lab in set(g.values()):
            rows = [idx[r] for r, lb in g.items() if lb == lab and r in idx]
            if not rows:
                continue
            M = C[np.ix_(rows, kc)]
            for j in range(M.shape[1]):
                col = M[:, j]
                col = col[col >= 0]
                if len(col) < 2:
                    continue
                maj = np.bincount(col, minlength=4).max()
                dis += len(col) - maj
                tot_obs += len(col)
        out["discord"] = dis / tot_obs if tot_obs else np.nan
    return out


def one_batch(args):
    wd, crIDs, cfg = args
    containers = consensus.load_crs_containers_from_db(path_db=wd / "crs_containers.db", crIDs=crIDs)
    items = sorted(containers.items(), key=lambda it: (
        min(cr.chr for cr in it[1]["crs"]), min(cr.referenceStart for cr in it[1]["crs"]), it[0]))
    cache = read_cache_mod.ReadSequenceCache(path_alignments=Path(cfg["alignments"]))
    ref = pysam.FastaFile(cfg["reference"])
    rows = []
    try:
        for rep, container in items:
            crs = {cr.crID: cr for cr in container["crs"]}
            chrom = min(cr.chr for cr in crs.values())
            cn = max(2, copynumber_tracks.query_copynumber_from_regions(
                bgzip_bed=wd / "copy_number_tracks.bed.gz",
                regions=[(c.chr, c.referenceStart, c.referenceEnd) for c in crs.values()]))
            if cn > MAX_CN:
                cache.advance(window_start=float("inf"))
                continue
            alns, recs = {}, {}
            for cr in crs.values():
                a, s = cache.fetch_for_cr(cr)
                alns[cr.crID] = a
                recs.update(s)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns, buffer_clipped_length=BUFFER_CLIPPED)
            cutreads = consensus.trim_reads(
                dict_alignments=alns, intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns, buffer_clipped_length=BUFFER_CLIPPED,
                flank=PHASING_FLANK)
            preads = consensus.trim_reads(
                dict_alignments=alns, intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs)
            preads = {rn: r for rn, r in preads.items() if rn in cutreads}
            on_chr = [c for c in crs.values() if c.chr == chrom]
            t0 = time.perf_counter()
            row = {"crID": rep}
            row.update(features(preads, [a for al in alns.values() for a in al], ref, chrom,
                                min(c.referenceStart for c in on_chr),
                                max(c.referenceEnd for c in on_chr)))
            row["feat_secs"] = time.perf_counter() - t0
            rows.append(row)
            cache.advance(window_start=float("inf"))
    finally:
        cache.close()
        ref.close()
    return rows


def main():
    from multiprocessing import Pool

    p = argparse.ArgumentParser()
    p.add_argument("workdir", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--procs", type=int, default=8)
    p.add_argument("--limit", type=int, default=0)
    a = p.parse_args()
    cfg = json.load(open(a.workdir / "config.json"))
    jobs = []
    with open(a.workdir / "consensus_batches.tsv") as f:
        ic = f.readline().rstrip("\n").split("\t").index("crIDs")
        for line in f:
            crIDs = [int(x) for x in line.rstrip("\n").split("\t")[ic].split(",")]
            jobs.append(crIDs)
    if a.limit:
        jobs = jobs[: a.limit]
    jobs = [(a.workdir, c[i : i + 10], cfg) for c in jobs for i in range(0, len(c), 10)]
    with Pool(a.procs) as pool:
        rows = [r for rs in pool.imap_unordered(one_batch, jobs) for r in rs]
    pd.DataFrame(rows).sort_values("crID").to_csv(a.out, sep="\t", index=False)
    print(f"{len(rows)} containers -> {a.out}")


if __name__ == "__main__":
    main()
