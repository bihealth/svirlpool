"""Ablation: read phasing with sites from all-vs-all alignments (arm A,
``read_phasing.phase_reads``) vs. from the reads' reference alignments (arm B,
``ref_read_phasing.phase_reads_reference``). Everything else is shared.

Per container, on the reads the consensus stage uses (as kmeans_gate.py):
status, allele count, site counts, low-quality reads and time of each arm; the
raw read groups (read level: trio pair accuracy of the assigned reads) and the
partition production would assemble (groups when phased, else one group of
all reads but the low-quality ones; allele level vs T2TQ100, as
../ava_phasing/allele_eval.py).

usage: ref_vs_ava_eval.py <workdir> <out.tsv> [--procs 8] [--timeout 300] [--limit N]
                          [--truth container_truth.tsv] [--trio trio_read_labels.tsv.gz]
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import pandas as pd
import pysam
from kmeans_gate import (
    BUFFER_CLIPPED,
    MAX_CN,
    PHASING_FLANK,
    TRUTH,
    allele_eval,
)
from phasing_eval import load_trio, pair_accuracy

from svirlpool.localassembly import consensus, read_phasing, ref_read_phasing
from svirlpool.localassembly import read_cache as read_cache_mod
from svirlpool.signalprocessing import copynumber_tracks

logging.disable(logging.WARNING)


def _production(res: read_phasing.PhasingResult, cutreads) -> dict[str, int]:
    if res.status == "phased":
        return dict(res.groups)
    lowq = set(res.low_quality)
    return {rn: 0 for rn in cutreads if rn not in lowq}


def one_batch(args):
    wd, crIDs, cfg, timeout = args
    containers = consensus.load_crs_containers_from_db(
        path_db=wd / "crs_containers.db", crIDs=crIDs
    )
    items = sorted(
        containers.items(),
        key=lambda it: (
            min(cr.chr for cr in it[1]["crs"]),
            min(cr.referenceStart for cr in it[1]["crs"]),
            it[0],
        ),
    )
    cache = read_cache_mod.ReadSequenceCache(path_alignments=Path(cfg["alignments"]))
    ref = pysam.FastaFile(cfg["reference"])
    rows = []
    try:
        for rep, container in items:
            crs = {cr.crID: cr for cr in container["crs"]}
            chrom = min(cr.chr for cr in crs.values())
            row = {"crID": rep, "chr": chrom,
                   "start": min(cr.referenceStart for cr in crs.values())}
            cn = max(2, copynumber_tracks.query_copynumber_from_regions(
                bgzip_bed=wd / "copy_number_tracks.bed.gz",
                regions=[(c.chr, c.referenceStart, c.referenceEnd) for c in crs.values()],
            ))
            row["skipped"] = cn > MAX_CN
            if row["skipped"]:
                rows.append(row)
                cache.advance(window_start=float("inf"))
                continue
            alns, recs = {}, {}
            for cr in crs.values():
                a, s = cache.fetch_for_cr(cr)
                alns[cr.crID] = a
                recs.update(s)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns,
                buffer_clipped_length=BUFFER_CLIPPED,
            )
            cutreads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            signals = [s for cr in crs.values() for s in cr.sv_signals]
            row["sig_trf_frac"] = (
                sum(s.repeatID != -1 for s in signals) / len(signals) if signals else 0.0
            )
            row["n_reads"] = len(cutreads)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns,
                buffer_clipped_length=BUFFER_CLIPPED, flank=PHASING_FLANK,
            )
            preads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            preads = {rn: r for rn, r in preads.items() if rn in cutreads}

            # arm A: as consensus_while_phasing (reads on the reference strand)
            t0 = time.perf_counter()
            res_a = read_phasing.phase_reads(
                reads=consensus.orient_reads_to_reference(preads, alns),
                threads=1, timeout=timeout,
            )
            row["A_secs"] = time.perf_counter() - t0
            # arm B: same read set, sites from the reference alignments
            on_chr = [c for c in crs.values() if c.chr == chrom]
            t0 = time.perf_counter()
            res_b = ref_read_phasing.phase_reads_reference(
                reads=preads,
                alns=[a for al in alns.values() for a in al],
                ref=ref, chrom=chrom,
                start=min(c.referenceStart for c in on_chr),
                end=max(c.referenceEnd for c in on_chr),
                flank=PHASING_FLANK,
            )
            row["B_secs"] = time.perf_counter() - t0
            for arm, res in (("A", res_a), ("B", res_b)):
                row[f"{arm}_status"] = res.status
                row[f"{arm}_n_alleles"] = res.n_alleles
                row[f"{arm}_n_snv"] = res.n_snv_sites
                row[f"{arm}_n_sv"] = res.n_sv_sites
                row[f"{arm}_n_lowq"] = len(res.low_quality)
                row[f"{arm}_groups"] = json.dumps(res.groups, sort_keys=True)
                row[f"{arm}_prod"] = json.dumps(_production(res, cutreads), sort_keys=True)
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
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--limit", type=int, default=0, help="first N batches only")
    p.add_argument("--truth", type=Path, default=TRUTH)
    p.add_argument("--trio", type=Path, default=None)
    a = p.parse_args()
    cfg = json.load(open(a.workdir / "config.json"))
    jobs = []
    with open(a.workdir / "consensus_batches.tsv") as f:
        ic = f.readline().rstrip("\n").split("\t").index("crIDs")
        for line in f:
            crIDs = [int(x) for x in line.rstrip("\n").split("\t")[ic].split(",")]
            jobs.append((a.workdir, crIDs, cfg, a.timeout))
    if a.limit:
        jobs = jobs[: a.limit]
    jobs = [(wd, c[i : i + 10], cf, t) for wd, c, cf, t in jobs for i in range(0, len(c), 10)]
    t0 = time.perf_counter()
    with Pool(a.procs) as pool:
        rows = [r for rs in pool.imap_unordered(one_batch, jobs) for r in rs]
    print(f"wall {time.perf_counter() - t0:.0f} s, {len(rows)} containers")

    d = pd.DataFrame(rows).sort_values("crID")
    truth = pd.read_csv(a.truth, sep="\t").set_index("crID")
    trio = load_trio(a.trio) if a.trio else load_trio()
    ev = []
    for r in d.itertuples():
        out = {"crID": r.crID}
        t_reads = trio.get(r.crID, {})
        if r.crID in truth.index:
            t = truth.loc[r.crID]
            out["truth_locus_ok"] = t.chr == r.chr and abs(int(t.start) - int(r.start)) < 5000
            out["cat"] = t.T2T_cat
            out["bench"] = t.T2T_bench_frac >= 1
            out["trf"] = bool(t.trf)
        if not r.skipped:
            out["trio_total"] = len(t_reads)
            for arm in ("A", "B"):
                m, acc = pair_accuracy(json.loads(getattr(r, f"{arm}_groups")), t_reads)
                out[f"{arm}_trio_assigned"] = m
                out[f"{arm}_acc"] = acc
                if "cat" in out:
                    e = allele_eval(json.loads(getattr(r, f"{arm}_prod")), t_reads, out["cat"])
                    for k, v in (e or {}).items():
                        out[f"{arm}_{k}"] = v
        ev.append(out)
    d = d.merge(pd.DataFrame(ev), on="crID", how="left")
    d.to_csv(a.out, sep="\t", index=False)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
