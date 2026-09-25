"""Dump the read set that ``consensus_while_phasing`` phases, per container.

Cuts the reads exactly as the consensus stage does (CR +- CR_CUT_FLANK for the
consensus, CR +- phasing_flank for phasing, restricted to the consensus reads)
and writes ``<out_dir>/<crID>.fa`` for every container, so phasing variants
can be run and timed without the BAM (``phasing_eval.py``).

usage: dump_phasing_reads.py <workdir> <out_dir> [--flank 10000] [--procs 12]
"""

from __future__ import annotations

import argparse
import json
import logging
from multiprocessing import Pool
from pathlib import Path

from Bio import SeqIO

from svirlpool.localassembly import consensus
from svirlpool.localassembly import read_cache as read_cache_mod

logging.disable(logging.INFO)


def dump_batch(args):
    wd, out_dir, crIDs, flank, cfg, orient = args
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
    n = 0
    try:
        for rep, container in items:
            crs = {cr.crID: cr for cr in container["crs"]}
            alns, recs = {}, {}
            for cr in crs.values():
                a, s = cache.fetch_for_cr(cr)
                alns[cr.crID] = a
                recs.update(s)
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()), dict_alignments=alns, buffer_clipped_length=500
            )
            cutreads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            iv = consensus.get_read_alignment_intervals_in_cr(
                crs=list(crs.values()),
                dict_alignments=alns,
                buffer_clipped_length=500,
                flank=flank,
            )
            reads = consensus.trim_reads(
                dict_alignments=alns,
                intervals=consensus.get_max_extents_of_read_alignments_on_cr(iv),
                read_records=recs,
            )
            reads = {rn: r for rn, r in reads.items() if rn in cutreads}
            if orient:
                reads = consensus.orient_reads_to_reference(reads, alns)
            with open(out_dir / f"{rep}.fa", "w") as f:
                SeqIO.write([reads[rn] for rn in sorted(reads)], f, "fasta")
            n += 1
            cache.advance(window_start=float("inf"))
    finally:
        cache.close()
    return n


def main():
    p = argparse.ArgumentParser()
    p.add_argument("workdir", type=Path)
    p.add_argument("out_dir", type=Path)
    p.add_argument("--flank", type=int, default=10000)
    p.add_argument("--procs", type=int, default=12)
    p.add_argument(
        "--orient", action="store_true", help="reads on the reference strand"
    )
    a = p.parse_args()
    cfg = json.load(open(a.workdir / "config.json"))
    a.out_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    with open(a.workdir / "consensus_batches.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        ic = header.index("crIDs")
        for line in f:
            crIDs = [int(x) for x in line.rstrip("\n").split("\t")[ic].split(",")]
            jobs.append((a.workdir, a.out_dir, crIDs, a.flank, cfg, a.orient))
    with Pool(a.procs) as pool:
        print("containers dumped:", sum(pool.imap_unordered(dump_batch, jobs)))


if __name__ == "__main__":
    main()
