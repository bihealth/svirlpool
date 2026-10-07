"""Step 3 of the read-selection study: the consensus of a container built from
the top-k crossing reads of each CR only (locus_reads.py: the reads crossing
the CR in one alignment or collinear split alignments, ranked by how far they
run beyond it on their shorter side).

Runs svirlpool's batch driver (crs_containers_to_consensus, the pipeline's
escalation and settings) on the studied containers with the read cache
filtered to the selected reads. A CR without crossing reads gets no reads.
Writes <out>/consensus/0/consensus.batch_<i>.{jsonl,log}, the layout
container_value.py and consensus_q100.py read.

usage: readsel_consensus.py <locus dir> <out dir> --k 40 [--groups slow,mid]
         [--workers 6] [--threads 4] [--work <svirlpool workdir of the run>]
"""

from __future__ import annotations

import argparse
import json
import logging
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from svirlpool.localassembly import consensus as C
from svirlpool.localassembly import read_cache as RC

WORK = "/home/mayv_c/development/svp_tiered15/results/cr_base/20x/HG002/work"
FMT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


class SelectedReads:
    """ReadSequenceCache that returns only the selected reads of each CR."""

    keep: dict[int, set[str]] = {}

    def __init__(self, **kw):
        self.inner = RC.ReadSequenceCache(**kw)

    def fetch_for_cr(self, cr):
        alns, seqs = self.inner.fetch_for_cr(cr)
        k = self.keep.get(cr.crID, set())
        return ([a for a in alns if a.query_name in k],
                {r: s for r, s in seqs.items() if r in k})

    def __getattr__(self, name):
        return getattr(self.inner, name)


def run(job):
    i, cids, keep, out, work, threads = job
    cfg = json.load(open(Path(work) / "config.json"))
    d = Path(out) / "consensus" / "0"
    d.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    h = logging.FileHandler(d / f"consensus.batch_{i}.log", mode="w")
    h.setFormatter(logging.Formatter(FMT))
    root.addHandler(h)
    root.setLevel(logging.INFO)
    C.log.setLevel(logging.INFO)
    SelectedReads.keep = keep
    C.read_cache_mod.ReadSequenceCache = SelectedReads
    C.crs_containers_to_consensus(
        samplename=cfg["samplename"], input=Path(work) / "crs_containers.db",
        copy_number_tracks=Path(work) / "copy_number_tracks.bed.gz",
        output=d / f"consensus.batch_{i}.jsonl", lamassemble_mat=cfg["lamassemble_mat"],
        path_alignments=Path(cfg["alignments"]), threads=threads,
        buffer_clipped_sequence=500, consensus_method=cfg["consensus_method"],
        reference=Path(cfg["reference"]),
        escalation=C.parse_escalation(cfg["consensus_escalation"]), crIDs=cids,
        max_padding_size=cfg["max_padding_size"],
        max_copy_number_threshold=cfg["max_consensus_copy_number"],
        clustering_mode=cfg["consensus_clustering_mode"], phasing_flank=cfg["phasing_flank"],
        phasing_fallback=cfg["phasing_fallback"],
        clustering_strategy=cfg["clustering_strategy"],
    )
    return i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("locus", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--groups", default="slow,mid")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", default=WORK)
    a = ap.parse_args()
    cont = pd.read_csv(a.locus / "containers.tsv", sep="\t")
    cont = cont[cont.group.isin(a.groups.split(",")) & cont.db.str.startswith(a.work)]
    reads = pd.read_csv(a.locus / "reads.tsv", sep="\t")
    reads = reads[reads.db.str.startswith(a.work) & (reads.anchor >= 0)]
    reads = reads.sort_values("anchor", ascending=False).groupby("crID").head(a.k)
    keep = {int(cr): set(g.read) for cr, g in reads.groupby("crID")}
    # longest-first over the workers, by the pipeline's seconds
    cids = cont.sort_values("seconds", ascending=False).container.astype(int).tolist()
    chunks = [cids[w::a.workers] for w in range(a.workers)]
    jobs = [(w, ch, keep, str(a.out), a.work, a.threads) for w, ch in enumerate(chunks) if ch]
    a.out.mkdir(parents=True, exist_ok=True)
    with Pool(len(jobs), maxtasksperchild=1) as pool:
        for i in pool.imap_unordered(run, jobs):
            print(f"worker {i} done", flush=True)


if __name__ == "__main__":
    main()
