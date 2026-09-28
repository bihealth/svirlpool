"""Step through ``read_phasing.phase_reads`` on one container at added error rates.

Per rate: SNV / SV sites (and how many split the reads along the trio
haplotypes), mean pair weight within / between haplotypes, the fraction of
positive between-haplotype and negative within-haplotype weights, and the
cluster sizes after correlation clustering and after refinement.

usage: phase_steps.py <reads_dir> <crID> [rates ...]
"""

import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO

from svirlpool.localassembly import read_phasing as rp
from svirlpool.localassembly.consensus import add_read_errors

sys.path.insert(0, str(Path(__file__).parent))
import eval_rates  # noqa: E402

reads_dir, crID = Path(sys.argv[1]), int(sys.argv[2])
rates = [float(x) for x in sys.argv[3:]] or [0.0, 0.04, 0.06, 0.08]
trio = eval_rates.load_trio().get(crID, {})
orig = {r.id: r for r in SeqIO.parse(reads_dir / f"{crID}.fa", "fasta")}
params = rp.PhasingParams()


def site_truth(s):
    a = {trio[r] for r in s.agree if r in trio}
    d = {trio[r] for r in s.disagree if r in trio}
    if not a or not d:
        return "?"
    return "hap" if len(a) == 1 and len(d) == 1 and a != d else "mixed"


for rate in rates:
    reads = {n: add_read_errors(r, 0, rate) for n, r in orig.items()} if rate else orig
    keep = sorted(reads, key=lambda n: (-len(reads[n].seq), n))[: params.max_reads]
    reads = {n: reads[n] for n in sorted(keep)}
    seqs = {n: str(r.seq) for n, r in reads.items()}
    with tempfile.TemporaryDirectory() as tmp:
        alns = rp.run_ava(reads, Path(tmp), params, threads=4, timeout=120)
        arrs = rp.ReadArrays(seqs)
        pairs = [p for a in alns for p in rp.parse_pair_both(a, arrs)]
    lowq = rp.low_quality_reads(pairs, params)
    pairs = [p for p in pairs if p.q not in lowq and p.t not in lowq]
    by_t = defaultdict(list)
    for p in pairs:
        by_t[p.t].append(p)
    raw = rp.snv_sites(by_t, seqs, params)
    snv = rp.recurrent_sites(raw, min_support=params.recurrence)
    sv = rp.sv_sites(by_t, seqs, params)
    names = sorted(r for r in reads if r not in lowq)
    W = rp.pair_weights(snv + sv, names, params.error_rate)
    hap = np.array([trio.get(n, "?") for n in names])
    lab = hap != "?"
    same = (hap[:, None] == hap[None, :]) & lab[:, None] & lab[None, :]
    diff = (hap[:, None] != hap[None, :]) & lab[:, None] & lab[None, :]
    np.fill_diagonal(same, False)
    cc = rp.correlation_cluster(W)
    ref = rp.refine_clusters(cc, names, snv + sv, W, params)
    aln_len = np.mean([p.t_end - p.t_start for p in pairs]) if pairs else 0
    print(
        f"== rate {rate}: reads {len(names)} (lowq {len(lowq)}), pairs {len(pairs)}, mean aligned {aln_len:.0f} bp"
    )
    print(f"   snv raw {len(raw)} {dict(Counter(map(site_truth, raw)))}")
    print(f"   snv recurrent {len(snv)} {dict(Counter(map(site_truth, snv)))}")
    print(f"   sv {len(sv)} {dict(Counter(map(site_truth, sv)))}")
    print(
        f"   W within hap: mean {W[same].mean():.1f}, <0 {np.mean(W[same] < 0):.2f}, ==0 {np.mean(W[same] == 0):.2f} | "
        f"between: mean {W[diff].mean():.1f}, >0 {np.mean(W[diff] > 0):.2f}"
    )
    print(
        f"   clusters after correlation {sorted(Counter(cc).values(), reverse=True)[:8]}  after refine {dict(Counter(ref))}"
    )
