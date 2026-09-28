"""Time the phasing all-vs-all (joint mode, as run_ava) of dumped containers
under minimap2 parameter variants.

usage: [AVA_VARIANTS="-U1,20 ..."] ava_params.py <dump_dir> <out.tsv> <crID> [...]
(AVA_VARIANTS: extra minimap2 options per variant, added to "base")
(reads as dumped by dump_phasing_reads.py --orient; capped to the 50 longest)
"""

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pysam
from Bio import SeqIO

TIMEOUT = 300
# -U a,b: the minimizer occurrence threshold is max(a, min(b, the -f
# threshold)); minimizers occurring more often in the index are not seeds
VARIANTS = {
    "base": "-k15 -w5 -e0 -m100 -r2k",
    "f001": "-k15 -w5 -e0 -m100 -r2k -f0.001",
    "U20,30": "-k15 -w5 -e0 -m100 -r2k -U20,30",
    "U15,20": "-k15 -w5 -e0 -m100 -r2k -U15,20",
    "U10,15": "-k15 -w5 -e0 -m100 -r2k -U10,15",
    "U5,10": "-k15 -w5 -e0 -m100 -r2k -U5,10",
}

if os.environ.get("AVA_VARIANTS"):
    VARIANTS = {
        v: VARIANTS["base"] + ("" if v == "base" else f" {v}")
        for v in os.environ["AVA_VARIANTS"].split()
    }


dump, out_tsv = Path(sys.argv[1]), Path(sys.argv[2])
with open(out_tsv, "w") as out:
    out.write("crID\treads\tvariant\tsecs\tpairs\n")
    for crID in sys.argv[3:]:
        reads = list(SeqIO.parse(dump / f"{crID}.fa", "fasta"))
        reads = sorted(reads, key=lambda r: (-len(r.seq), r.id))[:50]
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "reads.fa"
            SeqIO.write(reads, fa, "fasta")
            for name, params in VARIANTS.items():
                sam = Path(tmp) / f"{name}.sam"
                cmd = (
                    f"minimap2 -a {params} -D --dual=no --eqx -t 1 --secondary=yes "
                    f"-N {len(reads)} -p 0.05"
                ).split() + [str(fa), str(fa)]
                t0 = time.perf_counter()
                try:
                    with open(sam, "w") as f:
                        subprocess.run(
                            cmd, stdout=f, stderr=subprocess.DEVNULL, timeout=TIMEOUT
                        )
                except subprocess.TimeoutExpired:
                    out.write(f"{crID}\t{len(reads)}\t{name}\t{TIMEOUT}\t-1\n")
                    out.flush()
                    continue
                dt = time.perf_counter() - t0
                pairs = set()
                with pysam.AlignmentFile(str(sam), "r", check_sq=False) as f:
                    for a in f:
                        if a.is_unmapped or a.is_supplementary:
                            continue
                        pairs.add(tuple(sorted((a.query_name, a.reference_name))))
                out.write(f"{crID}\t{len(reads)}\t{name}\t{dt:.2f}\t{len(pairs)}\n")
                out.flush()
