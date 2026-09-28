"""Time the phasing all-vs-all of dumped containers in two minimap2 modes and
count the read pairs each one aligns.

  joint:  one index of all reads, every other read a secondary hit of the
          query (-N <reads> -p 0.05; the current run_ava)
  pertarget: one minimap2 call per target read i with the reads j > i as
          queries, -N 0: one best alignment per pair, no secondaries (a
          multi-part index, -I, would do the same in one call, but minimap2
          then aligns every read to itself and both directions of each pair)

usage: ava_modes.py <dump_dir> <out.tsv> <crID> [...]
(reads as dumped by dump_phasing_reads.py --orient; capped to the 50 longest)
"""

import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pysam
from Bio import SeqIO

PARAMS = "-k15 -w5 -e0 -m100 -r2k --eqx -t 1"
TIMEOUT = 600


def count(sams: list[Path]):
    records, pairs = 0, set()
    for sam in sams:
        with pysam.AlignmentFile(str(sam), "r", check_sq=False) as f:
            for a in f:
                if (
                    a.is_unmapped
                    or a.is_supplementary
                    or a.query_name == a.reference_name
                ):
                    continue
                records += 1
                pairs.add(tuple(sorted((a.query_name, a.reference_name))))
    return records, len(pairs)


def minimap(args: list[str], sam: Path) -> None:
    with open(sam, "w") as out:
        subprocess.run(
            args, stdout=out, stderr=subprocess.DEVNULL, check=True, timeout=TIMEOUT
        )


def joint(reads, tmp: Path):
    fa = tmp / "reads.fa"
    SeqIO.write(reads, fa, "fasta")
    sam = tmp / "joint.sam"
    cmd = f"minimap2 -a {PARAMS} -D --dual=no --secondary=yes -N {len(reads)} -p 0.05"
    minimap(cmd.split() + [str(fa), str(fa)], sam)
    return [sam]


def pertarget(reads, tmp: Path):
    reads = sorted(reads, key=lambda r: r.id)
    sams = []
    for i, target in enumerate(reads[:-1]):
        t_fa, q_fa = tmp / f"t{i}.fa", tmp / f"q{i}.fa"
        SeqIO.write([target], t_fa, "fasta")
        SeqIO.write(reads[i + 1 :], q_fa, "fasta")
        sam = tmp / f"pt{i}.sam"
        minimap(f"minimap2 -a {PARAMS} -N 0".split() + [str(t_fa), str(q_fa)], sam)
        sams.append(sam)
    return sams


def run(mode, reads, tmp: Path):
    t0 = time.perf_counter()
    try:
        sams = mode(reads, tmp)
    except subprocess.TimeoutExpired:
        return TIMEOUT, -1, -1
    dt = time.perf_counter() - t0
    return (dt, *count(sams))


dump, out_tsv = Path(sys.argv[1]), Path(sys.argv[2])
with open(out_tsv, "w") as out:
    out.write("crID\treads\tmode\tsecs\trecords\tpairs\n")
    for crID in sys.argv[3:]:
        reads = list(SeqIO.parse(dump / f"{crID}.fa", "fasta"))
        reads = sorted(reads, key=lambda r: (-len(r.seq), r.id))[:50]
        for name, mode in (("joint", joint), ("pertarget", pertarget)):
            with tempfile.TemporaryDirectory() as tmp:
                dt, records, pairs = run(mode, reads, Path(tmp))
            out.write(f"{crID}\t{len(reads)}\t{name}\t{dt:.2f}\t{records}\t{pairs}\n")
            out.flush()
