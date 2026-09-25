"""Exact check of the coordinate mirroring used by parse_pair_both: every '='
column must map q's forward base (complemented for reverse alignments) onto
t's base, and every inverted mismatch must hold a base that differs from q's.

usage: check_inversion_exact.py <reads.fa> [...]
"""

import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

from svirlpool.localassembly import read_phasing as rp

RC = str.maketrans("ACGTN", "TGCAN")
for path in sys.argv[1:]:
    reads = {r.id: r for r in SeqIO.parse(path, "fasta")}
    seqs = {n: str(r.seq) for n, r in reads.items()}
    with tempfile.TemporaryDirectory() as tmp:
        alns = rp.run_ava(reads, Path(tmp), rp.PhasingParams(), threads=1, timeout=600)
    n_eq = bad_eq = n_x = bad_x = n_rev = 0
    for a in alns:
        q, t = seqs[a.query_name], seqs[a.reference_name]
        lq, rev = len(q), a.is_reverse
        n_rev += rev
        ct = a.cigartuples
        rpos, qpos = a.reference_start, (ct[0][1] if ct[0][0] == 5 else 0)
        for op, ln in ct:
            if op == 7:
                for k in range(ln):
                    fq = lq - 1 - (qpos + k) if rev else qpos + k
                    qb = q[fq].translate(RC) if rev else q[fq]
                    n_eq += 1
                    bad_eq += qb != t[rpos + k]
            if op in (7, 8):
                rpos += ln
                qpos += ln
            elif op == 2:
                rpos += ln
            elif op in (1, 4):
                qpos += ln
        _, inv = rp.parse_pair_both(a, seqs)
        for pos, b in inv.mism.items():
            n_x += 1
            bad_x += q[pos] == b
    print(
        f"{Path(path).name}: alignments {len(alns)} (reverse {n_rev})  "
        f"'=' columns {n_eq} mismatching under the mapping {bad_eq}  "
        f"inverted X {n_x} equal to q's base {bad_x}"
    )
