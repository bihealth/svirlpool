"""CPU profile of phase_reads' Python side on dumped read sets.

The all-vs-all alignments are computed once up front and replayed (run_ava is
patched), so the profile holds only the Python work, measured in process CPU
time (robust to other load on the machine).

usage: profile_phasing.py <reads_dir> <n_containers> [param=value ...]
"""

import cProfile
import pstats
import sys
import tempfile
import time
from pathlib import Path

from Bio import SeqIO

from svirlpool.localassembly import read_phasing as rp

reads_dir, n = Path(sys.argv[1]), int(sys.argv[2])
kw = {k: eval(v) for k, v in (x.split("=", 1) for x in sys.argv[3:])}
params = rp.PhasingParams(**kw)
files = sorted(reads_dir.glob("*.fa"), key=lambda f: int(f.stem))
step = max(1, len(files) // n)
files = files[::step][:n]

orig_run_ava = rp.run_ava
cache = {}
sets = []
for f in files:
    reads = {r.id: r for r in SeqIO.parse(f, "fasta")}
    sets.append((f.stem, reads))
    with tempfile.TemporaryDirectory() as tmp:
        cache[f.stem] = orig_run_ava(reads, Path(tmp), params, threads=1, timeout=60)

current = {}
rp.run_ava = lambda reads, tmp_dir, params, threads, timeout: cache[current["id"]]
prof = cProfile.Profile(time.process_time)
t0 = time.process_time()
prof.enable()
for cid, reads in sets:
    current["id"] = cid
    rp.phase_reads(reads, params=params)
prof.disable()
print(f"{len(sets)} containers, python cpu {time.process_time() - t0:.2f} s")
pstats.Stats(prof).sort_stats("tottime").print_stats(18)
