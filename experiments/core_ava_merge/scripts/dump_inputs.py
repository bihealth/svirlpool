"""Dump the pre-merge SVcomposites of the trio and the cores of their consensuses.

Outputs (proto dir): all_svComposites.pkl, cores.pkl ({(sample, cid): (core_start, seq)}).
"""

import logging
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from svirlpool.localassembly import svirltile
from svirlpool.svcalling import core_alignment_merge as cam
from svirlpool.svcalling.multisample_sv_calling import (
    SUPPORTED_SV_TYPES,
    generate_svComposites_from_dbs,
)

OUT = Path("/home/mayv_c/development/svp_mendel/proto")
DBS = [
    f"/home/mayv_c/development/svp_mendel/inputs/dbs/{s}/svirltile.db"
    for s in ("HG002", "HG003", "HG004")
]


def load(args):
    db, wanted = args
    t = time.time()
    cores, _ = cam.load_core_sequences(db, wanted)
    return db, cores, time.time() - t


def main():
    logging.basicConfig(level=logging.INFO)
    print("svirlpool from", cam.__file__, flush=True)
    t = time.time()
    if (OUT / "all_svComposites.pkl").exists():
        comps = pickle.load(open(OUT / "all_svComposites.pkl", "rb"))
    else:
        comps = generate_svComposites_from_dbs(
            input=DBS, sv_types=set(SUPPORTED_SV_TYPES), collapse_repeats="proximal", collapse_max_gap=50
        )
        with open(OUT / "all_svComposites.pkl", "wb") as f:
            pickle.dump(comps, f)
    print(f"{len(comps)} composites in {time.time() - t:.0f}s", flush=True)
    wanted: dict[str, set[str]] = {}
    for svc in comps:
        for p in svc.svPatterns:
            wanted.setdefault(p.samplename, set()).add(p.consensusID)
    names = {db: svirltile.get_metadata(Path(db))["samplename"] for db in DBS}
    cores = {}
    with ProcessPoolExecutor(3) as ex:
        for db, c, dt in ex.map(load, [(db, wanted[names[db]]) for db in DBS]):
            print(f"{names[db]}: {len(c)} cores, {sum(len(s) for _, s in c.values()) / 1e6:.1f} Mb in {dt:.0f}s", flush=True)
            for cid, v in c.items():
                cores[(names[db], cid)] = v
    with open(OUT / "cores.pkl", "wb") as f:
        pickle.dump(cores, f)
    print("done", flush=True)


if __name__ == "__main__":
    sys.exit(main())
