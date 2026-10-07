"""The sample's median depth, from the candidate-region containers: per
container the highest median signal coverage (ExtendedSVsignal.coverage) of
its CRs, the median over all containers. Most containers are ordinary loci,
so this is the typical read depth at a candidate region. Used by the consensus
stage's --read-selection-factor.

usage: python -m svirlpool.candidateregions.container_depth -i crs_containers.db -o depth.txt
"""

from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


def median_container_depth(path_db: Path) -> float:
    depths = []
    conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
    try:
        for (data,) in conn.execute("SELECT data FROM containers"):
            per_cr = [
                np.median([s["coverage"] for s in cr["sv_signals"]])
                for cr in json.loads(data).get("crs", [])
                if cr["sv_signals"]
            ]
            if per_cr:
                depths.append(max(per_cr))
    finally:
        conn.close()
    return float(np.median(depths)) if depths else 0.0


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("-i", "--input", type=Path, required=True)
    p.add_argument("-o", "--output", type=Path, required=True)
    a = p.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    depth = median_container_depth(a.input)
    log.info(f"median container depth: {depth:g}")
    a.output.write_text(f"{depth:g}\n")


if __name__ == "__main__":
    main()
