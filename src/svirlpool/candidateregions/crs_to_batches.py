"""
Group candidate-region containers into adjacent batches.

This module is run after :mod:`svirlpool.candidateregions.crs_to_containers_db`
in the Snakemake workflow. It reads the containers database, sorts the
containers by genomic coordinate and greedily packs them into batches
that a single ``consensus`` job can process while keeping a shared read
cache. Containers whose CRs span more than ``--max-batch-span`` bases or
multiple chromosomes are emitted as their own single-container batches
(this preserves the simple sliding-window cache semantics for the common
case while still allowing BND-connected wide containers to be processed
correctly).

Outputs:
- ``consensus_batches.tsv``: TSV with header
  ``batch_id\tchr\tstart\tend\tn_containers\tcrIDs``
  where ``crIDs`` is a comma-separated list of representative container
  crIDs (matching keys in the containers database).
- ``batch_ids.txt``: newline-separated batch ids, used by the Snakemake
  checkpoint expansion to fan out one consensus job per batch.
"""

from __future__ import annotations

import argparse
import json
import logging
import sqlite3
from pathlib import Path

from ..util import datatypes  # noqa: F401 -- imported for cattrs side-effects

log = logging.getLogger(__name__)


def _load_container_extents(
    path_db: Path,
) -> dict[int, tuple[list[str], int, int]]:
    """Return ``{rep_crID: (chrs_sorted, min_start, max_end)}`` for every container."""
    if not path_db.exists():
        raise FileNotFoundError(f"Containers database {path_db} does not exist.")
    extents: dict[int, tuple[list[str], int, int]] = {}
    conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
    try:
        c = conn.cursor()
        c.execute("SELECT crID, data FROM containers")
        for rep_crID, data in c.fetchall():
            payload = json.loads(data)
            crs = payload.get("crs", [])
            if not crs:
                continue
            chrs: list[str] = sorted({str(cr["chr"]) for cr in crs})
            min_start = min(int(cr["referenceStart"]) for cr in crs)
            max_end = max(int(cr["referenceEnd"]) for cr in crs)
            extents[int(rep_crID)] = (chrs, min_start, max_end)
    finally:
        conn.close()
    return extents


def build_batches(
    extents: dict[int, tuple[list[str], int, int]],
    batch_size: int,
    max_batch_span: int,
) -> list[dict]:
    """Greedy pack containers into batches.

    Wide containers (multi-chromosome or span > ``max_batch_span``)
    become single-container batches. The remaining containers are sorted
    by ``(chr, min_start, rep_crID)`` and packed in order with a hard
    cap of ``batch_size`` containers per batch and a hard cap of
    ``max_batch_span`` on the total spanned range within one batch.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be > 0")
    if max_batch_span <= 0:
        raise ValueError("max_batch_span must be > 0")

    wide: list[tuple[int, list[str], int, int]] = []
    local: list[tuple[str, int, int, int]] = []  # (chr, start, end, rep_crID)
    for rep_crID, (chrs, start, end) in extents.items():
        span = end - start
        if len(chrs) > 1 or span > max_batch_span:
            wide.append((rep_crID, chrs, start, end))
        else:
            local.append((chrs[0], start, end, rep_crID))

    local.sort()

    batches: list[dict] = []
    next_id = 0

    # Each wide container becomes its own batch.
    for rep_crID, chrs, start, end in wide:
        batches.append({
            "batch_id": next_id,
            "chr": ",".join(chrs),
            "start": start,
            "end": end,
            "crIDs": [rep_crID],
        })
        next_id += 1

    # Greedy pack local containers.
    current: list[tuple[str, int, int, int]] = []
    current_chr: str | None = None
    current_start: int | None = None
    current_end: int | None = None

    def _flush() -> None:
        nonlocal current, current_chr, current_start, current_end, next_id
        if not current:
            return
        batches.append({
            "batch_id": next_id,
            "chr": current_chr,
            "start": current_start,
            "end": current_end,
            "crIDs": [rep for _, _, _, rep in current],
        })
        next_id += 1
        current = []
        current_chr = None
        current_start = None
        current_end = None

    for chr_, start, end, rep_crID in local:
        if current_chr is None:
            current_chr = chr_
            current_start = start
            current_end = end
        new_start = start if current_start is None else min(current_start, start)
        new_end = end if current_end is None else max(current_end, end)
        span = new_end - new_start
        if chr_ != current_chr or len(current) >= batch_size or span > max_batch_span:
            _flush()
            current_chr = chr_
            current_start = start
            current_end = end
        else:
            current_start = new_start
            current_end = new_end
        current.append((chr_, start, end, rep_crID))
    _flush()

    return batches


def write_outputs(batches: list[dict], tsv_path: Path, ids_path: Path) -> None:
    with open(tsv_path, "w") as f:
        f.write("batch_id\tchr\tstart\tend\tn_containers\tcrIDs\n")
        for b in batches:
            crIDs_csv = ",".join(str(c) for c in b["crIDs"])
            f.write(
                f"{b['batch_id']}\t{b['chr']}\t{b['start']}\t{b['end']}\t{len(b['crIDs'])}\t{crIDs_csv}\n"
            )
    with open(ids_path, "w") as f:
        for b in batches:
            f.write(f"{b['batch_id']}\n")


def crs_containers_to_batches(
    path_containers: Path,
    tsv_out: Path,
    ids_out: Path,
    batch_size: int = 100,
    max_batch_span: int = 20_000_000,
) -> list[dict]:
    extents = _load_container_extents(path_containers)
    log.info(f"Loaded {len(extents)} containers from {path_containers}.")
    batches = build_batches(
        extents=extents, batch_size=batch_size, max_batch_span=max_batch_span
    )
    log.info(
        f"Built {len(batches)} batches (batch_size={batch_size}, max_batch_span={max_batch_span})."
    )
    write_outputs(batches=batches, tsv_path=tsv_out, ids_path=ids_out)
    log.info(f"Wrote {tsv_out} and {ids_out}.")
    return batches


def run(args, **kwargs):
    log_level = (
        getattr(logging, args.log_level) if hasattr(args, "log_level") else logging.INFO
    )
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        force=True,
    )
    crs_containers_to_batches(
        path_containers=args.input,
        tsv_out=args.tsv,
        ids_out=args.ids,
        batch_size=args.batch_size,
        max_batch_span=args.max_batch_span,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Groups candidate-region containers (output of crs_to_containers_db) "
            "into batches of genomically adjacent containers for the consensus stage."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the crs containers SQLite database.",
    )
    parser.add_argument(
        "-t",
        "--tsv",
        type=Path,
        required=True,
        help="Output TSV path (one row per batch).",
    )
    parser.add_argument(
        "--ids",
        type=Path,
        required=True,
        help="Output text file with one batch id per line (Snakemake checkpoint input).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Maximum number of containers per batch (default: 100).",
    )
    parser.add_argument(
        "--max-batch-span",
        type=int,
        default=20_000_000,
        help="Maximum total reference span (bp) per batch (default: 20Mb).",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
