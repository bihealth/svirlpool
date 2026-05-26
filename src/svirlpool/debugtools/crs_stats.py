#!/usr/bin/env python3
"""Load candidate regions and print summary statistics to the terminal."""

import argparse
import logging
from collections import Counter
from pathlib import Path

from ..candidateregions import signalstrength_to_crs

log = logging.getLogger(__name__)


def compute_percentiles(values: list[int], steps: int = 10) -> dict[int, int]:
    if not values:
        return {}

    sorted_values = sorted(values)
    n = len(sorted_values)
    percentiles = {}
    for p in range(0, 101, steps):
        index = min(max(int(round(p / 100 * (n - 1))), 0), n - 1)
        percentiles[p] = sorted_values[index]
    return percentiles


def crs_stats(crs: list[signalstrength_to_crs.datatypes.CandidateRegion]) -> None:
    region_sizes = [cr.referenceEnd - cr.referenceStart for cr in crs]
    chromosomes = [cr.chr for cr in crs]
    read_counts = [len(cr.get_read_names()) for cr in crs]

    log.info("Candidate regions loaded: %d", len(crs))

    print("=== CRS Summary ===")
    print(f"Total candidate regions: {len(crs)}")
    print("")

    print("Region size distribution (percentiles):")
    for p, value in compute_percentiles(region_sizes, steps=10).items():
        print(f"  {p:3d}th percentile: {value}")
    print("")

    print("Signal read counts per region:")
    print(f"  min: {min(read_counts) if read_counts else 0}")
    print(f"  max: {max(read_counts) if read_counts else 0}")
    print(f"  mean: {sum(read_counts) / len(read_counts) if read_counts else 0:.2f}")
    print(f"  median: {compute_percentiles(read_counts, steps=50).get(50, 0)}")
    print("")

    print("Regions per chromosome:")
    for chrom, count in Counter(chromosomes).most_common():
        print(f"  {chrom}: {count}")
    print("")

    print("20 largest candidate regions:")
    largest = sorted(crs, key=lambda x: x.referenceEnd - x.referenceStart, reverse=True)[:20]
    for cr in largest:
        size = cr.referenceEnd - cr.referenceStart
        print(
            f"  {cr.chr}\t{cr.referenceStart}\t{cr.referenceEnd}\t{size}\t{cr.crID}"
        )


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Print statistics for a candidate regions database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the candidate regions database file (e.g. crs.db).",
        metavar="DB",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser


def main() -> int:
    parser = get_parser()
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if not args.input.exists():
        log.error("Input database file not found: %s", args.input)
        return 1

    try:
        crs = signalstrength_to_crs.load_crs_from_db(path_db=args.input)
    except Exception as exc:
        log.error("Failed to load candidate regions: %s", exc)
        return 1

    crs_stats(crs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())