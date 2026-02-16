#!/usr/bin/env python3
"""
This script reads the log output from multisample_sv_calling and queries it by a given region.
It then prints all relevant decisions that regard this region.

Input arguments:
- log-file path
- at least one of:
    - region (chr:start-end)
    - crIDs (comma-separated list of candidate region IDs)
- region-margin (int, optional, default: 0) - margin to extend the region by when querying the log file

Procedure:
Both a region and crIDs can be provided. Both are used to find relevant lines in the log file.
If a region is provided:
  - Parse all results from the log file and store them in a dictionary of intervaltrees (chr:IntervalTree)
  - The intervaltree has as the data point the index of the line in the log file
If crIDs are provided:
  - Parse the log file line by line and store all lines that match the given crIDs
Print the relevant lines to the console in the order they appear in the log file.
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

from intervaltree import IntervalTree


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Trace SV calling decisions from log files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "log_file",
        type=Path,
        help="Path to the log file from multisample_sv_calling",
    )
    parser.add_argument(
        "--region",
        type=str,
        help="Genomic region to query (format: chr:start-end, e.g., 1:1000000-2000000)",
    )
    parser.add_argument(
        "--crIDs",
        type=str,
        help="Comma-separated list of candidate region IDs to query (e.g., 1,5,10)",
    )
    parser.add_argument(
        "--region-margin",
        type=int,
        default=0,
        help="Margin in bp to extend the region by when querying (default: 0)",
    )

    args = parser.parse_args()

    if not args.region and not args.crIDs:
        parser.error("At least one of --region or --crIDs must be provided")

    return args


def parse_region(region_str: str) -> Tuple[str, int, int]:
    """
    Parse a region string in the format chr:start-end.
    Handles commas in numbers (e.g., 1:1,000,000-2,000,000).

    Returns:
        Tuple of (chromosome, start, end)
    """
    # Match pattern with optional commas in numbers
    match = re.match(r"^(.+?):[\d,]+\-[\d,]+$", region_str)
    if not match:
        raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")

    # Split on colon to get chromosome and coordinates
    parts = region_str.split(":")
    if len(parts) != 2:
        raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")

    chrom = parts[0]
    coords = parts[1].split("-")
    if len(coords) != 2:
        raise ValueError(f"Invalid region format: {region_str}. Expected chr:start-end")

    # Remove commas and convert to int
    start = int(coords[0].replace(",", ""))
    end = int(coords[1].replace(",", ""))

    return chrom, start, end


def extract_region_from_log_line(line: str) -> List[Tuple[str, int, int]]:
    """
    Extract genomic region(s) from a log line.

    Returns:
        List of tuples (chr, start, end) found in the line
    """
    regions = []

    # Pattern: region=chr:start-end or chr_start_end:start-end (for candidate region format)
    # Note: These patterns don't use commas, as log files typically use unformatted numbers
    region_patterns = [
        r"region=([^:\s]+):(\d+)-(\d+)",  # region=chr:start-end
        r"region=(\d+)_(\d+)_(\d+):(\d+)-(\d+)",  # region=chr_refstart_refend:start-end
        r"regions=([^:\s]+):(\d+)-(\d+)",  # regions=chr:start-end
        r"\b(\d+)_(\d+)_(\d+):(\d+)-(\d+)",  # chr_refstart_refend:start-end (without prefix)
    ]

    for pattern in region_patterns:
        matches = re.finditer(pattern, line)
        for match in matches:
            groups = match.groups()
            if len(groups) == 3:
                chrom, start, end = groups
                regions.append((chrom, int(start), int(end)))
            elif len(groups) == 5:
                # Format: chr_refstart_refend:start-end
                chrom, ref_start, ref_end, start, end = groups
                regions.append((chrom, int(ref_start), int(ref_end)))

    return regions


def extract_crIDs_from_log_line(line: str) -> Set[int]:
    """
    Extract candidate region IDs (crIDs) from a log line.

    Returns:
        Set of crIDs found in the line
    """
    crIDs = set()

    # Pattern: crID=N or crIDs={N1,N2,...} or crIDs=[N1,N2,...]
    patterns = [
        r"crID=(\d+)",
        r"crIDs=\{([0-9,\s]+)\}",
        r"crIDs=\[([0-9,\s]+)\]",
    ]

    for pattern in patterns:
        matches = re.finditer(pattern, line)
        for match in matches:
            value = match.group(1)
            if "," in value:
                # Multiple crIDs
                for crid_str in value.split(","):
                    crid_str = crid_str.strip()
                    if crid_str:
                        crIDs.add(int(crid_str))
            else:
                # Single crID
                crIDs.add(int(value))

    return crIDs


def build_region_index(log_file: Path) -> Tuple[Dict[str, IntervalTree], List[str]]:
    """
    Build an interval tree index from the log file based on genomic regions.

    Returns:
        Tuple of (interval_trees dict, all_lines list)
    """
    interval_trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
    all_lines: List[str] = []

    with open(log_file, "r") as f:
        for line_idx, line in enumerate(f):
            all_lines.append(line.rstrip())

            # Extract regions from the line
            regions = extract_region_from_log_line(line)

            for chrom, start, end in regions:
                # Handle point intervals (start == end) by extending them by 1
                # IntervalTree doesn't allow null intervals
                if start == end:
                    end = start + 1

                # Add to interval tree with line index as data
                interval_trees[chrom].addi(start, end, line_idx)

    return interval_trees, all_lines


def query_by_region(log_file: Path, region: str, margin: int) -> List[Tuple[int, str]]:
    """
    Query log file by genomic region.

    Returns:
        List of (line_index, line_content) tuples
    """
    chrom, start, end = parse_region(region)

    # Apply margin
    start = max(0, start - margin)
    end = end + margin

    print(f"Querying region: {chrom}:{start}-{end}", file=sys.stderr)

    # Build interval tree index
    interval_trees, all_lines = build_region_index(log_file)

    if chrom not in interval_trees:
        print(f"No entries found for chromosome: {chrom}", file=sys.stderr)
        return []

    # Query the interval tree
    overlapping_intervals = interval_trees[chrom].overlap(start, end)

    # Collect unique line indices
    line_indices = set()
    for interval in overlapping_intervals:
        line_indices.add(interval.data)

    # Sort by line index to maintain order
    sorted_indices = sorted(line_indices)

    return [(idx, all_lines[idx]) for idx in sorted_indices]


def query_by_crIDs(log_file: Path, crIDs: Set[int]) -> List[Tuple[int, str]]:
    """
    Query log file by candidate region IDs.

    Returns:
        List of (line_index, line_content) tuples
    """
    print(f"Querying crIDs: {sorted(crIDs)}", file=sys.stderr)

    matching_lines = []

    with open(log_file, "r") as f:
        for line_idx, line in enumerate(f):
            line_crIDs = extract_crIDs_from_log_line(line)

            # Check if any of the query crIDs are in this line
            if line_crIDs.intersection(crIDs):
                matching_lines.append((line_idx, line.rstrip()))

    return matching_lines


def main():
    """Main function."""
    args = parse_arguments()

    if not args.log_file.exists():
        print(f"Error: Log file not found: {args.log_file}", file=sys.stderr)
        sys.exit(1)

    # Collect matching lines
    matching_lines: Dict[int, str] = {}  # line_index -> line_content

    # Query by region if provided
    if args.region:
        region_matches = query_by_region(args.log_file, args.region, args.region_margin)
        for idx, line in region_matches:
            matching_lines[idx] = line
        print(f"Found {len(region_matches)} lines matching region", file=sys.stderr)

    # Query by crIDs if provided
    if args.crIDs:
        crID_set = {int(x.strip()) for x in args.crIDs.split(",")}
        crID_matches = query_by_crIDs(args.log_file, crID_set)
        for idx, line in crID_matches:
            matching_lines[idx] = line
        print(f"Found {len(crID_matches)} lines matching crIDs", file=sys.stderr)

    # Print matching lines in order
    if not matching_lines:
        print("No matching lines found", file=sys.stderr)
        return

    print(f"\n{'=' * 80}", file=sys.stderr)
    print(f"Total unique matching lines: {len(matching_lines)}", file=sys.stderr)
    print(f"{'=' * 80}\n", file=sys.stderr)

    for idx in sorted(matching_lines.keys()):
        print(matching_lines[idx])


if __name__ == "__main__":
    main()
