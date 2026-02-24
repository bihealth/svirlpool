#!/usr/bin/env python
"""
Search log files for DROPPED or TRANSFORMED lines matching given crIDs or regions.

This script reads a log file (produced with logging level DEBUG) and searches for lines
starting with DROPPED or TRANSFORMED. It then checks whether any of the user-specified
crIDs or regions appear in those lines.

crIDs are identified in log lines by the patterns:
    crID=<int>
    crIDs={<int>,<int>,...}
    crIDs=[<int>, <int>, ...]

Regions are identified by the patterns:
    region=<chr>:<start>-<end>
    regions=<chr>:<start>-<end>;<chr>:<start>-<end>
    regions=<chr>:<start>-<end>,<chr>:<start>-<end>

Matching logic:
    - crID match: set intersection between query crIDs and line crIDs
    - Region match: interval overlap between any query region and any line region
      (same chromosome, intervals overlap)

If at least one crID or region matches, the line is printed.

Usage:
    python -m svirlpool.debugtools.log_search -l <logfile> [--crIDs 1 2 3] [--regions chr1:100-200 chr2:300-400]
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Region:
    """A genomic region with chromosome, start, and end."""

    chrom: str
    start: int
    end: int

    def overlaps(self, other: Region) -> bool:
        """Check if two regions on the same chromosome overlap."""
        if self.chrom != other.chrom:
            return False
        return self.start < other.end and other.start < self.end

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def parse_region(region_str: str) -> Region:
    """Parse a region string like chr1:100-200 into a Region object.

    Handles commas in coordinates (e.g., chr1:1,000-2,000).
    """
    # Split on the last colon to handle chromosome names with colons
    colon_idx = region_str.rfind(":")
    if colon_idx == -1:
        raise ValueError(f"Invalid region format (no colon): {region_str}")
    chrom = region_str[:colon_idx]
    coords = region_str[colon_idx + 1 :]
    if "-" not in coords:
        raise ValueError(
            f"Invalid region format (no dash in coordinates): {region_str}"
        )
    start_str, end_str = coords.split("-", 1)
    start = int(start_str.replace(",", ""))
    end = int(end_str.replace(",", ""))
    return Region(chrom=chrom, start=start, end=end)


# Regex patterns for extracting crIDs from log lines
# Matches: crID=123
_RE_CRID_SINGLE = re.compile(r"crID=(\d+)")
# Matches: crIDs={1,2,3} or crIDs=[1, 2, 3] or crIDs={1, 2, 3}
_RE_CRIDS_SET = re.compile(r"crIDs=[\[{]([\d,\s]+)[\]}]")

# Regex patterns for extracting regions from log lines
# Matches: region=chr1:100-200
_RE_REGION_SINGLE = re.compile(r"region=([\w._-]+:\d[\d,]*-\d[\d,]*)")
# Matches: regions=chr1:100-200;chr2:300-400 or regions=chr1:100-200,chr2:300-400
# Also handles the format: regions=chr1:100-200 chr2:300-400 (space-separated with region= prefix)
_RE_REGIONS_BLOCK = re.compile(
    r"regions=([\w._-]+:\d[\d,]*-\d[\d,]*(?:[;,]\s*[\w._-]+:\d[\d,]*-\d[\d,]*)*)"
)


def extract_crIDs_from_line(line: str) -> set[int]:
    """Extract all crIDs mentioned in a log line."""
    crIDs: set[int] = set()

    # Extract crIDs from set/list notation: crIDs={1,2,3} or crIDs=[1, 2, 3]
    for match in _RE_CRIDS_SET.finditer(line):
        for num_str in match.group(1).split(","):
            num_str = num_str.strip()
            if num_str.isdigit():
                crIDs.add(int(num_str))

    # Extract single crID= values
    for match in _RE_CRID_SINGLE.finditer(line):
        crIDs.add(int(match.group(1)))

    return crIDs


def extract_regions_from_line(line: str) -> list[Region]:
    """Extract all regions mentioned in a log line."""
    regions: list[Region] = []
    seen: set[str] = set()

    # Extract regions= blocks (multi-region)
    for match in _RE_REGIONS_BLOCK.finditer(line):
        block = match.group(1)
        # Split on ; or , but be careful: commas can appear inside coordinates
        # Strategy: split on ; first. If that gives only one item, try splitting
        # by looking for chrom:start-end patterns.
        parts = re.findall(r"[\w._-]+:\d[\d,]*-\d[\d,]*", block)
        for part in parts:
            if part not in seen:
                try:
                    regions.append(parse_region(part))
                    seen.add(part)
                except ValueError:
                    pass

    # Extract single region= values
    for match in _RE_REGION_SINGLE.finditer(line):
        region_str = match.group(1)
        if region_str not in seen:
            try:
                regions.append(parse_region(region_str))
                seen.add(region_str)
            except ValueError:
                pass

    return regions


def line_matches(line: str, query_crIDs: set[int], query_regions: list[Region]) -> bool:
    """Check if a log line matches any of the query crIDs or regions.

    A match occurs if:
    - Any query crID is found in the line's crIDs (set intersection), OR
    - Any query region overlaps with any region in the line (interval overlap)
    """
    # If no queries specified, nothing matches
    if not query_crIDs and not query_regions:
        return False

    # Check crID match (set intersection)
    if query_crIDs:
        line_crIDs = extract_crIDs_from_line(line)
        if line_crIDs & query_crIDs:
            return True

    # Check region overlap
    if query_regions:
        line_regions = extract_regions_from_line(line)
        for q_region in query_regions:
            for l_region in line_regions:
                if q_region.overlaps(l_region):
                    return True

    return False


def is_relevant_line(line: str) -> bool:
    """Check if a line is a DROPPED or TRANSFORMED log line.

    The DROPPED/TRANSFORMED keyword can appear after the logging prefix
    (e.g. timestamp, logger name, level), so we search for it anywhere in the line.
    """
    return "DROPPED" in line or "TRANSFORMED" in line


def search_log(
    log_path: Path,
    query_crIDs: set[int],
    query_regions: list[Region],
    show_all_relevant: bool = False,
) -> list[str]:
    """Search a log file for matching DROPPED/TRANSFORMED lines.

    Args:
        log_path: Path to log file (plain text or gzipped)
        query_crIDs: Set of crIDs to search for
        query_regions: List of regions to search for (interval overlap)
        show_all_relevant: If True, show all DROPPED/TRANSFORMED lines, not just matches

    Returns:
        List of matching lines
    """
    matches: list[str] = []

    if str(log_path).endswith(".gz"):
        fh: io.TextIOBase = gzip.open(log_path, "rt")  # type: ignore[assignment]
    else:
        fh = open(log_path, "r")  # type: ignore[assignment]

    with fh:
        for raw_line in fh:
            line = str(raw_line).rstrip("\n")
            if not is_relevant_line(line):
                continue
            if show_all_relevant or line_matches(line, query_crIDs, query_regions):
                matches.append(line)

    return matches


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Search log files for DROPPED/TRANSFORMED lines matching given crIDs or regions. "
            "Lines are matched if they contain an overlapping crID (set overlap) or an overlapping "
            "genomic region (interval overlap on the same chromosome)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Search by crIDs\n"
            "  python -m svirlpool.debugtools.log_search -l run.log --crIDs 5 19 42\n\n"
            "  # Search by region\n"
            "  python -m svirlpool.debugtools.log_search -l run.log --regions chr1:1000-2000\n\n"
            "  # Search by both crIDs and regions\n"
            "  python -m svirlpool.debugtools.log_search -l run.log --crIDs 5 --regions chr1:1000-2000 chr2:500-800\n\n"
            "  # Search a gzipped log\n"
            "  python -m svirlpool.debugtools.log_search -l run.log.gz --crIDs 5\n\n"
            "  # Show all DROPPED/TRANSFORMED lines (no filtering)\n"
            "  python -m svirlpool.debugtools.log_search -l run.log --all\n"
        ),
    )

    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        required=True,
        help="Path to log file (plain text or gzipped)",
    )
    parser.add_argument(
        "--crIDs",
        type=int,
        nargs="+",
        default=[],
        help="crID(s) to search for (integers). A line matches if any of its crIDs overlap with these.",
    )
    parser.add_argument(
        "--regions",
        type=str,
        nargs="+",
        default=[],
        help="Region(s) to search for, formatted as chr:start-end. A line matches if any of its regions overlap with these.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        dest="show_all",
        help="Show all DROPPED/TRANSFORMED lines without filtering by crID or region.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Path to output file. If not provided, results are printed to stdout.",
    )

    return parser


def main() -> None:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args()

    # Parse query regions
    query_regions: list[Region] = []
    for region_str in args.regions:
        try:
            query_regions.append(parse_region(region_str))
        except ValueError as e:
            print(f"Error parsing region '{region_str}': {e}", file=sys.stderr)
            sys.exit(1)

    query_crIDs: set[int] = set(args.crIDs)

    if not query_crIDs and not query_regions and not args.show_all:
        print(
            "Error: At least one of --crIDs, --regions, or --all must be provided.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check that the log file exists
    if not args.log.exists():
        print(f"Error: Log file not found: {args.log}", file=sys.stderr)
        sys.exit(1)

    # Print search parameters
    print(f"Searching log: {args.log}", file=sys.stderr)
    if query_crIDs:
        print(f"  Query crIDs: {sorted(query_crIDs)}", file=sys.stderr)
    if query_regions:
        print(f"  Query regions: {[str(r) for r in query_regions]}", file=sys.stderr)
    if args.show_all:
        print("  Showing all DROPPED/TRANSFORMED lines", file=sys.stderr)

    # Search
    matches = search_log(
        log_path=args.log,
        query_crIDs=query_crIDs,
        query_regions=query_regions,
        show_all_relevant=args.show_all,
    )

    # Output
    output_file = open(args.output, "w") if args.output else sys.stdout
    try:
        for line in matches:
            output_file.write(line + "\n")

        print(f"\nFound {len(matches)} matching lines.", file=sys.stderr)
    finally:
        if args.output:
            output_file.close()


if __name__ == "__main__":
    main()
