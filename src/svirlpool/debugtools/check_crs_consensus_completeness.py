#!/usr/bin/env python
"""
Check candidate regions consensus completeness.

This script takes a CRs BED file and a CRS container results file to check for completeness,
meaning it identifies which candidate regions successfully produced consensus sequences
and which ones failed to do so.

The BED file should have columns: chrom, start, end, crID
where each line represents a single candidate region.

The CRS container results file contains CrsContainerResult objects (one per line as JSON),
which store Consensus objects. Each Consensus object has a crIDs field listing all
candidate region IDs that contributed to creating that consensus object.

This allows us to track:
- Which candidate regions produced consensus sequences
- Which candidate regions were merged together (multiple crIDs in one Consensus)
- Which candidate regions failed to produce any consensus
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

# Import the necessary classes and functions from the svirlpool package
from ..localassembly import consensus_align

log = logging.getLogger(__name__)


@dataclass
class CRConsensusMatch:
    """Results of matching a candidate region against consensus sequences."""

    cr_id: int
    cr_region: str
    found: bool
    consensus_count: int
    consensus_ids: list[str]
    merged_with_crs: set[int]  # Other CR IDs that were merged with this one

    def __str__(self) -> str:
        """Format the match result as a string."""
        status = "FOUND" if self.found else "MISSING"

        # Format merged CRs
        merged_str = "N/A"
        if self.merged_with_crs:
            merged_str = ",".join(str(cid) for cid in sorted(self.merged_with_crs))

        return (
            f"{self.cr_id}\t{self.cr_region}\t{status}\t"
            f"{self.consensus_count}\t{merged_str}\t"
            f"{','.join(self.consensus_ids) if self.consensus_ids else 'N/A'}"
        )


@dataclass
class CandidateRegion:
    """Simple representation of a candidate region from BED file."""

    cr_id: int
    chrom: str
    start: int
    end: int

    def region_string(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def load_candidate_regions_from_bed(bed_path: Path) -> list[CandidateRegion]:
    """Load candidate regions from BED file.

    Expected BED format: chrom, start, end, crID
    """
    crs: list[CandidateRegion] = []

    with open(bed_path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 4:
                log.warning(
                    f"Line {line_num}: Expected at least 4 fields, got {len(fields)}"
                )
                continue

            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                cr_id = int(fields[3])

                crs.append(
                    CandidateRegion(
                        cr_id=cr_id,
                        chrom=chrom,
                        start=start,
                        end=end,
                    )
                )
            except (ValueError, IndexError) as e:
                log.warning(f"Line {line_num}: Failed to parse: {e}")
                continue

    log.info(f"Loaded {len(crs)} candidate regions from BED file")
    return crs


def load_consensus_from_container_results(
    container_results_path: Path,
) -> tuple[dict[int, list[str]], dict[str, set[int]]]:
    """Load consensus objects from CRS container results file.

    Returns:
        Tuple of:
        - Dictionary mapping CR ID to list of consensus IDs
        - Dictionary mapping consensus ID to set of all CR IDs that contributed to it
    """
    consensus_by_cr: dict[int, list[str]] = {}
    consensus_to_crs: dict[str, set[int]] = {}

    log.info(f"Parsing CRS container results from {container_results_path}")

    for crs_container in consensus_align.parse_crs_container_results(
        container_results_path
    ):
        # Process each Consensus object in this container
        for consensus_id, consensus_obj in crs_container.consensus_dicts.items():
            # Record the consensus ID and which CRs contributed to it
            consensus_to_crs[consensus_id] = set(consensus_obj.crIDs)

            # For each CR ID that contributed, record this consensus
            for cr_id in consensus_obj.crIDs:
                if cr_id not in consensus_by_cr:
                    consensus_by_cr[cr_id] = []
                consensus_by_cr[cr_id].append(consensus_id)

    total_consensus = sum(len(v) for v in consensus_by_cr.values())
    unique_consensus = len(consensus_to_crs)

    log.info(
        f"Loaded {unique_consensus} unique consensus sequences "
        f"(total {total_consensus} CR-to-consensus mappings) "
        f"for {len(consensus_by_cr)} candidate regions"
    )
    return consensus_by_cr, consensus_to_crs


def compare_crs_to_consensus(
    crs_bed_path: Path, container_results_path: Path
) -> tuple[list[CRConsensusMatch], dict[str, int]]:
    """Compare candidate regions to consensus sequences and return detailed results.

    This function is suitable for use in unit tests.

    Args:
        crs_bed_path: Path to candidate regions BED file
        container_results_path: Path to CRS container results file

    Returns:
        Tuple of (results, statistics) where:
        - results: List of CRConsensusMatch objects for each candidate region
        - statistics: Dict with keys 'total', 'found', 'missing', 'total_consensus', 'merged'
    """
    candidate_regions = load_candidate_regions_from_bed(crs_bed_path)
    consensus_by_cr, consensus_to_crs = load_consensus_from_container_results(
        container_results_path
    )

    results: list[CRConsensusMatch] = []

    for cr in candidate_regions:
        cr_id = cr.cr_id
        cr_region = cr.region_string()

        # Get consensus IDs for this CR
        consensus_ids = consensus_by_cr.get(cr_id, [])
        found = len(consensus_ids) > 0

        # Determine which other CRs this one was merged with
        merged_with_crs = set()
        for consensus_id in consensus_ids:
            # Get all CRs that contributed to this consensus
            contributing_crs = consensus_to_crs.get(consensus_id, set())
            # Add all other CRs (excluding the current one)
            merged_with_crs.update(crid for crid in contributing_crs if crid != cr_id)

        match = CRConsensusMatch(
            cr_id=cr_id,
            cr_region=cr_region,
            found=found,
            consensus_count=len(consensus_ids),
            consensus_ids=consensus_ids,
            merged_with_crs=merged_with_crs,
        )
        results.append(match)

    # Calculate statistics
    total = len(results)
    found = sum(1 for r in results if r.found)
    missing = total - found
    total_consensus = sum(r.consensus_count for r in results)
    merged = sum(1 for r in results if r.merged_with_crs)

    statistics = {
        "total": total,
        "found": found,
        "missing": missing,
        "total_consensus": total_consensus,
        "merged": merged,
    }

    return results, statistics


def check_crs_consensus_completeness(
    crs_bed_path: Path, container_results_path: Path, output_path: Path | None = None
) -> None:
    """Check candidate regions consensus completeness."""
    log.info(f"Reading candidate regions from {crs_bed_path}")
    log.info(f"Reading consensus from container results {container_results_path}")

    results, statistics = compare_crs_to_consensus(crs_bed_path, container_results_path)

    total = statistics["total"]
    found = statistics["found"]
    missing = statistics["missing"]
    total_consensus = statistics["total_consensus"]
    merged = statistics["merged"]

    # Output results
    output_file = open(output_path, "w") if output_path else sys.stdout

    try:
        # Write header
        output_file.write(
            "CR_ID\tREGION\tSTATUS\tCONSENSUS_COUNT\tMERGED_WITH_CRS\tCONSENSUS_IDS\n"
        )

        # Write results
        for result in results:
            output_file.write(str(result) + "\n")
            # Print missing CRs to terminal
            if not result.found:
                log.warning(
                    f"Missing consensus for CR {result.cr_id}: {result.cr_region}"
                )
            # Print merged CRs to terminal
            # elif result.merged_with_crs:
            #     merged_str = ",".join(str(cid) for cid in sorted(result.merged_with_crs))
            #     log.info(
            #         f"CR {result.cr_id} was merged with CR(s) {merged_str}"
            #     )

        # Summary statistics
        output_file.write("\n# SUMMARY\n")
        output_file.write(f"# Total candidate regions: {total}\n")
        found_percentage = (found / total * 100) if total > 0 else 0
        output_file.write(f"# With consensus: {found} ({found_percentage:.1f}%)\n")
        missing_percentage = (missing / total * 100) if total > 0 else 0
        output_file.write(
            f"# Without consensus: {missing} ({missing_percentage:.1f}%)\n"
        )
        merged_percentage = (merged / total * 100) if total > 0 else 0
        output_file.write(
            f"# Merged with other CRs: {merged} ({merged_percentage:.1f}%)\n"
        )
        output_file.write(f"# Total consensus sequences: {total_consensus}\n")
        average_consensus = (total_consensus / found) if found > 0 else 0
        output_file.write(
            f"# Average consensus per CR (with consensus): {average_consensus:.2f}\n"
        )

        log.info(f"Candidate regions with consensus: {found}/{total}")
        log.info(f"Candidate regions merged with others: {merged}")
        log.info(f"Total consensus sequences: {total_consensus}")

    finally:
        if output_path:
            output_file.close()


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Check candidate regions consensus completeness.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--crs-bed",
        type=Path,
        required=True,
        help="Path to candidate regions BED file (chrom, start, end, crID)",
    )
    parser.add_argument(
        "-c",
        "--container-results",
        type=Path,
        required=True,
        help="Path to CRS container results file (JSON lines format)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Path to output file (default: stdout)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    return parser


def main() -> None:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    check_crs_consensus_completeness(
        crs_bed_path=args.crs_bed,
        container_results_path=args.container_results,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
