#!/usr/bin/env python
"""
Check candidate regions consensus completeness.

This script takes a CRs BED file and a consensus FASTA file and checks for completeness,
meaning it counts the consensus sequences for each candidate region and lists the candidate
regions that did not produce any consensus.

The BED file should have columns: chrom, start, end, cr_ids
where cr_ids can be comma-separated when multiple CRs are merged into one region.

The consensus sequences are represented in the record name of each record in the FASTA,
with the prefix being the candidate region ID. The suffix is an index of the consensus
sequence for that candidate region.

Example:
    A candidate region with ID 13 can have consensus sequences with IDs:
    "13.0", "13.1", "13.2", etc.

    Note: Some candidate regions may be merged, so multiple CR IDs can share the same
    consensus sequence. In the BED file, this is represented as comma-separated IDs
    (e.g., "12,13").
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

log = logging.getLogger(__name__)


@dataclass
class CRConsensusMatch:
    """Results of matching a candidate region against consensus sequences."""

    cr_ids: list[int]  # Can be multiple IDs if merged
    cr_region: str
    found: bool
    consensus_count: int
    consensus_ids: list[str]

    def __str__(self) -> str:
        """Format the match result as a string."""
        status = "FOUND" if self.found else "MISSING"
        cr_ids_str = ",".join(str(cid) for cid in self.cr_ids)
        return (
            f"{cr_ids_str}\t{self.cr_region}\t{status}\t"
            f"{self.consensus_count}\t{','.join(self.consensus_ids) if self.consensus_ids else 'N/A'}"
        )


@dataclass
class CandidateRegion:
    """Simple representation of a candidate region from BED file."""

    cr_ids: list[int]  # Can be multiple IDs if merged
    chrom: str
    start: int
    end: int

    def region_string(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    def ids_string(self) -> str:
        return ",".join(str(cid) for cid in self.cr_ids)


def load_candidate_regions_from_bed(bed_path: Path) -> list[CandidateRegion]:
    """Load candidate regions from BED file.

    Expected BED format: chrom, start, end, cr_ids (comma-separated)
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
                # Parse comma-separated CR IDs
                cr_ids_str = fields[3]
                cr_ids = [int(cid.strip()) for cid in cr_ids_str.split(",")]

                crs.append(
                    CandidateRegion(
                        cr_ids=cr_ids,
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


def load_consensus_sequences(consensus_fasta_path: Path) -> dict[int, list[str]]:
    """Load consensus sequences from FASTA file and group by CR ID.

    Returns:
        Dictionary mapping CR ID to list of consensus sequence IDs
    """
    consensus_by_cr: dict[int, list[str]] = {}

    # Handle both gzipped and regular FASTA files
    if str(consensus_fasta_path).endswith(".gz"):
        handle = gzip.open(consensus_fasta_path, "rt")
    else:
        handle = open(consensus_fasta_path, "r")

    try:
        for record in SeqIO.parse(handle, "fasta"):
            # Parse consensus ID (e.g., "13.0" -> CR ID 13, consensus index 0)
            consensus_id = record.id
            try:
                # Split by '.' to get CR ID
                cr_id_str = consensus_id.split(".")[0]
                cr_id = int(cr_id_str)

                if cr_id not in consensus_by_cr:
                    consensus_by_cr[cr_id] = []
                consensus_by_cr[cr_id].append(consensus_id)
            except (ValueError, IndexError) as e:
                log.warning(f"Could not parse consensus ID '{consensus_id}': {e}")
                continue
    finally:
        handle.close()

    log.info(
        f"Loaded {sum(len(v) for v in consensus_by_cr.values())} consensus sequences "
        f"for {len(consensus_by_cr)} candidate regions"
    )
    return consensus_by_cr


def compare_crs_to_consensus(
    crs_bed_path: Path, consensus_fasta_path: Path
) -> tuple[list[CRConsensusMatch], dict[str, int]]:
    """Compare candidate regions to consensus sequences and return detailed results.

    This function is suitable for use in unit tests.

    Args:
        crs_bed_path: Path to candidate regions BED file
        consensus_fasta_path: Path to consensus FASTA file

    Returns:
        Tuple of (results, statistics) where:
        - results: List of CRConsensusMatch objects for each candidate region
        - statistics: Dict with keys 'total', 'found', 'missing', 'total_consensus'
    """
    candidate_regions = load_candidate_regions_from_bed(crs_bed_path)
    consensus_by_cr = load_consensus_sequences(consensus_fasta_path)

    results: list[CRConsensusMatch] = []

    for cr in candidate_regions:
        cr_ids = cr.cr_ids
        cr_region = cr.region_string()

        # Collect consensus IDs for all CR IDs in this region (in case of merged CRs)
        all_consensus_ids = []
        for cr_id in cr_ids:
            all_consensus_ids.extend(consensus_by_cr.get(cr_id, []))

        # Remove duplicates while preserving order
        seen = set()
        unique_consensus_ids = []
        for cid in all_consensus_ids:
            if cid not in seen:
                seen.add(cid)
                unique_consensus_ids.append(cid)

        found = len(unique_consensus_ids) > 0

        match = CRConsensusMatch(
            cr_ids=cr_ids,
            cr_region=cr_region,
            found=found,
            consensus_count=len(unique_consensus_ids),
            consensus_ids=unique_consensus_ids,
        )
        results.append(match)

    # Calculate statistics
    total = len(results)
    found = sum(1 for r in results if r.found)
    missing = total - found
    total_consensus = sum(r.consensus_count for r in results)

    statistics = {
        "total": total,
        "found": found,
        "missing": missing,
        "total_consensus": total_consensus,
    }

    return results, statistics


def check_crs_consensus_completeness(
    crs_bed_path: Path, consensus_fasta_path: Path, output_path: Path | None = None
) -> None:
    """Check candidate regions consensus completeness."""
    log.info(f"Reading candidate regions from {crs_bed_path}")
    log.info(f"Reading consensus sequences from {consensus_fasta_path}")

    results, statistics = compare_crs_to_consensus(crs_bed_path, consensus_fasta_path)

    total = statistics["total"]
    found = statistics["found"]
    missing = statistics["missing"]
    total_consensus = statistics["total_consensus"]

    # Output results
    output_file = open(output_path, "w") if output_path else sys.stdout

    try:
        # Write header
        output_file.write("CR_ID\tREGION\tSTATUS\tCONSENSUS_COUNT\tCONSENSUS_IDS\n")

        # Write results
        for result in results:
            output_file.write(str(result) + "\n")
            # Print missing CRs to terminal
            if not result.found:
                cr_ids_str = ",".join(str(cid) for cid in result.cr_ids)
                log.warning(
                    f"Missing consensus for CR(s) {cr_ids_str}: {result.cr_region}"
                )

        # Summary statistics
        output_file.write("\n# SUMMARY\n")
        output_file.write(f"# Total candidate regions: {total}\n")
        output_file.write(f"# With consensus: {found} ({found / total * 100:.1f}%)\n")
        output_file.write(
            f"# Without consensus: {missing} ({missing / total * 100:.1f}%)\n"
        )
        output_file.write(f"# Total consensus sequences: {total_consensus}\n")
        output_file.write(
            f"# Average consensus per CR: {total_consensus / total:.2f}\n"
        )

        log.info(f"Candidate regions with consensus: {found}/{total}")
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
        help="Path to candidate regions BED file (chrom, start, end, cr_ids)",
    )
    parser.add_argument(
        "-c",
        "--consensus-fasta",
        type=Path,
        required=True,
        help="Path to consensus FASTA file (can be gzipped)",
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
        consensus_fasta_path=args.consensus_fasta,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
