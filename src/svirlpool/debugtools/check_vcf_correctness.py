#!/usr/bin/env python
"""
Check VCF correctness against expected variants.

This script takes as input a (bgzipped) VCF file and compares its variants to the expected
variants from a regions + descriptors file.

Expected variants file format:
    chrom:start-end    SVTYPE    size    HET/HOM

Example:
    R8_14_91775678_92309099:247,964-248,188    DEL    98    HOM
    R8_14_91775678_92309099:325,398-325,407    INS    2313    HET

Within each region, the script checks if a variant is present, and if this variant is of
the expected type, the expected size (within 20% margin), and the expected genotype (HET/HOM).
For each expected variant, it outputs whether it was found in the VCF and whether it matches
the expected type, size, and genotype.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import vcfpy

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class ExpectedVariant:
    """Represents an expected variant from the regions file."""

    chrom: str
    start: int
    end: int
    svtype: str
    size: int
    genotype: str  # HET or HOM

    @classmethod
    def from_line(cls, line: str) -> ExpectedVariant:
        """Parse an expected variant from a line in the regions file."""
        # Split on tabs first, but if that doesn't work, split on whitespace
        fields = line.strip().split("\t")
        if len(fields) < 4:
            # Try splitting on any whitespace
            fields = line.strip().split()
        if len(fields) < 4:
            raise ValueError(f"Invalid line format: {line}")

        # Parse chrom:start-end
        region_str = fields[0]
        chrom, coords = region_str.rsplit(":", 1)
        start_str, end_str = coords.split("-")
        # Remove commas from coordinates
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))

        svtype = fields[1]
        size = int(fields[2])
        genotype = fields[3]

        return cls(
            chrom=chrom,
            start=start,
            end=end,
            svtype=svtype,
            size=size,
            genotype=genotype,
        )


@dataclass
class VariantMatch:
    """Results of matching an expected variant against VCF."""

    expected: ExpectedVariant
    found: bool
    type_match: bool
    size_match: bool
    genotype_match: bool
    vcf_variant: vcfpy.Record | None = None

    def __str__(self) -> str:
        """Format the match result as a string."""
        status = "FOUND" if self.found else "MISSING"
        details = []

        if self.found:
            if self.type_match:
                details.append("TYPE:OK")
            else:
                details.append("TYPE:MISMATCH")

            if self.size_match:
                details.append("SIZE:OK")
            else:
                details.append("SIZE:MISMATCH")

            if self.genotype_match:
                details.append("GT:OK")
            else:
                details.append("GT:MISMATCH")

        region = f"{self.expected.chrom}:{self.expected.start}-{self.expected.end}"
        details_str = ",".join(details) if details else "N/A"

        return (
            f"{region}\t{self.expected.svtype}\t{self.expected.size}\t"
            f"{self.expected.genotype}\t{status}\t{details_str}"
        )


def parse_expected_variants(path: Path) -> list[ExpectedVariant]:
    """Parse expected variants from the regions file."""
    variants: list[ExpectedVariant] = []

    # First pass: count tabs in each line to find the majority
    tab_counts: dict[int, list[int]] = {}  # tab_count -> [line_numbers]
    lines_data: list[tuple[int, str]] = []  # (line_number, line_content)

    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            lines_data.append((line_num, line))
            tab_count = line.count("\t")
            if tab_count not in tab_counts:
                tab_counts[tab_count] = []
            tab_counts[tab_count].append(line_num)

    # Find the majority tab count
    if tab_counts:
        majority_tab_count = max(tab_counts.keys(), key=lambda k: len(tab_counts[k]))

        # Report lines that deviate from the majority
        for tab_count, line_numbers in tab_counts.items():
            if tab_count != majority_tab_count:
                log.warning(
                    f"Lines with {tab_count} tabs (expected {majority_tab_count}): "
                    f"{', '.join(map(str, line_numbers))}"
                )

    # Second pass: parse the variants
    for line_num, line in lines_data:
        try:
            variant = ExpectedVariant.from_line(line)
            variants.append(variant)
        except Exception as e:
            log.warning(f"Failed to parse line {line_num}: {line}. Error: {e}")
            continue

    log.info(f"Parsed {len(variants)} expected variants")
    return variants


def get_genotype_from_vcf(variant: vcfpy.Record) -> str:
    """Extract genotype (HET/HOM) from a VCF record.

    Returns "HET" for heterozygous, "HOM" for homozygous alt, "REF" for homozygous ref.
    """
    # Get the first sample (assuming single-sample VCF)
    if not variant.calls:
        return "UNKNOWN"

    call = variant.calls[0]
    try:
        gt = call.data["GT"]  # type: ignore
    except (KeyError, AttributeError):
        return "UNKNOWN"

    if not gt:
        return "UNKNOWN"

    # Parse the genotype string
    if "/" in gt:
        alleles = gt.split("/")
    elif "|" in gt:
        alleles = gt.split("|")
    else:
        return "UNKNOWN"

    # Filter out missing alleles (.)
    alleles = [a for a in alleles if a != "."]

    if not alleles:
        return "UNKNOWN"

    # Count reference (0) and alt (non-0) alleles
    alt_count = sum(1 for a in alleles if a != "0")
    ref_count = sum(1 for a in alleles if a == "0")

    if alt_count == 0:
        return "REF"
    elif ref_count > 0:
        return "HET"
    else:
        return "HOM"


def variant_overlaps_region(
    variant: vcfpy.Record, chrom: str, start: int, end: int
) -> bool:
    """Check if a variant overlaps with the given region."""
    if variant.CHROM != chrom:
        return False

    var_start = variant.POS
    var_end = variant.INFO.get("END", var_start + 1)
    if isinstance(var_end, list):
        var_end = var_end[0]

    # Check if intervals overlap
    return var_start < end and var_end > start


def find_matching_variants(
    expected: ExpectedVariant, vcf_variants: list[vcfpy.Record]
) -> VariantMatch:
    """Find and compare variants in the VCF that match the expected variant."""
    # Filter variants that overlap with the expected region
    overlapping_variants = [
        v
        for v in vcf_variants
        if variant_overlaps_region(v, expected.chrom, expected.start, expected.end)
    ]

    if not overlapping_variants:
        return VariantMatch(
            expected=expected,
            found=False,
            type_match=False,
            size_match=False,
            genotype_match=False,
        )

    # Find the best matching variant
    best_match = None
    best_score = 0

    for vcf_var in overlapping_variants:
        # Check SVTYPE
        vcf_svtype = vcf_var.INFO.get("SVTYPE", "UNKNOWN")
        type_match = vcf_svtype == expected.svtype

        # Check size (within 20% margin)
        vcf_svlen = vcf_var.INFO.get("SVLEN", 0)
        if isinstance(vcf_svlen, list):
            vcf_svlen = vcf_svlen[0]
        vcf_size = abs(int(vcf_svlen))

        size_margin = expected.size * 0.2
        size_match = abs(vcf_size - expected.size) <= size_margin

        # Check genotype
        vcf_genotype = get_genotype_from_vcf(vcf_var)
        genotype_match = vcf_genotype == expected.genotype

        # Score this match
        score = sum([type_match, size_match, genotype_match])

        if score > best_score:
            best_score = score
            best_match = VariantMatch(
                expected=expected,
                found=True,
                type_match=type_match,
                size_match=size_match,
                genotype_match=genotype_match,
                vcf_variant=vcf_var,
            )

    if best_match is None:
        return VariantMatch(
            expected=expected,
            found=False,
            type_match=False,
            size_match=False,
            genotype_match=False,
        )

    return best_match


def compare_vcf_to_expected(
    vcf_path: Path, expected_path: Path
) -> tuple[list[VariantMatch], dict[str, int]]:
    """Compare VCF variants to expected variants and return detailed results.

    This function is suitable for use in unit tests.

    Args:
        vcf_path: Path to VCF file
        expected_path: Path to expected variants file

    Returns:
        Tuple of (results, statistics) where:
        - results: List of VariantMatch objects for each expected variant
        - statistics: Dict with keys 'total', 'found', 'type_ok', 'size_ok', 'gt_ok', 'all_ok'
    """
    expected_variants = parse_expected_variants(expected_path)
    vcf_reader = vcfpy.Reader.from_path(str(vcf_path))

    # Read all VCF variants and group by chromosome
    vcf_variants_by_chrom: dict[str, list[vcfpy.Record]] = {}
    for variant in vcf_reader:
        if variant is None:
            continue
        chrom = variant.CHROM
        if chrom not in vcf_variants_by_chrom:
            vcf_variants_by_chrom[chrom] = []
        vcf_variants_by_chrom[chrom].append(variant)

    results: list[VariantMatch] = []

    for expected in expected_variants:
        chrom_variants = vcf_variants_by_chrom.get(expected.chrom, [])
        match = find_matching_variants(expected, chrom_variants)
        results.append(match)

    # Calculate statistics
    total = len(results)
    found = sum(1 for r in results if r.found)
    type_ok = sum(1 for r in results if r.found and r.type_match)
    size_ok = sum(1 for r in results if r.found and r.size_match)
    gt_ok = sum(1 for r in results if r.found and r.genotype_match)
    all_ok = sum(
        1
        for r in results
        if r.found and r.type_match and r.size_match and r.genotype_match
    )

    statistics = {
        "total": total,
        "found": found,
        "type_ok": type_ok,
        "size_ok": size_ok,
        "gt_ok": gt_ok,
        "all_ok": all_ok,
    }

    return results, statistics


def check_vcf_correctness(
    vcf_path: Path, expected_path: Path, output_path: Path | None = None
) -> None:
    """Check VCF correctness against expected variants."""
    log.info(f"Reading expected variants from {expected_path}")
    results, statistics = compare_vcf_to_expected(vcf_path, expected_path)

    total = statistics["total"]
    found = statistics["found"]
    type_ok = statistics["type_ok"]
    size_ok = statistics["size_ok"]
    gt_ok = statistics["gt_ok"]
    all_ok = statistics["all_ok"]

    # Output results
    output_file = open(output_path, "w") if output_path else sys.stdout

    try:
        # Write header
        output_file.write("REGION\tEXP_TYPE\tEXP_SIZE\tEXP_GT\tSTATUS\tDETAILS\n")

        # Write results
        for result in results:
            output_file.write(str(result) + "\n")
            # Print incorrect variants to terminal
            if result.found and not (
                result.type_match and result.size_match and result.genotype_match
            ):
                log.warning(f"Incorrect variant: {result}")
            elif not result.found:
                log.warning(f"Missing variant: {result}")

        output_file.write("\n# SUMMARY\n")
        output_file.write(f"# Total expected variants: {total}\n")
        output_file.write(f"# Found: {found} ({found / total * 100:.1f}%)\n")
        output_file.write(f"# Type match: {type_ok} ({type_ok / total * 100:.1f}%)\n")
        output_file.write(f"# Size match: {size_ok} ({size_ok / total * 100:.1f}%)\n")
        output_file.write(f"# Genotype match: {gt_ok} ({gt_ok / total * 100:.1f}%)\n")
        output_file.write(f"# All correct: {all_ok} ({all_ok / total * 100:.1f}%)\n")

        log.info(f"Found {found}/{total} expected variants")
        log.info(f"All correct: {all_ok}/{total}")

    finally:
        if output_path:
            output_file.close()


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Check VCF correctness against expected variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--vcf",
        type=Path,
        required=True,
        help="Path to VCF file (bgzipped and tabix-indexed)",
    )
    parser.add_argument(
        "-e",
        "--expected",
        type=Path,
        required=True,
        help="Path to expected variants file (region:start-end SVTYPE size HET/HOM)",
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

    check_vcf_correctness(
        vcf_path=args.vcf, expected_path=args.expected, output_path=args.output
    )


if __name__ == "__main__":
    main()
