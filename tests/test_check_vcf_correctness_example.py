"""
Example test for check_vcf_correctness module.

This demonstrates how to use compare_vcf_to_expected() in unit tests.
"""

from pathlib import Path

from svirlpool.debugtools.check_vcf_correctness import compare_vcf_to_expected


def vcf_correctness():
    """Test that VCF contains all expected variants with correct properties."""
    vcf_path = Path("path/to/your/output.vcf.gz")
    expected_path = Path("path/to/expected_variants.bed")

    # Get detailed results
    results, statistics = compare_vcf_to_expected(vcf_path, expected_path)

    # Collect failures for detailed error messages
    missing = [r for r in results if not r.found]
    incorrect = [
        r
        for r in results
        if r.found and not (r.type_match and r.size_match and r.genotype_match)
    ]

    # Build informative error message
    error_msg_parts = []

    if missing:
        error_msg_parts.append(f"\nMissing {len(missing)} variants:")
        for m in missing:
            error_msg_parts.append(
                f"  - {m.expected.chrom}:{m.expected.start}-{m.expected.end} "
                f"{m.expected.svtype} size={m.expected.size} gt={m.expected.genotype}"
            )

    if incorrect:
        error_msg_parts.append(f"\nIncorrect {len(incorrect)} variants:")
        for inc in incorrect:
            mismatches = []
            if not inc.type_match:
                mismatches.append("TYPE")
            if not inc.size_match:
                mismatches.append("SIZE")
            if not inc.genotype_match:
                mismatches.append("GT")
            error_msg_parts.append(
                f"  - {inc.expected.chrom}:{inc.expected.start}-{inc.expected.end} "
                f"{inc.expected.svtype} size={inc.expected.size} gt={inc.expected.genotype} "
                f"[{','.join(mismatches)} mismatch]"
            )

    # Assert with detailed message
    assert statistics["all_ok"] == statistics["total"], (
        f"Expected all {statistics['total']} variants to be correct, "
        f"but only {statistics['all_ok']} matched perfectly. "
        + "\n".join(error_msg_parts)
    )

    # Alternative: check specific thresholds
    assert statistics["found"] >= statistics["total"] * 0.95, (
        f"Expected at least 95% variants found, got {statistics['found']}/{statistics['total']}"
    )

    assert statistics["all_ok"] >= statistics["total"] * 0.90, (
        f"Expected at least 90% variants correct, got {statistics['all_ok']}/{statistics['total']}"
    )
