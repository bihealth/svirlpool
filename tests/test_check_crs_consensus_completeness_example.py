"""
Example test for check_crs_consensus_completeness module.

This demonstrates how to use compare_crs_to_consensus() in unit tests.
"""

from pathlib import Path

from svirlpool.debugtools.check_crs_consensus_completeness import (
    compare_crs_to_consensus,
)


def test_crs_consensus_completeness():
    """Test that all candidate regions have consensus sequences."""
    crs_bed_path = Path("path/to/crs.bed")
    consensus_fasta_path = Path("path/to/consensus.fasta.gz")

    # Get detailed results
    results, statistics = compare_crs_to_consensus(crs_bed_path, consensus_fasta_path)

    # Collect CRs without consensus
    missing_crs = [r for r in results if not r.found]

    # Build informative error message
    error_msg_parts = []

    if missing_crs:
        error_msg_parts.append(
            f"\nMissing consensus for {len(missing_crs)} candidate regions:"
        )
        for m in missing_crs:
            error_msg_parts.append(f"  - CR {m.cr_id}: {m.cr_region}")

    # Assert with detailed message
    assert statistics["missing"] == 0, (
        f"Expected all {statistics['total']} candidate regions to have consensus, "
        f"but {statistics['missing']} are missing. " + "\n".join(error_msg_parts)
    )

    # Alternative: check specific thresholds
    assert statistics["found"] >= statistics["total"] * 0.95, (
        f"Expected at least 95% CRs with consensus, "
        f"got {statistics['found']}/{statistics['total']}"
    )

    # Check average consensus per CR
    avg_consensus = statistics["total_consensus"] / statistics["total"]
    assert avg_consensus >= 1.0, (
        f"Expected at least 1 consensus per CR on average, got {avg_consensus:.2f}"
    )


def test_specific_crs_have_consensus():
    """Test that specific important candidate regions have consensus sequences."""
    crs_bed_path = Path("path/to/crs.bed")
    consensus_fasta_path = Path("path/to/consensus.fasta.gz")

    # Get detailed results
    results, statistics = compare_crs_to_consensus(crs_bed_path, consensus_fasta_path)

    # Create lookup dictionary
    results_by_id = {r.cr_id: r for r in results}

    # Check specific CRs that are known to be important
    important_cr_ids = [10, 42, 123]  # Example IDs

    for cr_id in important_cr_ids:
        assert cr_id in results_by_id, f"CR {cr_id} not found in database"
        result = results_by_id[cr_id]
        assert result.found, (
            f"Important CR {cr_id} ({result.cr_region}) has no consensus sequences"
        )
        assert result.consensus_count > 0, (
            f"Important CR {cr_id} ({result.cr_region}) has 0 consensus sequences"
        )


def test_no_crs_with_excessive_consensus():
    """Test that no candidate region has an excessive number of consensus sequences."""
    crs_bed_path = Path("path/to/crs.bed")
    consensus_fasta_path = Path("path/to/consensus.fasta.gz")

    # Get detailed results
    results, statistics = compare_crs_to_consensus(crs_bed_path, consensus_fasta_path)

    max_consensus_per_cr = 10  # Reasonable threshold

    excessive_crs = [r for r in results if r.consensus_count > max_consensus_per_cr]

    assert len(excessive_crs) == 0, (
        f"Found {len(excessive_crs)} CRs with more than {max_consensus_per_cr} consensus sequences: "
        + ", ".join(f"CR {r.cr_id} ({r.consensus_count})" for r in excessive_crs[:5])
    )
