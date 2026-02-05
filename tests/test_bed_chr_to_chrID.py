import os
import tempfile
from pathlib import Path

import pytest

from svirlpool.util import util

DATA_DIR = Path(__file__).parent / "data" / "util"


def test_matching_chromosomes():
    """Test with BED file that has all chromosomes present in the reference"""
    reference_fai = DATA_DIR / "ref.fa"
    matching_bed = DATA_DIR / "test_matching.bed"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as tmp:
        output_path = Path(tmp.name)

    try:
        # Run the conversion
        util.bed_chr_to_chrID(
            input=matching_bed, reference=reference_fai, output=output_path, add_id=True
        )

        # Read the output
        with open(output_path, "r") as f:
            lines = [line.strip().split("\t") for line in f]

        # Verify the output
        # chr1 should be ID 0, chr2 -> 1, chr3 -> 2, chrX -> 4
        assert len(lines) == 4, "Should have 4 output lines"

        # Check first line (chr1 -> 0)
        assert lines[0][0] == "0", "chr1 should map to ID 0"
        assert lines[0][1] == "1000", "Start position preserved"
        assert lines[0][2] == "2000", "End position preserved"
        assert lines[0][3] == "region1", "Region name preserved"
        assert lines[0][4] == "100", "Score preserved"
        assert lines[0][5] == "0", "Line ID added"

        # Check second line (chr2 -> 1)
        assert lines[1][0] == "1", "chr2 should map to ID 1"
        assert lines[1][5] == "1", "Line ID should be 1"

        # Check third line (chr3 -> 2)
        assert lines[2][0] == "2", "chr3 should map to ID 2"
        assert lines[2][5] == "2", "Line ID should be 2"

        # Check fourth line (chrX -> 4)
        assert lines[3][0] == "4", "chrX should map to ID 4"
        assert lines[3][5] == "3", "Line ID should be 3"

    finally:
        # Clean up
        if output_path.exists():
            os.unlink(output_path)


def test_mismatching_chromosomes():
    """Test with BED file that has some chromosomes not in the reference"""
    reference_fai = DATA_DIR / "ref.fa"
    mismatching_bed = DATA_DIR / "test_mismatching.bed"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as tmp:
        output_path = Path(tmp.name)

    try:
        # Run the conversion (non-pedantic mode should skip unknown chromosomes)
        util.bed_chr_to_chrID(
            input=mismatching_bed,
            reference=reference_fai,
            output=output_path,
            add_id=True,
            pedantic=False,
        )

        # Read the output
        with open(output_path, "r") as f:
            lines = [line.strip().split("\t") for line in f]

        # Verify the output - should only contain chr1 and chrX (chr5, chr9, chrZ should be skipped)
        assert len(lines) == 2, "Should have 2 output lines (chr5, chr9, chrZ skipped)"

        # Check first line (chr1 -> 0)
        assert lines[0][0] == "0", "chr1 should map to ID 0"
        assert lines[0][3] == "region1", "Region name preserved"
        assert lines[0][5] == "0", "Line ID should be 0 (original line 0)"

        # Check second line (chrX -> 4)
        assert lines[1][0] == "4", "chrX should map to ID 4"
        assert lines[1][3] == "region4", "Region name preserved"
        assert lines[1][5] == "3", "Line ID should be 3 (original line 3)"

    finally:
        # Clean up
        if output_path.exists():
            os.unlink(output_path)


def test_mismatching_chromosomes_pedantic():
    """Test that pedantic mode raises an error for unknown chromosomes"""
    reference_fai = DATA_DIR / "ref.fa"
    mismatching_bed = DATA_DIR / "test_mismatching.bed"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as tmp:
        output_path = Path(tmp.name)

    try:
        # Run the conversion in pedantic mode - should raise KeyError
        with pytest.raises(KeyError, match="Chromosome name.*not found in reference"):
            util.bed_chr_to_chrID(
                input=mismatching_bed,
                reference=reference_fai,
                output=output_path,
                add_id=True,
                pedantic=True,
            )

    finally:
        # Clean up
        if output_path.exists():
            os.unlink(output_path)


def test_without_add_id():
    """Test conversion without adding line IDs"""
    reference_fai = DATA_DIR / "ref.fa"
    matching_bed = DATA_DIR / "test_matching.bed"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as tmp:
        output_path = Path(tmp.name)

    try:
        # Run the conversion without adding IDs
        util.bed_chr_to_chrID(
            input=matching_bed,
            reference=reference_fai,
            output=output_path,
            add_id=False,
        )

        # Read the output
        with open(output_path, "r") as f:
            lines = [line.strip().split("\t") for line in f]

        # Verify the output
        assert len(lines) == 4, "Should have 4 output lines"

        # Check that no extra column was added
        assert len(lines[0]) == 5, "Should have 5 columns (no ID column added)"
        assert lines[0][0] == "0", "chr1 should map to ID 0"
        assert lines[0][4] == "100", "Last column should be score, not ID"

    finally:
        # Clean up
        if output_path.exists():
            os.unlink(output_path)
