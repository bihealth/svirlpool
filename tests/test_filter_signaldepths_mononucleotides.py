import os
import tempfile
from pathlib import Path

from svirlpool.signalprocessing import filter_signaldepths_mononucleotides

DATA_DIR = Path(__file__).parent / "data"


def test_filter_matching_chromosomes():
    """Test with mononucleotides BED file that has all chromosomes present in the reference"""
    reference_fai = DATA_DIR / "util" / "ref.fa"
    matching_bed = DATA_DIR / "util" / "test_matching.bed"
    signaldepths = DATA_DIR / "signalprocessing" / "test_signaldepths_matching.tsv"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv.bgz") as tmp:
        output_path = Path(tmp.name)

    try:
        # Run the filter function
        filter_signaldepths_mononucleotides.filter_mononucleotide_exclusive_deletions(
            mononucleotides=matching_bed,
            signaldepths=signaldepths,
            reference=reference_fai,
            output=output_path,
            margin=5,
            threads=1,
        )

        # Verify the output file was created
        assert output_path.exists(), "Output file should be created"

        # The output should be a bgzipped file
        # We can verify it exists and is not empty
        assert output_path.stat().st_size > 0, "Output file should not be empty"

    finally:
        # Clean up
        if output_path.exists():
            os.unlink(output_path)
        # Also clean up the index file if it exists
        index_path = Path(str(output_path) + ".tbi")
        if index_path.exists():
            os.unlink(index_path)


def test_filter_mismatching_chromosomes():
    """Test with mononucleotides BED file that has some chromosomes not in the reference"""
    reference_fai = DATA_DIR / "util" / "ref.fa"
    mismatching_bed = DATA_DIR / "util" / "test_mismatching.bed"
    signaldepths = DATA_DIR / "signalprocessing" / "test_signaldepths_mismatching.tsv"

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv.bgz") as tmp:
        output_path = Path(tmp.name)

    try:
        # Run the filter function - should skip chromosomes not in reference
        filter_signaldepths_mononucleotides.filter_mononucleotide_exclusive_deletions(
            mononucleotides=mismatching_bed,
            signaldepths=signaldepths,
            reference=reference_fai,
            output=output_path,
            margin=5,
            threads=1,
        )

        # Verify the output file was created
        assert output_path.exists(), "Output file should be created"

        # The output should be a bgzipped file
        assert output_path.stat().st_size > 0, "Output file should not be empty"

    finally:
        # Clean up
        if output_path.exists():
            os.unlink(output_path)
        # Also clean up the index file if it exists
        index_path = Path(str(output_path) + ".tbi")
        if index_path.exists():
            os.unlink(index_path)


def test_filter_with_different_margins():
    """Test that different margin values work correctly"""
    reference_fai = DATA_DIR / "util" / "ref.fa"
    matching_bed = DATA_DIR / "util" / "test_matching.bed"
    signaldepths = DATA_DIR / "signalprocessing" / "test_signaldepths_matching.tsv"

    for margin in [0, 5, 10]:
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".tsv.bgz"
        ) as tmp:
            output_path = Path(tmp.name)

        try:
            # Run the filter function with different margins
            filter_signaldepths_mononucleotides.filter_mononucleotide_exclusive_deletions(
                mononucleotides=matching_bed,
                signaldepths=signaldepths,
                reference=reference_fai,
                output=output_path,
                margin=margin,
                threads=1,
            )

            # Verify the output file was created
            assert output_path.exists(), (
                f"Output file should be created for margin={margin}"
            )
            assert output_path.stat().st_size > 0, (
                f"Output file should not be empty for margin={margin}"
            )

        finally:
            # Clean up
            if output_path.exists():
                os.unlink(output_path)
            # Also clean up the index file if it exists
            index_path = Path(str(output_path) + ".tbi")
            if index_path.exists():
                os.unlink(index_path)
