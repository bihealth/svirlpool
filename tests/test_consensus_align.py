import tempfile
from pathlib import Path

from svirlpool.localassembly import consensus_align


def test_process_partition_for_trf_overlaps_basic():
    """Test basic TRF overlap detection with simple overlapping intervals."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Create a dummy reference fasta and index
        ref_fasta = tmp_path / "reference.fa"
        ref_fai = tmp_path / "reference.fa.fai"

        with open(ref_fasta, "w") as f:
            f.write(">chr1\n")
            f.write("A" * 10000 + "\n")

        with open(ref_fai, "w") as f:
            f.write("chr1\t10000\t6\t10000\t10001\n")

        # Create a TRF bed file with some repeats
        trf_bed = tmp_path / "trf.bed"
        with open(trf_bed, "w") as f:
            f.write("chr1\t100\t200\n")  # TRF at 100-200
            f.write("chr1\t500\t600\n")  # TRF at 500-600
            f.write("chr1\t1500\t1600\n")  # TRF at 1500-1600

        # Core intervals: consensusID mapped to list of (alignment_idx, chrom, start, end)
        core_intervals = {
            "consensus1": [
                (0, "chr1", 50, 150),  # Overlaps TRF at 100-200
                (1, "chr1", 1000, 2000),  # Overlaps TRF at 1500-1600
            ],
            "consensus2": [
                (0, "chr1", 450, 650)  # Overlaps TRF at 500-600
            ],
            "consensus3": [
                (0, "chr1", 3000, 4000)  # No overlaps
            ],
        }

        # Run the function
        result = consensus_align._process_partition_for_trf_overlaps(
            partition_idx=0,
            core_intervals=core_intervals,
            input_trf=trf_bed,
            reference=ref_fasta,
            tmp_dir=tmp_path,
        )

        # Verify results
        assert "consensus1" in result
        assert "consensus2" in result
        assert "consensus3" in result

        # consensus1 should have TRF overlaps from both alignments (deduplicated)
        assert len(result["consensus1"]) == 2
        assert ("chr1", 100, 200, 0) in result["consensus1"]  # First TRF (line 0)
        assert ("chr1", 1500, 1600, 2) in result["consensus1"]  # Third TRF (line 2)

        # consensus2 should have 1 TRF overlap
        assert len(result["consensus2"]) == 1
        assert ("chr1", 500, 600, 1) in result["consensus2"]  # Second TRF (line 1)

        # consensus3 should have no overlaps
        assert len(result["consensus3"]) == 0


def test_process_partition_for_trf_overlaps_multiple_overlaps():
    """Test TRF overlap detection when one core interval overlaps multiple TRFs."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Create a dummy reference fasta and index
        ref_fasta = tmp_path / "reference.fa"
        ref_fai = tmp_path / "reference.fa.fai"

        with open(ref_fasta, "w") as f:
            f.write(">chr1\n")
            f.write("A" * 10000 + "\n")

        with open(ref_fai, "w") as f:
            f.write("chr1\t10000\t6\t10000\t10001\n")

        # Create a TRF bed file with multiple overlapping repeats
        trf_bed = tmp_path / "trf.bed"
        with open(trf_bed, "w") as f:
            f.write("chr1\t100\t200\n")
            f.write("chr1\t250\t350\n")
            f.write("chr1\t400\t500\n")

        # Core interval that overlaps all three TRFs
        core_intervals = {
            "consensus1": [
                (0, "chr1", 50, 600)  # Overlaps all three TRFs
            ]
        }

        # Run the function
        result = consensus_align._process_partition_for_trf_overlaps(
            partition_idx=0,
            core_intervals=core_intervals,
            input_trf=trf_bed,
            reference=ref_fasta,
            tmp_dir=tmp_path,
        )

        # Verify results
        assert "consensus1" in result
        assert len(result["consensus1"]) == 3

        # Check all three TRF overlaps are present with repeat_ids (line numbers)
        assert ("chr1", 100, 200, 0) in result["consensus1"]
        assert ("chr1", 250, 350, 1) in result["consensus1"]
        assert ("chr1", 400, 500, 2) in result["consensus1"]


def test_process_partition_for_trf_overlaps_no_overlaps():
    """Test TRF overlap detection when there are no overlaps."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Create a dummy reference fasta and index
        ref_fasta = tmp_path / "reference.fa"
        ref_fai = tmp_path / "reference.fa.fai"

        with open(ref_fasta, "w") as f:
            f.write(">chr1\n")
            f.write("A" * 10000 + "\n")

        with open(ref_fai, "w") as f:
            f.write("chr1\t10000\t6\t10000\t10001\n")

        # Create a TRF bed file
        trf_bed = tmp_path / "trf.bed"
        with open(trf_bed, "w") as f:
            f.write("chr1\t1000\t1100\n")
            f.write("chr1\t2000\t2100\n")

        # Core intervals that don't overlap any TRFs
        core_intervals = {"consensus1": [(0, "chr1", 100, 200), (1, "chr1", 300, 400)]}

        # Run the function
        result = consensus_align._process_partition_for_trf_overlaps(
            partition_idx=0,
            core_intervals=core_intervals,
            input_trf=trf_bed,
            reference=ref_fasta,
            tmp_dir=tmp_path,
        )

        # Verify results - should have empty list
        assert "consensus1" in result
        assert len(result["consensus1"]) == 0


def test_process_partition_for_trf_overlaps_multiple_chromosomes():
    """Test TRF overlap detection with multiple chromosomes."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Create a dummy reference fasta and index with multiple chromosomes
        ref_fasta = tmp_path / "reference.fa"
        ref_fai = tmp_path / "reference.fa.fai"

        with open(ref_fasta, "w") as f:
            f.write(">chr1\n")
            f.write("A" * 10000 + "\n")
            f.write(">chr2\n")
            f.write("A" * 10000 + "\n")

        with open(ref_fai, "w") as f:
            f.write("chr1\t10000\t6\t10000\t10001\n")
            f.write("chr2\t10000\t10012\t10000\t10001\n")

        # Create a TRF bed file with repeats on both chromosomes
        trf_bed = tmp_path / "trf.bed"
        with open(trf_bed, "w") as f:
            f.write("chr1\t100\t200\n")
            f.write("chr2\t500\t600\n")

        # Core intervals on different chromosomes
        core_intervals = {
            "consensus1": [
                (0, "chr1", 50, 150),  # Overlaps chr1 TRF
                (1, "chr2", 450, 650),  # Overlaps chr2 TRF
            ]
        }

        # Run the function
        result = consensus_align._process_partition_for_trf_overlaps(
            partition_idx=0,
            core_intervals=core_intervals,
            input_trf=trf_bed,
            reference=ref_fasta,
            tmp_dir=tmp_path,
        )

        # Verify results
        assert "consensus1" in result
        assert len(result["consensus1"]) == 2

        # Should have TRFs from both chromosomes with repeat_ids
        assert ("chr1", 100, 200, 0) in result["consensus1"]
        assert ("chr2", 500, 600, 1) in result["consensus1"]


def test_process_partition_for_trf_overlaps_empty_input():
    """Test TRF overlap detection with empty core intervals."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Create a dummy reference fasta and index
        ref_fasta = tmp_path / "reference.fa"
        ref_fai = tmp_path / "reference.fa.fai"

        with open(ref_fasta, "w") as f:
            f.write(">chr1\n")
            f.write("A" * 10000 + "\n")

        with open(ref_fai, "w") as f:
            f.write("chr1\t10000\t6\t10000\t10001\n")

        # Create a TRF bed file
        trf_bed = tmp_path / "trf.bed"
        with open(trf_bed, "w") as f:
            f.write("chr1\t100\t200\n")

        # Empty core intervals
        core_intervals = {}

        # Run the function
        result = consensus_align._process_partition_for_trf_overlaps(
            partition_idx=0,
            core_intervals=core_intervals,
            input_trf=trf_bed,
            reference=ref_fasta,
            tmp_dir=tmp_path,
        )

        # Verify results - should be empty
        assert result == {}
