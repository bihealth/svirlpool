#%%
import tempfile
from pathlib import Path

from svirlpool.localassembly import SVprimitives, consensus_align
from svirlpool.localassembly.consensus_align_lib import \
    parse_sv_signals_from_consensus
from svirlpool.util import datatypes

#%%

DATADIR = Path(__file__).parent / "data" / "consensus_align"

#%%

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

#%%

def load_parse_sv_signals_from_consensus_test_data(path: Path | str) -> dict[str,object]:
    """Loads test data for parse_sv_signals_from_consensus from a json file."""
    if isinstance(path, str):
        path = Path(path)
    import gzip
    import json

    import cattrs

    with gzip.open(path, "rt", encoding="utf-8") as f:
        input_data = json.load(f)
    # reconstruct pysam alignment from unstructured datatypes.Alignment
    pysam_alignment = cattrs.structure(input_data["pysam_alignment"], datatypes.Alignment).to_pysam()
    return {
        "samplename": input_data["samplename"],
        "reference_name": input_data["reference_name"],
        "consensus_sequence": input_data["consensus_sequence"],
        "pysam_alignment": pysam_alignment,
        "interval_core": tuple(input_data["interval_core"]),
        "trf_intervals": [tuple(trf) for trf in input_data["trf_intervals"]],
        "min_signal_size": input_data["min_signal_size"],
        "min_bnd_size": input_data["min_bnd_size"]}

def test_generate_SVprimitives_inv15():
    """Test parse_sv_signals_from_consensus with inversion and reverse complement cases."""
    # three alignments are part of the INV 15 locus:
    # parse_sv_signals_from_consensus_INV15.168996671.json.gz
    # parse_sv_signals_from_consensus_INV15.169093632.json.gz
    # parse_sv_signals_from_consensus_INV15.169094647.json.gz
    input_files = [
        DATADIR / "parse_sv_signals_from_consensus_INV15.168996671.json.gz",
        DATADIR / "parse_sv_signals_from_consensus_INV15.169093632.json.gz",
        DATADIR / "parse_sv_signals_from_consensus_INV15.169094647.json.gz",
    ]
    # parse all of them to 
    inputs = [load_parse_sv_signals_from_consensus_test_data(f) for f in input_files]

    results = []
    for test_data in inputs:
        merged_svs = parse_sv_signals_from_consensus(
            samplename=test_data["samplename"],
            reference_name=test_data["reference_name"],
            consensus_sequence=test_data["consensus_sequence"],
            pysam_alignment=test_data["pysam_alignment"],
            interval_core=test_data["interval_core"],
            trf_intervals=test_data["trf_intervals"],
            min_signal_size=test_data["min_signal_size"],
            min_bnd_size=test_data["min_bnd_size"],
        )
        for sv in merged_svs:
            results.append(sv)
            
    # test consensus_align.generate_SVprimitives
    svPrimitives = SVprimitives.generate_SVprimitives(results)
    
    
    svPrimitives.extend(
        SVprimitives.generate_SVprimitives(
            samplename=params.samplename,
            mergedSVs=mergedSVs,
            consensus=consensus_obj,
            consensus_alignment=alignment,
            alignmentID=alignment_idx,
            core_interval=core_interval_by_idx[alignment_idx],
        )

    svPatterns = consensus_align.svPrimitives_to_svPatterns(
        SVprimitives=svPrimitives, max_del_size=params.max_del_size
    )

# %%
