# %%
import logging
import tempfile
from pathlib import Path

from svirlpool.localassembly import (
    SVpatterns,
    SVprimitives,
    consensus_align,
    consensus_class,
)

# %%

DATADIR = Path(__file__).parent / "data"

# %%


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
                (0, "chr1", 50, 150, 0, 0),  # Overlaps TRF at 100-200
                (1, "chr1", 1000, 2000, 0, 0),  # Overlaps TRF at 1500-1600
            ],
            "consensus2": [
                (0, "chr1", 450, 650, 0, 0)  # Overlaps TRF at 500-600
            ],
            "consensus3": [
                (0, "chr1", 3000, 4000, 0, 0)  # No overlaps
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
                (0, "chr1", 50, 600, 0, 0)  # Overlaps all three TRFs
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
        core_intervals = {
            "consensus1": [(0, "chr1", 100, 200, 0, 0), (1, "chr1", 300, 400, 0, 0)]
        }

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
                (0, "chr1", 50, 150, 0, 0),  # Overlaps chr1 TRF
                (1, "chr2", 450, 650, 0, 0),  # Overlaps chr2 TRF
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


# %%

# data was saved like this:
# # DEBUG START
# # if the consensusID of one of the SVprimitives is 15.0 then save the parameters to json to debug and test
# if any(svp.consensusID == "15.0" for svp in SVprimitives):
#     import json
#     from gzip import open as gzip_open
#     data = {
#         "SVprimitives": [svp.unstructure() for svp in SVprimitives],
#         "max_del_size": max_del_size,
#     }
#     with gzip_open("/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_align/svPrimitives_to_svPatterns.INV15.json.gz", "wt") as debug_f:
#         json.dump(data, debug_f, indent=4)
# # DEBUG END


def load_svprimitives(path: Path) -> list[SVprimitives.SVprimitive]:
    """Load SVprimitives from a gzipped json file for testing."""
    import json
    from gzip import open as gzip_open

    import cattrs

    with gzip_open(path, "rt") as f:
        data = json.load(f)
    sv_primitives = [
        cattrs.structure(svp_dict, SVprimitives.SVprimitive)
        for svp_dict in data["SVprimitives"]
    ]
    return sv_primitives


def load_max_del_size(path: Path) -> int:
    """Load max_del_size from a gzipped json file for testing."""
    import json
    from gzip import open as gzip_open

    with gzip_open(path, "rt") as f:
        data = json.load(f)
    return data["max_del_size"]


def test_svPrimitives_to_svPatterns_inv15():
    """Test svPrimitives_to_svPatterns function with SVprimitives for INV15."""
    # Load SVprimitives and max_del_size for INV15
    svp_path = DATADIR / "consensus_align" / "svPrimitives_to_svPatterns.INV15.json.gz"
    sv_primitives = load_svprimitives(svp_path)
    max_del_size = load_max_del_size(svp_path)

    group = [svp for svp in sv_primitives if svp.consensusID == "15.0"]

    result = SVpatterns.parse_SVprimitives_to_SVpatterns(
        SVprimitives=group, max_del_size=max_del_size, log_level_override=logging.DEBUG
    )

    assert 1 == len(result)
    assert type(result[0]) == SVpatterns.SVpatternInversion
    assert 4 == len(result[0].SVprimitives)
    assert 13 == result[0].SVprimitives[0].get_total_support()


# %%
# ---------------------------------------------------------------------------
# add_consensus_sequences_to_svPatterns on a genuinely matched
# (Consensus, SVprimitives) pair taken from ONE real run, so the objects under
# test are objects the pipeline can actually produce. It was generated with:
#
#   from svirlpool.localassembly import consensus_align, SVpatterns
#   cons = {}
#   for ccr in consensus_align.parse_crs_container_results(
#           RUN / "wd/consensus/0/consensus.batch_0.jsonl"):
#       cons.update(ccr.consensus_dicts)
#   c = cons["0.1"]
#   p = [p for p in SVpatterns.read_svPatterns_from_db(RUN / "wd/svirltile.db")
#        if p.consensusID == "0.1" and p.ref_start == 157299125][0]
#   json.dump({"consensus": c.unstructure(),
#              "SVprimitives": [s.unstructure() for s in p.SVprimitives],
#              "max_del_size": 100_000}, gzip.open(OUT, "wt"), indent=1)
#
# where RUN is a svirlpool run on HG002 chr6:157,290,937-157,340,937 (GRCh38,
# production defaults). Do NOT pair fixtures by consensus ID alone: IDs are
# "<crID>.<subID>" and are only unique within a run, so two files from
# different datasets can share an ID while describing loci megabases apart.
#
# The consensus record predates the removal of the per-read size distortions
# and still carries ``cut_read_alignment_signals``; loading it is also the
# backward-compatibility check for old consensus JSONL.
# ---------------------------------------------------------------------------

MATCHED_PAIR_FIXTURE = (
    DATADIR / "consensus_align" / "matched_pair.chr6_157299125.json.gz"
)


def load_matched_pattern_and_consensus() -> tuple[
    SVpatterns.SVpatternType, consensus_class.Consensus
]:
    """Load a genuinely matched SVpattern/Consensus pair from one real run.

    The SVpattern is rebuilt from its SVprimitives the way production does, so it
    arrives without any consensus-derived sequence — the state
    ``add_consensus_sequences_to_svPatterns`` expects.

    The assertions below are the point of this helper: they pin the invariants a
    matched pair must satisfy, so a mismatched pair can never silently be used as
    test input again.
    """
    import json
    from gzip import open as gzip_open

    import cattrs

    with gzip_open(MATCHED_PAIR_FIXTURE, "rt") as f:
        data = json.load(f)

    consensus = cattrs.structure(data["consensus"], consensus_class.Consensus)
    sv_primitives = [
        cattrs.structure(svp, SVprimitives.SVprimitive) for svp in data["SVprimitives"]
    ]
    patterns = SVpatterns.parse_SVprimitives_to_SVpatterns(
        SVprimitives=sv_primitives, max_del_size=data["max_del_size"]
    )
    assert 1 == len(patterns), (
        f"fixture must describe exactly one SVpattern: {patterns}"
    )
    pattern = patterns[0]

    # --- invariants of a matched pair -------------------------------------
    assert consensus.ID == pattern.consensusID, (
        f"consensus {consensus.ID} does not belong to pattern {pattern.consensusID}"
    )

    supporting_reads = set(pattern.get_supporting_reads())
    consensus_reads = consensus.get_used_readnames()
    assert supporting_reads, "the pattern must have supporting reads"
    assert supporting_reads <= consensus_reads, (
        "every supporting read of the pattern must be a read of the consensus; "
        f"{len(supporting_reads - consensus_reads)} of {len(supporting_reads)} are not — "
        "the two fixtures do not belong together"
    )

    assert consensus.consensus_padding is not None
    padding_left = consensus.consensus_padding.padding_size_left
    core_start = pattern.read_start - padding_left
    core_end = pattern.read_end - padding_left
    assert 0 <= core_start <= core_end <= len(consensus.consensus_sequence), (
        f"the pattern's consensus-local interval [{core_start}, {core_end}] falls "
        f"outside the consensus sequence [0, {len(consensus.consensus_sequence)}] — "
        "the two fixtures do not belong together"
    )

    return pattern, consensus


def test_matched_pair_fixture_is_a_coherent_production_object():
    """The fixture really is one pattern and one consensus from the same run."""
    pattern, consensus = load_matched_pattern_and_consensus()

    assert "0.1" == consensus.ID
    assert ("chr6", 157299125, 157299125) == pattern.get_reference_region()
    assert isinstance(pattern, SVpatterns.SVpatternInsertion)
    assert pattern.inserted_sequence is None

    # the invariants load_matched_pattern_and_consensus() asserts, restated here
    # so that a regression names the broken invariant rather than a helper.
    supporting_reads = set(pattern.get_supporting_reads())
    assert 19 == len(supporting_reads)
    assert supporting_reads <= consensus.get_used_readnames()
    padding_left = consensus.consensus_padding.padding_size_left
    assert (250, 280) == (
        pattern.read_start - padding_left,
        pattern.read_end - padding_left,
    )
    assert 530 == len(consensus.consensus_sequence)


def test_an_old_consensus_record_with_read_signals_still_loads():
    """Consensus JSON written before the size distortions were removed loads.

    The dropped ``cut_read_alignment_signals`` key is ignored; everything else,
    in particular the cut-read intervals, is kept.
    """
    import json
    from gzip import open as gzip_open

    import cattrs

    with gzip_open(MATCHED_PAIR_FIXTURE, "rt") as f:
        data = json.load(f)["consensus"]
    assert data["cut_read_alignment_signals"], "fixture must carry the old field"

    consensus = cattrs.structure(data, consensus_class.Consensus)
    assert not hasattr(consensus, "cut_read_alignment_signals")
    assert len(data["intervals_cutread_alignments"]) == len(
        consensus.intervals_cutread_alignments
    )
    assert "cut_read_alignment_signals" not in consensus.unstructure()


def test_add_consensus_sequences_sets_the_inserted_sequence():
    """The insertion receives its inserted sequence from the consensus."""
    pattern, consensus = load_matched_pattern_and_consensus()
    expected = pattern.get_sequence_from_consensus(consensus=consensus)
    assert 30 == len(expected)

    processed = consensus_align.add_consensus_sequences_to_svPatterns(
        consensus_objects={consensus.ID: consensus},
        svPatterns=[pattern],
    )

    assert [pattern] == processed
    assert expected == processed[0].get_sequence()
    assert processed[0].get_sequence_complexity() is not None


def test_add_consensus_sequences_keeps_patterns_without_a_consensus():
    """A pattern whose consensus is not in the batch is passed through unchanged."""
    pattern, _consensus = load_matched_pattern_and_consensus()

    processed = consensus_align.add_consensus_sequences_to_svPatterns(
        consensus_objects={},
        svPatterns=[pattern],
    )

    assert [pattern] == processed
    assert processed[0].get_sequence() is None
