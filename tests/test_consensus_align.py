# %%
import logging
import tempfile
from pathlib import Path

import pytest

from svirlpool.localassembly import (
    SVpatterns,
    SVprimitives,
    consensus_align,
    consensus_class,
)
from svirlpool.util import datatypes

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
# F7: the "no distortions" guard in
# add_consensus_sequence_and_size_distortions_to_svPatterns must inspect the
# *values* of the distortion dict, not only its length.
#
# distortions_by_svPattern never returns an empty dict for a pattern that has
# supporting reads: when it finds no distortion it can weight, it returns
# dict.fromkeys(supporting_reads, 0.0). The length-only guard cannot see that.
#
# The fixture below is a genuinely matched (Consensus, SVprimitives) pair taken
# from ONE real run, so the objects under test are objects the pipeline can
# actually produce. It was generated with:
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
# ---------------------------------------------------------------------------

MATCHED_PAIR_FIXTURE = (
    DATADIR / "consensus_align" / "all_zero_size_distortions.chr6_157299125.json.gz"
)


@pytest.fixture(autouse=True)
def _reset_all_zero_distortion_warning_counter():
    """Reset the per-process warning rate limiter between tests."""
    getattr(consensus_align, "_reset_all_zero_distortion_warnings", lambda: None)()
    yield


def load_matched_pattern_and_consensus() -> tuple[
    SVpatterns.SVpatternType, consensus_class.Consensus
]:
    """Load a genuinely matched SVpattern/Consensus pair from one real run.

    The SVpattern is rebuilt from its SVprimitives the way production does, so it
    arrives with ``size_distortions is None`` — the state
    ``add_consensus_sequence_and_size_distortions_to_svPatterns`` expects.

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


def attach_distortions_to_consensus(
    consensus: consensus_class.Consensus,
    distortions: dict[str, list[tuple[int, int, int]]],
) -> consensus_class.Consensus:
    """Attach cut-read alignment signals to a Consensus.

    ``distortions`` maps a read name to a list of ``(position, size, sv_type)``
    triples in **core-consensus** coordinates — the space in which
    ``Consensus.get_consensus_distortions`` surfaces ``signal.ref_start``.
    """
    core_length = len(consensus.consensus_sequence)
    for signals in distortions.values():
        for position, _size, _sv_type in signals:
            assert 0 <= position < core_length, (
                f"distortion position {position} is outside the consensus "
                f"[0, {core_length})"
            )

    consensus.cut_read_alignment_signals = [
        datatypes.ReadAlignmentSignals(
            samplename="testsample",
            read_name=readname,
            reference_name=consensus.ID,
            alignment_forward=True,
            SV_signals=[
                datatypes.SVsignal(
                    ref_start=position,
                    ref_end=position + max(size, 1),
                    read_start=position,
                    read_end=position + max(size, 1),
                    size=size,
                    sv_type=sv_type,
                )
                for position, size, sv_type in signals
            ],
        )
        for readname, signals in distortions.items()
    ]
    return consensus


def _spread_distortions(
    readnames, core_length: int
) -> dict[str, list[tuple[int, int, int]]]:
    """One insertion distortion per read, spread across the core consensus."""
    readnames = list(readnames)
    step = max(core_length // (len(readnames) + 2), 1)
    return {
        readname: [(step * (i + 1), 30 + i, 0)] for i, readname in enumerate(readnames)
    }


def test_matched_pair_fixture_is_a_coherent_production_object():
    """The fixture really is one pattern and one consensus from the same run."""
    pattern, consensus = load_matched_pattern_and_consensus()

    assert "0.1" == consensus.ID
    assert ("chr6", 157299125, 157299125) == pattern.get_reference_region()
    assert isinstance(pattern, SVpatterns.SVpatternInsertion)
    assert pattern.size_distortions is None

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


def test_real_locus_without_read_distortions_yields_an_all_zero_dict():
    """The state the F7 guard has to be able to see, straight from real data.

    At this locus no cut read carries an indel signal, so
    ``distortions_by_svPattern`` returns ``dict.fromkeys(supporting_reads, 0.0)``:
    one entry per supporting read, every value exactly 0.0. The dict is not
    empty, so the old ``len(...) > 0`` guard passes.

    This route to an all-zero dict does not depend on the distance weighting at
    all — there is nothing to weight — so it holds both before and after F1.
    """
    pattern, consensus = load_matched_pattern_and_consensus()
    supporting_reads = pattern.get_supporting_reads()
    assert [] == [
        d
        for d in consensus.get_consensus_distortions()
        if d.readname in set(supporting_reads)
    ]

    result = SVpatterns.distortions_by_svPattern(
        svPattern=pattern,
        consensus=consensus,
        distance_scale=5000.0,
        falloff=1.0,
    )
    assert len(supporting_reads) == len(result)
    assert all(value == 0.0 for value in result.values())


def test_guard_detects_all_zero_size_distortions(caplog):
    """All-zero distortions must be detected, counted and warned about.

    On unmodified ``main`` the guard only checks ``len(...) > 0``, so this
    pathological locus passes silently and nothing at all is logged. The
    warning assertion therefore comes first: it is the substantive failure.
    """
    pattern, consensus = load_matched_pattern_and_consensus()
    chrom, start, _end = pattern.get_reference_region()

    with caplog.at_level(logging.WARNING, logger=consensus_align.log.name):
        result = (
            consensus_align.add_consensus_sequence_and_size_distortions_to_svPatterns(
                consensus_objects={consensus.ID: consensus},
                svPatterns=[pattern],
                distance_scale=5000.0,
                falloff=1.0,
            )
        )

    messages = [rec.getMessage() for rec in caplog.records]
    assert any(consensus.ID in msg and f"{chrom}:{start}" in msg for msg in messages), (
        "an all-zero size-distortion dict must be reported with the consensus ID "
        f"and the genomic coordinate; logged instead: {messages}"
    )

    processed, stats = result
    # observability only: the pattern is still processed and still carries the
    # (all-zero) distortions it had before.
    assert 1 == len(processed)
    assert processed[0].size_distortions is not None
    assert all(v == 0.0 for v in processed[0].size_distortions.values())

    assert 1 == stats["n_all_zero_distortions"]
    assert 1 == stats["n_patterns"]


def test_guard_still_raises_on_a_genuinely_empty_distortion_dict():
    """An empty dict means the pattern has no supporting reads at all.

    That is a real invariant violation and stays a hard failure.
    """
    pattern, consensus = load_matched_pattern_and_consensus()
    for svp in pattern.SVprimitives:
        svp.genotypeMeasurement.supporting_reads_start = []
        svp.genotypeMeasurement.supporting_reads_end = []
    assert [] == pattern.get_supporting_reads()

    with pytest.raises(ValueError, match="No size distortions"):
        consensus_align.add_consensus_sequence_and_size_distortions_to_svPatterns(
            consensus_objects={consensus.ID: consensus},
            svPatterns=[pattern],
            distance_scale=5000.0,
            falloff=1.0,
        )


def test_guard_is_silent_for_a_healthy_all_nonzero_distortion_dict(caplog):
    """A dict with only non-zero values is normal: no warning, no error."""
    pattern, consensus = load_matched_pattern_and_consensus()
    supporting_reads = pattern.get_supporting_reads()
    consensus = attach_distortions_to_consensus(
        consensus,
        _spread_distortions(supporting_reads, len(consensus.consensus_sequence)),
    )

    # A distance scale large enough that the weight cannot underflow whichever
    # coordinate space the distance is measured in, so this test says nothing
    # about F1 either way.
    with caplog.at_level(logging.WARNING, logger=consensus_align.log.name):
        processed, stats = (
            consensus_align.add_consensus_sequence_and_size_distortions_to_svPatterns(
                consensus_objects={consensus.ID: consensus},
                svPatterns=[pattern],
                distance_scale=1e9,
                falloff=1.0,
            )
        )

    assert 1 == len(processed)
    distortions = processed[0].size_distortions
    assert distortions is not None
    assert len(distortions) == len(supporting_reads)
    assert all(value != 0.0 for value in distortions.values())
    assert 0 == stats["n_all_zero_distortions"]
    assert not [
        rec for rec in caplog.records if "ALL_ZERO_DISTORTIONS" in rec.getMessage()
    ]


def test_guard_treats_a_partly_zero_distortion_dict_as_normal(caplog):
    """Only an *all*-zero dict is pathological.

    Reads without a distortion legitimately contribute 0.0 — a locus where some
    reads are distorted and others are not is exactly what the noise model is
    meant to measure, so it must not be flagged.
    """
    pattern, consensus = load_matched_pattern_and_consensus()
    supporting_reads = sorted(pattern.get_supporting_reads())
    assert len(supporting_reads) >= 4
    distorted = supporting_reads[: len(supporting_reads) // 2]
    consensus = attach_distortions_to_consensus(
        consensus, _spread_distortions(distorted, len(consensus.consensus_sequence))
    )

    with caplog.at_level(logging.WARNING, logger=consensus_align.log.name):
        processed, stats = (
            consensus_align.add_consensus_sequence_and_size_distortions_to_svPatterns(
                consensus_objects={consensus.ID: consensus},
                svPatterns=[pattern],
                distance_scale=1e9,
                falloff=1.0,
            )
        )

    distortions = processed[0].size_distortions
    assert distortions is not None
    zeros = [v for v in distortions.values() if v == 0.0]
    nonzeros = [v for v in distortions.values() if v != 0.0]
    assert len(zeros) > 0 and len(nonzeros) > 0, "test setup must be partly zero"
    assert 0 == stats["n_all_zero_distortions"]
    assert not [
        rec for rec in caplog.records if "ALL_ZERO_DISTORTIONS" in rec.getMessage()
    ]


def test_all_zero_distortion_warnings_are_rate_limited(caplog):
    """The warning is rate limited; the *count* is complete.

    While F1 is unfixed essentially every locus in the genome is all-zero, so an
    unlimited per-locus warning would emit millions of lines. Only the first few
    loci are warned about, but all of them are counted for the run summary.
    """
    import copy

    pattern, consensus = load_matched_pattern_and_consensus()
    n_loci = consensus_align.ALL_ZERO_DISTORTION_WARN_LIMIT + 7
    patterns = [copy.deepcopy(pattern) for _ in range(n_loci)]

    with caplog.at_level(logging.WARNING, logger=consensus_align.log.name):
        _processed, stats = (
            consensus_align.add_consensus_sequence_and_size_distortions_to_svPatterns(
                consensus_objects={consensus.ID: consensus},
                svPatterns=patterns,
                distance_scale=5000.0,
                falloff=1.0,
            )
        )

    assert n_loci == stats["n_all_zero_distortions"], "every locus must be counted"
    emitted = [
        rec for rec in caplog.records if "ALL_ZERO_DISTORTIONS" in rec.getMessage()
    ]
    assert len(emitted) <= consensus_align.ALL_ZERO_DISTORTION_WARN_LIMIT + 1, (
        f"expected at most the warn limit (+1 suppression notice), got {len(emitted)}"
    )
