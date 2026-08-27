"""Regression guards for run-to-run reproducibility of persisted output.

Python randomises string hashing per process (``PYTHONHASHSEED=random`` is the
default), so ``set`` iteration order over ``str`` elements differs between
interpreter processes.  Anywhere such an order is frozen into a list or into a
``dict`` that is later serialised, two runs of the pipeline on identical input
produce byte-different output and can no longer be regression-diffed.

Two complementary styles of test live here:

* **Cross-process tests** re-run a tiny snippet under several explicit
  ``PYTHONHASHSEED`` values and assert the results agree.  These are the ones
  that can actually *observe* the defect -- a single pytest process has one
  fixed hash seed, so an in-process test cannot.
* **Positive-contract tests** assert the stronger, durable property directly:
  the returned sequence equals its own sorted order, and is independent of the
  order in which the inputs were supplied.  These stay meaningful forever,
  including under ``PYTHONHASHSEED=0``.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from svirlpool.localassembly import SVpatterns
from svirlpool.localassembly.SVprimitives import SVprimitive
from svirlpool.svcalling import genotyping, multisample_sv_calling

# Readnames chosen so that their `set` iteration order is genuinely unstable
# across hash seeds (long, high-entropy strings similar to ONT read UUIDs).
READS_START = [
    "b2c5f0e1-1111-4a01-9c3e-000000000001",
    "0f9a7d34-2222-4a02-9c3e-000000000002",
    "7e13ab90-3333-4a03-9c3e-000000000003",
    "c48d6f22-4444-4a04-9c3e-000000000004",
    "1a2b3c4d-5555-4a05-9c3e-000000000005",
]
READS_END = [
    "9d8c7b6a-6666-4a06-9c3e-000000000006",
    "3f3f3f3f-7777-4a07-9c3e-000000000007",
    "b2c5f0e1-1111-4a01-9c3e-000000000001",  # deliberate overlap with START
]

HASH_SEEDS = ("1", "2", "3", "4", "5")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _src_root() -> str:
    """Directory that has to be on ``sys.path`` for ``import svirlpool`` to work.

    Derived from the *imported* package so that a subprocess tests the same
    working copy pytest is testing (worktrees vs. the editable install).
    """
    import svirlpool

    return str(Path(svirlpool.__file__).resolve().parent.parent)


def _run_under_hash_seed(seed: str, body: str) -> str:
    """Run ``body`` in a fresh interpreter with an explicit ``PYTHONHASHSEED``."""
    script = textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {_src_root()!r})
        {textwrap.indent(textwrap.dedent(body), "        ").lstrip()}
    """)
    proc = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env={
            "PYTHONHASHSEED": seed,
            "PATH": "/usr/bin:/bin",
            "HOME": "/tmp",
        },
        timeout=120,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"subprocess with PYTHONHASHSEED={seed} failed:\n{proc.stderr}"
        )
    return proc.stdout.strip()


def _make_svprimitive(
    *,
    reads_start: list[str],
    reads_end: list[str] | None = None,
    ref_start: int = 1000,
    ref_end: int = 1000,
) -> SVprimitive:
    return SVprimitive(
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=100,
        read_end=600,
        size=500,
        sv_type=0,  # INS
        chr="chr1",
        repeatIDs=[],
        original_alt_sequences=["A" * 500],
        original_ref_sequences=[],
        samplename="sample1",
        consensusID="1.0",
        alignmentID=0,
        svID=0,
        aln_is_reverse=False,
        consensus_aln_interval=("chr1", 500, 1500),
        genotypeMeasurement=genotyping.GenotypeMeasurement(
            start_on_consensus=100,
            supporting_reads_start=list(reads_start),
            end_on_consensus=None if reads_end is None else 600,
            supporting_reads_end=None if reads_end is None else list(reads_end),
        ),
    )


def _make_pattern(
    reads_start: list[str], reads_end: list[str] | None = None
) -> SVpatterns.SVpatternInsertion:
    return SVpatterns.SVpatternInsertion(
        SVprimitives=[_make_svprimitive(reads_start=reads_start, reads_end=reads_end)]
    )


# --------------------------------------------------------------------------- #
# SVpattern.get_supporting_reads
# --------------------------------------------------------------------------- #
def test_get_supporting_reads_is_sorted() -> None:
    """The returned sequence must equal its own sorted order."""
    pattern = _make_pattern(READS_START)
    reads = pattern.get_supporting_reads()
    assert reads == sorted(reads)


def test_get_supporting_reads_is_sorted_with_end_reads() -> None:
    """Same contract on the ``start + end`` branch of get_supporting_reads."""
    pattern = _make_pattern(READS_START, READS_END)
    reads = pattern.get_supporting_reads()
    assert reads == sorted(reads)
    # de-duplication must still happen: one readname is in both lists
    assert len(reads) == len(set(reads))
    assert set(reads) == set(READS_START) | set(READS_END)


@pytest.mark.parametrize("with_end", [False, True])
def test_get_supporting_reads_is_independent_of_insertion_order(
    with_end: bool,
) -> None:
    """Supplying the same readnames in different orders must give one result.

    Note on strength: insertion-order independence *alone* is a weak detector,
    because within one process ``set`` iteration over the same small element
    set usually lands in the same order whatever the insertion order -- whether
    it does depends on the hash seed of that particular pytest process.  The
    canonical-order assertion at the end is what makes this test fail on
    ``main`` deterministically.
    """
    permutations = [
        (READS_START, READS_END),
        (list(reversed(READS_START)), list(reversed(READS_END))),
        (READS_START[2:] + READS_START[:2], READS_END[1:] + READS_END[:1]),
        (sorted(READS_START), sorted(READS_END)),
    ]
    results = [
        _make_pattern(start, end if with_end else None).get_supporting_reads()
        for start, end in permutations
    ]
    assert all(r == results[0] for r in results), results
    assert results[0] == sorted(results[0])


def test_get_supporting_reads_is_stable_across_hash_seeds() -> None:
    """The defect itself: prove the order does not depend on the hash seed.

    Fails on ``main`` (``list(set(...))``); passes once the source sorts.
    """
    body = f"""
        from svirlpool.localassembly import SVpatterns
        from svirlpool.localassembly.SVprimitives import SVprimitive
        from svirlpool.svcalling import genotyping

        p = SVprimitive(
            ref_start=1000, ref_end=1000, read_start=100, read_end=600, size=500,
            sv_type=0, chr="chr1", repeatIDs=[], original_alt_sequences=["A" * 500],
            original_ref_sequences=[], samplename="sample1", consensusID="1.0",
            alignmentID=0, svID=0, aln_is_reverse=False,
            consensus_aln_interval=("chr1", 500, 1500),
            genotypeMeasurement=genotyping.GenotypeMeasurement(
                start_on_consensus=100,
                supporting_reads_start={READS_START!r},
                end_on_consensus=600,
                supporting_reads_end={READS_END!r},
            ),
        )
        pattern = SVpatterns.SVpatternInsertion(SVprimitives=[p])
        print(",".join(pattern.get_supporting_reads()))
    """
    outputs = {seed: _run_under_hash_seed(seed, body) for seed in HASH_SEEDS}
    distinct = set(outputs.values())
    assert len(distinct) == 1, (
        "SVpattern.get_supporting_reads() returned different orders under "
        f"different PYTHONHASHSEED values: {outputs}"
    )


def test_size_distortions_key_order_is_stable_across_hash_seeds() -> None:
    """``distortions_by_svPattern`` freezes the read order into a dict.

    That dict is stored as ``SVpattern.size_distortions`` and pickled into the
    svirltile DB, so its key order is persisted output.
    """
    body = f"""
        from svirlpool.localassembly import SVpatterns
        from svirlpool.localassembly.SVprimitives import SVprimitive
        from svirlpool.svcalling import genotyping

        class _NoDistortionConsensus:
            ID = "1.0"
            def get_consensus_distortions(self):
                return []

        p = SVprimitive(
            ref_start=1000, ref_end=1000, read_start=100, read_end=600, size=500,
            sv_type=0, chr="chr1", repeatIDs=[], original_alt_sequences=["A" * 500],
            original_ref_sequences=[], samplename="sample1", consensusID="1.0",
            alignmentID=0, svID=0, aln_is_reverse=False,
            consensus_aln_interval=("chr1", 500, 1500),
            genotypeMeasurement=genotyping.GenotypeMeasurement(
                start_on_consensus=100,
                supporting_reads_start={READS_START!r},
                end_on_consensus=600,
                supporting_reads_end={READS_END!r},
            ),
        )
        pattern = SVpatterns.SVpatternInsertion(SVprimitives=[p])
        distortions = SVpatterns.distortions_by_svPattern(
            svPattern=pattern,
            consensus=_NoDistortionConsensus(),
            distance_scale=5000.0,
            falloff=1.0,
        )
        print(",".join(distortions.keys()))
    """
    outputs = {seed: _run_under_hash_seed(seed, body) for seed in HASH_SEEDS}
    assert len(set(outputs.values())) == 1, (
        f"size_distortions key order depends on PYTHONHASHSEED: {outputs}"
    )
    # and the durable contract, checked in-process
    assert outputs[HASH_SEEDS[0]].split(",") == sorted(
        set(READS_START) | set(READS_END)
    )


# --------------------------------------------------------------------------- #
# VCF writer: the CONSENSUSIDs INFO field
# --------------------------------------------------------------------------- #
def _make_svcall(consensus_ids: list[str]) -> multisample_sv_calling.SVcall:
    import pickle

    return multisample_sv_calling.SVcall(
        genotypes={},
        passing=True,
        chrname="chr1",
        end=1000,
        start=1000,
        svtype="INS",
        svlen=42,
        pass_altreads=5,
        pass_gq=60,
        precise=True,
        mateid="",
        consensusIDs=list(consensus_ids),
        ref_sequence=pickle.dumps(""),
        alt_sequence=pickle.dumps("A" * 42),
    )


CONSENSUS_ID_PERMUTATIONS = [
    ["HG002:1.1", "HG002:1.0", "HG002:0.2"],
    ["HG002:0.2", "HG002:1.1", "HG002:1.0"],
    ["HG002:1.0", "HG002:0.2", "HG002:1.1"],
    ["HG002:0.2", "HG002:1.0", "HG002:1.1"],
]


def _consensusids_info_field(line: str) -> str:
    info = line.split("\t")[7]
    for field in info.split(";"):
        if field.startswith("CONSENSUSIDs="):
            return field.split("=", 1)[1]
    raise AssertionError(f"no CONSENSUSIDs field in INFO: {info}")


def test_vcf_consensusids_field_is_sorted() -> None:
    """``CONSENSUSIDs`` must be emitted in a canonical (sorted) order."""
    line = _make_svcall(CONSENSUS_ID_PERMUTATIONS[0]).to_vcf_line(
        vcfIDnumber=1,
        samplenames=[],
        covtrees={},
        refdict={"chr1:1000": "G"},
        symbolic_threshold=1000,
    )
    assert line is not None
    value = _consensusids_info_field(line)
    assert value.split(",") == sorted(value.split(","))


def test_vcf_consensusids_field_is_independent_of_input_order() -> None:
    """Whatever order the IDs were collected in, the VCF field is the same."""
    values = []
    for ids in CONSENSUS_ID_PERMUTATIONS:
        line = _make_svcall(ids).to_vcf_line(
            vcfIDnumber=1,
            samplenames=[],
            covtrees={},
            refdict={"chr1:1000": "G"},
            symbolic_threshold=1000,
        )
        assert line is not None
        values.append(_consensusids_info_field(line))
    assert len(set(values)) == 1, values
    assert values[0] == "HG002:0.2,HG002:1.0,HG002:1.1"


def test_svcall_to_log_id_consensusids_is_sorted() -> None:
    """``to_log_id`` is used for log-line diffing, so it must be canonical too."""
    log_ids = {_make_svcall(ids).to_log_id() for ids in CONSENSUS_ID_PERMUTATIONS}
    assert len(log_ids) == 1, log_ids
    assert "consensusIDs=HG002:0.2,HG002:1.0,HG002:1.1" in log_ids.pop()
