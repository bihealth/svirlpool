"""Tests for :mod:`svirlpool.candidateregions.crs_to_batches`."""

from __future__ import annotations

from pathlib import Path

from svirlpool.candidateregions import crs_to_batches


def _ext(chrs, start, end):
    return (chrs, start, end)


def test_build_batches_groups_adjacent_containers() -> None:
    extents = {
        1: _ext(["chr1"], 100, 200),
        2: _ext(["chr1"], 300, 400),
        3: _ext(["chr1"], 500, 600),
        4: _ext(["chr2"], 100, 200),
    }
    batches = crs_to_batches.build_batches(
        extents=extents, batch_size=100, max_batch_span=10_000
    )
    # chr1 has 3 containers, chr2 has 1; expect 2 batches total
    assert len(batches) == 2
    chr1_batch = [b for b in batches if b["chr"] == "chr1"][0]
    chr2_batch = [b for b in batches if b["chr"] == "chr2"][0]
    assert chr1_batch["crIDs"] == [1, 2, 3]
    assert chr2_batch["crIDs"] == [4]


def test_build_batches_respects_batch_size_cap() -> None:
    extents = {i: _ext(["chr1"], i * 100, i * 100 + 50) for i in range(1, 6)}
    batches = crs_to_batches.build_batches(
        extents=extents, batch_size=2, max_batch_span=10_000
    )
    sizes = sorted(len(b["crIDs"]) for b in batches)
    assert sizes == [1, 2, 2]


def test_build_batches_respects_max_span() -> None:
    extents = {
        1: _ext(["chr1"], 100, 200),
        2: _ext(["chr1"], 5_000, 5_100),
        3: _ext(["chr1"], 9_000, 9_100),
    }
    batches = crs_to_batches.build_batches(
        extents=extents, batch_size=100, max_batch_span=5_000
    )
    # 1 and 2 fit within 5kb span; 3 starts a new batch.
    crIDs_per_batch = sorted(tuple(b["crIDs"]) for b in batches)
    assert crIDs_per_batch == [(1, 2), (3,)]


def test_wide_container_becomes_singleton_batch() -> None:
    extents = {
        1: _ext(["chr1"], 100, 200),
        2: _ext(["chr1", "chr2"], 100, 200),  # multi-chromosome -> wide
        3: _ext(["chr3"], 0, 100_000_000),  # huge span -> wide
        4: _ext(["chr1"], 300, 400),
    }
    batches = crs_to_batches.build_batches(
        extents=extents, batch_size=100, max_batch_span=1_000_000
    )
    # 2 and 3 are wide -> their own batches; 1 and 4 group on chr1.
    wide_ids = {tuple(b["crIDs"]) for b in batches if len(b["crIDs"]) == 1}
    assert (2,) in wide_ids
    assert (3,) in wide_ids
    local = [b for b in batches if b["chr"] == "chr1" and len(b["crIDs"]) > 1]
    assert len(local) == 1
    assert local[0]["crIDs"] == [1, 4]


def test_build_batches_assigns_unique_sequential_ids() -> None:
    extents = {i: _ext(["chr1"], i * 100, i * 100 + 50) for i in range(1, 8)}
    batches = crs_to_batches.build_batches(
        extents=extents, batch_size=3, max_batch_span=10_000
    )
    ids = [b["batch_id"] for b in batches]
    assert ids == sorted(set(ids))
    assert ids[0] == 0


def test_write_outputs_roundtrips_through_consensus_loader(tmp_path: Path) -> None:
    # End-to-end: write the TSV and re-parse it via the consensus helper.
    from svirlpool.localassembly import consensus

    extents = {
        1: _ext(["chr1"], 100, 200),
        2: _ext(["chr1"], 300, 400),
        7: _ext(["chr2"], 100, 200),
    }
    batches = crs_to_batches.build_batches(
        extents=extents, batch_size=100, max_batch_span=10_000
    )
    tsv = tmp_path / "batches.tsv"
    ids_file = tmp_path / "ids.txt"
    crs_to_batches.write_outputs(batches=batches, tsv_path=tsv, ids_path=ids_file)

    expected_per_batch = {b["batch_id"]: b["crIDs"] for b in batches}
    for batch_id, expected_crIDs in expected_per_batch.items():
        got = consensus._load_crIDs_from_batch_tsv(tsv, batch_id)
        assert got == expected_crIDs

    with open(ids_file) as f:
        listed = [int(line.strip()) for line in f if line.strip()]
    assert sorted(listed) == sorted(expected_per_batch.keys())
