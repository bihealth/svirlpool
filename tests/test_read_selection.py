"""--read-selection-factor: the reads kept in a crowded candidate region."""

import json
import sqlite3
from types import SimpleNamespace

import pysam
import pytest

from svirlpool.candidateregions import container_depth
from svirlpool.localassembly import consensus, read_selection

HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 10_000_000}]})
S, E = 100_000, 102_000  # the CR


def aln(
    name, start, length, clip_left=0, clip_right=0, reverse=False, sa=None, supp=False
):
    a = pysam.AlignedSegment(HEADER)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = start
    cigar = []
    if clip_left:
        cigar.append((4, clip_left))
    cigar.append((0, length))
    if clip_right:
        cigar.append((4, clip_right))
    a.cigartuples = cigar
    a.query_sequence = "A" * (clip_left + length + clip_right)
    a.flag = (16 if reverse else 0) | (2048 if supp else 0)
    a.mapping_quality = 60
    if sa:
        a.set_tag("SA", sa)
    return a


def test_single_alignment_crossing():
    p = read_selection.read_placements(aln("r", S - 3000, 7000))
    assert read_selection.crossing(p, "chr1", S, E) == ("crossing", 2000)


def test_collinear_split_alignments_cross():
    # left part to S + 500, an expanded allele of 5 kb on the read, right part from E - 500
    left = aln("r", S - 4000, 4500, clip_right=12_000)
    sa = f"chr1,{E - 500 + 1},+,9500S6500M,60,0;"
    left.set_tag("SA", sa)
    p = read_selection.read_placements(left)
    assert read_selection.crossing(p, "chr1", S, E) == ("crossing", 4000)


def test_one_sided_and_inside():
    p = read_selection.read_placements(aln("r", S - 6000, 6500, clip_right=3000))
    assert read_selection.crossing(p, "chr1", S, E) == ("one-sided", 6000)
    p = read_selection.read_placements(
        aln("f", S + 200, 800, clip_left=5000, clip_right=5000)
    )
    assert read_selection.crossing(p, "chr1", S, E) == ("inside", -1)


def test_few_reads_are_all_kept():
    alns = [aln(f"r{i}", S + 100, 500) for i in range(5)]
    keep, n = read_selection.select_reads(alns, "chr1", S, E, k=10)
    assert keep is None and n["kept"] == 5


def test_crowded_cr_keeps_the_best_anchored_reads():
    alns = [
        aln(f"cross{i}", S - 1000 * (i + 1), E - S + 2000 * (i + 1)) for i in range(4)
    ]
    alns += [aln("left", S - 9000, 9500, clip_right=2000)]
    alns += [
        aln(f"frag{i}", S + 100, 900, clip_left=4000, clip_right=4000)
        for i in range(20)
    ]
    keep, n = read_selection.select_reads(alns, "chr1", S, E, k=3)
    assert keep == {"cross3", "cross2", "cross1"}
    assert n == {"reads": 25, "crossing": 4, "one_sided": 1, "kept": 3}
    keep, _ = read_selection.select_reads(alns, "chr1", S, E, k=6)
    assert keep == {"cross0", "cross1", "cross2", "cross3", "left"}  # fragments never


def test_crowded_cr_without_crossing_read_is_dropped():
    alns = [
        aln(f"frag{i}", S + 100, 900, clip_left=4000, clip_right=4000)
        for i in range(20)
    ]
    keep, n = read_selection.select_reads(alns, "chr1", S, E, k=5)
    assert keep == set() and n["kept"] == 0


def test_median_container_depth(tmp_path):
    db = tmp_path / "c.db"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE containers (crID INTEGER, data TEXT)")
    for i, covs in enumerate(([20, 21, 19], [18, 20], [80, 90], [22])):
        crs = [{"sv_signals": [{"coverage": c} for c in covs]}]
        con.execute(
            "INSERT INTO containers VALUES (?, ?)", (i, json.dumps({"crs": crs}))
        )
    con.commit()
    # container medians 20, 19, 85, 22
    assert container_depth.median_container_depth(db) == pytest.approx(21)


def _containers_db(path, depths):
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE containers (crID INTEGER, data TEXT)")
    for i, d in enumerate(depths):
        crs = [{"sv_signals": [{"coverage": d}]}]
        con.execute(
            "INSERT INTO containers VALUES (?, ?)", (i, json.dumps({"crs": crs}))
        )
    con.commit()


def test_k_from_args(tmp_path):
    f = tmp_path / "depth.txt"
    f.write_text("20\n")
    args = SimpleNamespace(
        read_selection_factor=2.0, median_depth_file=str(f), input=None
    )
    assert consensus.read_selection_k_from_args(args) == 40
    args.read_selection_factor = 0
    assert consensus.read_selection_k_from_args(args) == 0
    # without a depth file: the median of the containers database
    db = tmp_path / "c.db"
    _containers_db(db, [18, 30, 30])
    args = SimpleNamespace(read_selection_factor=2.0, median_depth_file=None, input=db)
    assert consensus.read_selection_k_from_args(args) == 60
    # an empty database: no depth, all reads
    empty = tmp_path / "e.db"
    _containers_db(empty, [])
    args.input = empty
    assert consensus.read_selection_k_from_args(args) == 0


def test_settings_are_on_by_default():
    from svirlpool.__main__ import get_parser

    run = get_parser().parse_args(
        [
            "run",
            "--samplename",
            "s",
            "--workdir",
            "w",
            "--alignments",
            "a.bam",
            "--reference",
            "r.fa",
            "--trf",
            "t.bed",
            "--mononucleotides",
            "m.bed",
            "--threads",
            "1",
        ]
    )
    assert run.read_selection_factor == 3
    assert run.container_time_limit == 120
    # reference-SNV phasing, else all-vs-all phasing (no KMeans route)
    assert run.clustering_strategy == "accurate"
    assert run.phasing_sites == "tiered"
    cons = consensus.get_consensus_parser().parse_args(
        [
            "-s",
            "s",
            "-i",
            "c.db",
            "-a",
            "a.bam",
            "-cn",
            "cn.bed.gz",
            "-o",
            "o",
            "-r",
            "r.fa",
        ]
    )
    assert cons.read_selection_factor == 3
    assert cons.container_time_limit == 120
    assert cons.clustering_strategy == "accurate"
    assert cons.phasing_sites == "tiered"


def test_workflow_always_passes_the_factor():
    """main.smk must pass --read-selection-factor even when it is 0: without
    it the consensus CLI applies its own default (2), and 0 could not switch
    the read selection off."""
    import re
    from pathlib import Path

    import svirlpool

    smk = (Path(svirlpool.__file__).parent / "workflows" / "main.smk").read_text()
    expr = re.search(r"read_selection_arg=\((.*?)\n        \),", smk, re.S).group(1)
    for factor, expected in (
        (0, "--read-selection-factor 0"),
        (2, "--read-selection-factor 2 --median-depth-file"),
    ):
        arg = eval("(" + expr + ")", {"read_selection_factor": factor})  # noqa: S307
        assert arg == expected or arg.startswith(expected + " ")
