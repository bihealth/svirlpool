"""--min-signal-support and --satellite-depth-factor of signalstrength_to_crs."""

import numpy as np

from svirlpool.candidateregions import signalstrength_to_crs as S2C


def support(rows):
    """rows: (start, size, sv_type, repeatID, read)"""
    starts, sizes, types, repeats, reads = (np.array(c) for c in zip(*rows, strict=True))
    return S2C.signal_support(
        starts=starts,
        sizes=sizes,
        sv_types=types,
        repeat_ids=repeats,
        reads=np.unique(reads, return_inverse=True)[1],
    ).tolist()


def test_reads_agreeing_on_one_sv_confirm_each_other():
    rows = [(1000 + 10 * i, 300, 0, -1, f"r{i}") for i in range(3)]
    assert support(rows) == [2, 2, 2]


def test_a_read_does_not_confirm_itself():
    rows = [(1000, 300, 0, -1, "r1"), (1020, 300, 0, -1, "r1")]
    assert support(rows) == [0, 0]


def test_dissimilar_sizes_far_signals_and_other_types_do_not_confirm():
    rows = [
        (1000, 300, 0, -1, "a"),
        (1010, 100, 0, -1, "b"),  # size ratio 1/3
        (1500, 300, 0, -1, "c"),  # 500 bp away > max(100, 150)
        (1000, 300, 1, -1, "d"),  # a deletion
    ]
    assert support(rows) == [0, 0, 0, 0]


def test_deletion_ends_count_as_one_type():
    rows = [(1000, 300, 1, -1, "a"), (1050, 300, 2, -1, "b")]
    assert support(rows) == [1, 1]


def test_the_same_tandem_repeat_counts_as_near():
    rows = [(1000, 300, 0, 7, "a"), (5000, 280, 0, 7, "b"), (9000, 300, 0, -1, "c")]
    assert support(rows) == [1, 1, 0]


def test_long_tandem_repeats_are_merged(tmp_path):
    bed = tmp_path / "chr1_repeats.bed"
    bed.write_text(
        "chr1\t0\t5000\tx\nchr1\t1000\t12000\tx\nchr1\t11000\t20000\tx\n"
        "chr1\t30000\t32000\tx\n"
    )
    intervals = S2C.load_long_tandem_repeats(bed, 6000)
    assert intervals.tolist() == [[1000, 20000]]
    assert S2C.load_long_tandem_repeats(tmp_path / "missing.bed", 6000).shape == (0, 2)


def test_covered_fraction():
    intervals = np.array([[1000, 2000], [3000, 4000]])
    assert S2C.covered_fraction(1500, 3500, intervals) == 0.5
    assert S2C.covered_fraction(1200, 1800, intervals) == 1.0
    assert S2C.covered_fraction(5000, 6000, intervals) == 0.0
    assert S2C.covered_fraction(5000, 6000, np.zeros((0, 2), dtype=np.int64)) == 0.0
