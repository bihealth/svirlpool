# %%
import json
from gzip import open as gzip_open
from pathlib import Path

import cattrs

from svirlpool.localassembly import consensus_class
from svirlpool.util.datatypes import Alignment
from svirlpool.util.util import (
    Direction,
    get_interval_on_ref_in_region,
    get_read_position_on_ref,
)

# %%

DATA_DIR = Path(__file__).parent / "data"

# %%

# files: INV.1.aln.json.gz  INV.1.consensus.json.gz  INV.2.aln.json.gz  INV.2.consensus.json.gz
# to save _process_alignment_file_for_core_intervals debug data:
# if consensus_obj.ID == "11.0":
#     cons = consensus_obj.unstructure()
#     alignment.unstructure()
#     consensus_path = "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.1.consensus.json"
#     with open(consensus_path, "w") as debug_f:
#         json.dump(cons, debug_f, indent=4)

# if consensus_obj.ID == "11.1":
#     cons = consensus_obj.unstructure()
#     alignment.unstructure()
#     consensus_path = "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.2.consensus.json"
#     with open(consensus_path, "w") as debug_f:
#         json.dump(cons, debug_f, indent=4)


def load_test_data(_name: str) -> consensus_class.Consensus:
    consensus_path = DATA_DIR / "consensus_class" / _name
    with gzip_open(consensus_path, "rt") as f:
        consensus_data = f.read()
        consensus: consensus_class.Consensus = cattrs.structure(
            json.loads(consensus_data), consensus_class.Consensus
        )
    return consensus


# This is used in the tests
def load_alignments(path: Path) -> list[Alignment]:
    with gzip_open(path, "rt") as f:
        data = json.load(f)
    alignments = cattrs.structure(data["alignments"], list[Alignment])
    return alignments


def generate_simulated_test_data():
    """Generate simulated test data and save as gzipped JSON files."""
    from svirlpool.util.util import (
        create_alignments_to_reference,
        delete_interval,
        generate_sequence,
        insertion,
        reverse_complement,
    )

    ref = generate_sequence(1000, seed=1)
    test_cases = {}

    # Simple forward read
    read_seq_forward = ref[100:800]
    alignments_forward = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_forward],
    )
    test_cases["simple_forward"] = {
        "alignments": [Alignment.from_pysam(aln) for aln in alignments_forward],
        "reference": ref,
    }

    # Reverse complemented read
    read_seq_rc = reverse_complement(ref[100:800])
    alignments_rc = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_rc],
    )
    test_cases["simple_reverse"] = {
        "alignments": [Alignment.from_pysam(aln) for aln in alignments_rc],
        "reference": ref,
    }

    # Read with deletion
    read_seq_del = delete_interval(ref, a=480, b=500)
    alignments_del = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_del],
    )
    test_cases["with_deletion"] = {
        "alignments": [Alignment.from_pysam(aln) for aln in alignments_del],
        "reference": ref,
    }

    # Read with deletion and reverse complement
    read_seq_del_rc = reverse_complement(delete_interval(ref, a=480, b=500))
    alignments_del_rc = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_del_rc],
    )
    test_cases["with_deletion_reverse"] = {
        "alignments": [Alignment.from_pysam(aln) for aln in alignments_del_rc],
        "reference": ref,
    }

    # Read with insertion
    read_seq_ins = insertion(ref, pos=500, sequence=["A"] * 20)
    alignments_ins = create_alignments_to_reference(reference=ref, reads=[read_seq_ins])
    test_cases["with_insertion"] = {
        "alignments": [Alignment.from_pysam(aln) for aln in alignments_ins],
        "reference": ref,
    }

    # Save each test case
    output_dir = DATA_DIR / "consensus_class"
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, data in test_cases.items():
        output_path = output_dir / f"simulated.{name}.json.gz"
        # Unstructure alignments
        unstructured_data = {
            "alignments": [aln.unstructure() for aln in data["alignments"]],
            "reference": data["reference"],
        }
        with gzip_open(output_path, "wt") as f:
            json.dump(unstructured_data, f, indent=2)
        print(f"Saved {name} to {output_path}")


# def debug_get_interval_on_ref_in_region_realdata():
#     # real data - the INV 15 data
#     cons = load_test_data("INV.15.consensus.json.gz")
#     print(
#         f"padding left: {cons.consensus_padding.padding_size_left}, padding right: {cons.consensus_padding.padding_size_right}, interval on padded seq: {cons.consensus_padding.consensus_interval_on_sequence_with_padding}"
#     )

#     alignments_path = DATA_DIR / "consensus_class" / "INV.15.alignments.json.gz"
#     alignments = [a.to_pysam() for a in load_alignments(alignments_path)]
#     # result_d = get_interval_on_ref_in_region(
#     #     a=alignments[0], start=107170497, end=107171181
#     # )

#     # this data is the bounds of the two candidate regions. It is wise to test if they can be traced back correctly from one of the three alignments.
#     cr_bounds = {
#         "end": {"read": 43639, "ref": 107171968},
#         "start": {"read": 48374, "ref": 107167422},
#     }
#     result_cr_bounds = {
#         i: get_interval_on_ref_in_region(
#             a=alignments[i],
#             start=cr_bounds["start"]["read"],
#             end=cr_bounds["end"]["read"],
#         )
#         for i in range(3)
#     }
#     # test for all three alignments if the ref positions match. At the end, both end positions should be found by at least one alignment.
#     #expected = cr_bounds["start"]["ref"], cr_bounds["end"]["ref"]
#     found_starts = set()
#     found_ends = set()
#     for i in range(3):
#         found_starts.add(result_cr_bounds[i][0])
#         found_ends.add(result_cr_bounds[i][1])
#     # assert expected[0] in found_starts, f"start position {expected[0]} not found in traced back positions {found_starts}"
#     # assert expected[1] in found_ends, f"end position {expected[1]} not found in traced back positions {found_ends}"

#     get_interval_on_ref_in_region(a=alignments[0], start=44536, end=44546)

#     # for the whole (39674, 41744) interval on the padded consensus, create a bed file to print the traced back reference positions
#     path_bed = DATA_DIR / "consensus_class" / "INV.15.consensus_interval.bed"
#     # chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
#     interval = cons.consensus_padding.consensus_interval_on_sequence_with_padding
#     colors = {0: "255,0,0", 1: "0,255,0", 2: "0,0,255"}
#     with open(path_bed, "w") as f:
#         for j in range(3):
#             strand = "+" if not alignments[j].is_reverse else "-"
#             isize = 1
#             for i in range(interval[0], interval[1], isize):
#                 start = i
#                 end = min(i + isize, interval[1])
#                 ref_start, ref_end = get_interval_on_ref_in_region(
#                     a=alignments[j], start=start, end=end
#                 )
#                 print(
#                     f"6\t{ref_start}\t{ref_end}\t{start}-{end}\t0\t{strand}\t{ref_start + 2}\t{ref_end - 2}\t{colors[j]}",
#                     file=f,
#                 )

#     # do a more precise test which checks all the individual positions in an interval on the alignments
#     from svirlpool.util.util import Direction, get_read_pitx_on_ref

#     path_bed_singles = (
#         DATA_DIR / "consensus_class" / "INV.15.consensus_interval_singles.bed"
#     )
#     with open(path_bed_singles, "w") as f:
#         for j in range(3):
#             strand = "+" if not alignments[j].is_reverse else "-"
#             for i in range(*interval, 5):
#                 ref_pos = get_read_pitx_on_ref(
#                     alignment=alignments[j], position=i, direction=Direction.NONE
#                 )[0]
#                 print(
#                     f"6\t{ref_pos}\t{ref_pos + 1}\t{i}\t0\t{strand}\t{ref_pos - 1}\t{ref_pos + 1}\t{colors[j]}",
#                     file=f,
#                 )


# %%


def test_get_interval_on_ref_in_region():
    alignments: list[Alignment] = load_alignments(
        DATA_DIR / "signalprocessing" / "alignments_to_rafs.dummy_inversion.json.gz"
    )
    # 2000M4000H - forward
    result_a = get_interval_on_ref_in_region(
        a=alignments[0].to_pysam(), start=1000, end=3000
    )
    expected_a = (1000, 2000)
    assert expected_a == result_a, f"got {result_a}, expected {expected_a}"
    # 2000H2000M2000H - reverse
    result_b = get_interval_on_ref_in_region(
        a=alignments[1].to_pysam(), start=1000, end=3000
    )
    expected_b = (3000, 4000)
    assert expected_b == result_b, f"got {result_b}, expected {expected_b}"
    # 4000S2000M - forward
    result_c = get_interval_on_ref_in_region(
        a=alignments[2].to_pysam(), start=3000, end=5000
    )
    expected_c = (4000, 5000)
    assert expected_c == result_c, f"got {result_c}, expected {expected_c}"


def test_get_read_position_on_ref_realdata():
    cons = load_test_data("INV.15.consensus.json.gz")
    alignments_path = DATA_DIR / "consensus_class" / "INV.15.alignments.json.gz"
    alignments = [a.to_pysam() for a in load_alignments(alignments_path)]

    results: dict[int, tuple[str, int, int]] = {
        i: consensus_class.get_consensus_core_alignment_interval_on_reference(
            consensus=cons, alignment=alignments[i]
        )
        for i in range(3)
    }
    expected = {
        0: ("6", 169094647, 169095161),
        1: ("6", 169093632, 169094647),
        2: ("6", 169093075, 169093632),
    }
    # assert results match expected
    for i in range(3):
        assert expected[i] == expected[i], f"got {results[i]}, expected {expected[i]}"


def test_get_read_position_on_ref_simulated():
    # Load pre-generated synthetic reads from files
    # Test simple forward read: pos 100 to 800 on reference
    alignments_forward = [
        aln.to_pysam()
        for aln in load_alignments(
            DATA_DIR / "consensus_class" / "simulated.simple_forward.json.gz"
        )
    ]
    result = get_read_position_on_ref(
        alignment=alignments_forward[0], position=100, direction=Direction.RIGHT
    )
    expected = 200
    assert result == expected, f"got {result}, expected {expected}"
    result_left = get_read_position_on_ref(
        alignment=alignments_forward[0], position=0, direction=Direction.RIGHT
    )
    expected_left = 100
    assert result_left == expected_left, f"got {result_left}, expected {expected_left}"
    result_right = get_read_position_on_ref(
        alignment=alignments_forward[0], position=766, direction=Direction.LEFT
    )
    expected_right = 800
    assert result_right == expected_right, (
        f"got {result_right}, expected {expected_right}"
    )

    # Test reverse complemented read
    alignments_rc = [
        aln.to_pysam()
        for aln in load_alignments(
            DATA_DIR / "consensus_class" / "simulated.simple_reverse.json.gz"
        )
    ]
    result_rc = get_read_position_on_ref(
        alignment=alignments_rc[0], position=100, direction=Direction.LEFT
    )
    expected_rc = 700
    assert result_rc == expected_rc, f"got {result_rc}, expected {expected_rc}"
    result_rc_left = get_read_position_on_ref(
        alignment=alignments_rc[0], position=0, direction=Direction.RIGHT
    )
    expected_rc_left = 800
    assert result_rc_left == expected_rc_left, (
        f"got {result_rc_left}, expected {expected_rc_left}"
    )
    result_rc_right = get_read_position_on_ref(
        alignment=alignments_rc[0], position=700, direction=Direction.LEFT
    )
    expected_rc_right = 100
    assert result_rc_right == expected_rc_right, (
        f"got {result_rc_right}, expected {expected_rc_right}"
    )

    # Test read with deletion (deleted positions 480-500)
    alignments_del = [
        aln.to_pysam()
        for aln in load_alignments(
            DATA_DIR / "consensus_class" / "simulated.with_deletion.json.gz"
        )
    ]
    for direction in [Direction.LEFT, Direction.RIGHT, Direction.NONE]:
        result_del_479 = get_read_position_on_ref(
            alignment=alignments_del[0], position=479, direction=direction
        )
        expected_del_479 = {
            Direction.LEFT: 479,
            Direction.RIGHT: 479,
            Direction.NONE: 479,
        }
        assert result_del_479 == expected_del_479[direction], (
            f"got {result_del_479}, expected {expected_del_479[direction]}, direction = {direction}"
        )
        result_del_480 = get_read_position_on_ref(
            alignment=alignments_del[0], position=480, direction=direction
        )
        expected_del_480 = {
            Direction.LEFT: 480,
            Direction.RIGHT: 500,
            Direction.NONE: 500,
        }
        assert result_del_480 == expected_del_480[direction], (
            f"got {result_del_480}, expected {expected_del_480[direction]}, direction = {direction}"
        )
        result_del_481 = get_read_position_on_ref(
            alignment=alignments_del[0], position=481, direction=direction
        )
        expected_del_481 = {
            Direction.LEFT: 501,
            Direction.RIGHT: 501,
            Direction.NONE: 501,
        }
        assert result_del_481 == expected_del_481[direction], (
            f"got {result_del_481}, expected {expected_del_481[direction]}, direction = {direction}"
        )

    # Test deletion with reverse complement
    alignments_del_rc = [
        aln.to_pysam()
        for aln in load_alignments(
            DATA_DIR / "consensus_class" / "simulated.with_deletion_reverse.json.gz"
        )
    ]
    for direction in [Direction.LEFT, Direction.RIGHT, Direction.NONE]:
        result_del_rc_500 = get_read_position_on_ref(
            alignment=alignments_del_rc[0], position=500, direction=direction
        )
        expected_del_rc_500 = {
            Direction.LEFT: 500,
            Direction.RIGHT: 480,
            Direction.NONE: 480,
        }
        assert expected_del_rc_500[direction] == result_del_rc_500, (
            f"got {result_del_rc_500}, expected {expected_del_rc_500[direction]}, direction = {direction}"
        )
        result_del_rc_499 = get_read_position_on_ref(
            alignment=alignments_del_rc[0], position=499, direction=direction
        )
        expected_del_rc_499 = {
            Direction.LEFT: 501,
            Direction.RIGHT: 501,
            Direction.NONE: 501,
        }
        assert expected_del_rc_499[direction] == result_del_rc_499, (
            f"got {result_del_rc_499}, expected {expected_del_rc_499[direction]}, direction = {direction}"
        )

    # Test read with insertion
    alignments_ins = [
        aln.to_pysam()
        for aln in load_alignments(
            DATA_DIR / "consensus_class" / "simulated.with_insertion.json.gz"
        )
    ]
    for direction in [Direction.LEFT, Direction.RIGHT, Direction.NONE]:
        result_ins_rc_501 = get_read_position_on_ref(
            alignment=alignments_ins[0], position=501, direction=direction
        )
        expected_ins_rc_501 = {
            Direction.LEFT: 499,
            Direction.RIGHT: 499,
            Direction.NONE: 499,
        }
        assert expected_ins_rc_501[direction] == result_ins_rc_501, (
            f"got {result_ins_rc_501}, expected {expected_ins_rc_501[direction]}, direction = {direction}"
        )
        result_ins_rc_500 = get_read_position_on_ref(
            alignment=alignments_ins[0], position=500, direction=direction
        )
        expected_ins_rc_500 = {
            Direction.LEFT: 499,
            Direction.RIGHT: 499,
            Direction.NONE: 499,
        }
        assert expected_ins_rc_500[direction] == result_ins_rc_500, (
            f"got {result_ins_rc_500}, expected {expected_ins_rc_500[direction]}, direction = {direction}"
        )
        result_ins_rc_499 = get_read_position_on_ref(
            alignment=alignments_ins[0], position=499, direction=direction
        )
        expected_ins_rc_499 = {
            Direction.LEFT: 499,
            Direction.RIGHT: 499,
            Direction.NONE: 499,
        }
        assert expected_ins_rc_499[direction] == result_ins_rc_499, (
            f"got {result_ins_rc_499}, expected {expected_ins_rc_499[direction]}, direction = {direction}"
        )


# %%


def test_get_consensus_core_alignment_interval_on_reference_inv():
    alignments_path = DATA_DIR / "consensus_class" / "INV.15.alignments.json.gz"

    alignments = [aln.to_pysam() for aln in load_alignments(alignments_path)]

    alignments_15 = [aln for aln in alignments if aln.query_name == "15.0"]

    cons1 = load_test_data("INV.15.consensus.json.gz")
    # for each alignment, find the core interval on reference
    results = {
        i: consensus_class.get_consensus_core_alignment_interval_on_reference(
            consensus=cons1, alignment=alignments_15[i]
        )
        for i in range(len(alignments_15))
    }

    expected = {
        0: ("6", 169094647, 169095161),
        1: ("6", 169093632, 169094647),
        2: ("6", 169093075, 169093632),
    }

    for i, res in results.items():
        assert res == expected[i], f"Alignment {i}: got {res}, expected {expected[i]}"
