import json
from gzip import open as gzip_open
from gzip import open as gzopen
from pathlib import Path

import cattrs
import pytest
from pysam import AlignmentFile

from svirlpool.localassembly import consensus_class
from svirlpool.util.datatypes import Alignment
from svirlpool.util.util import (Direction, get_interval_on_ref_in_region,
                                 get_read_position_on_ref)

DATA_DIR = Path(__file__).parent / "data"

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
    consensus_path = (
        DATA_DIR / "consensus_class" / _name
    )  # INV.110.consensus.json.gz or INV.111.consensus.json.gz
    # load the gzipped json dumped structured objects
    with gzip_open(consensus_path, "rt") as f:
        consensus_data = f.read()
        consensus: consensus_class.Consensus = cattrs.structure(
            json.loads(consensus_data), consensus_class.Consensus
        )
    return consensus


# This is used in the tests
def load_alignments(path: Path) -> list[Alignment]:
    with gzopen(path, "rt") as f:
        data = json.load(f)
    alignments = cattrs.structure(data["alignments"], list[Alignment])
    return alignments


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
    
def test_get_interval_on_ref_in_region_realdata():
    # real data - the INV 7 data
    cons7 = load_test_data("INV.7.consensus.json.gz")
    print(f"padding left: {cons7.consensus_padding.padding_size_left}, padding right: {cons7.consensus_padding.padding_size_right}, interval on padded seq: {cons7.consensus_padding.consensus_interval_on_sequence_with_padding}")
    
    alignments_path = DATA_DIR / "consensus_class" / "INV.7.alignments.json.gz"
    alignments = [a.to_pysam() for a in load_alignments(alignments_path)]
    result_d = get_interval_on_ref_in_region(a=alignments[0], start=107170497, end=107171181)
    
    # this data is the bounds of the two candidate regions. It is wise to test if they can be traced back correctly from one of the three alignments.
    cr_bounds = {"end":{"read":43639, "ref":107171968},
                  "start":{"read":48374, "ref": 107167422}}
    result_cr_bounds = {i: get_interval_on_ref_in_region(a=alignments[i], start=cr_bounds["start"]["read"], end=cr_bounds["end"]["read"]) for i in range(3)}
    # test for all three alignments if the ref positions match. At the end, both end positions should be found by at least one alignment.
    expected = cr_bounds["start"]["ref"], cr_bounds["end"]["ref"]
    found_starts = set()
    found_ends = set()
    for i in range(3):
        found_starts.add(result_cr_bounds[i][0])
        found_ends.add(result_cr_bounds[i][1])
    assert expected[0] in found_starts, f"start position {expected[0]} not found in traced back positions {found_starts}"
    assert expected[1] in found_ends, f"end position {expected[1]} not found in traced back positions {found_ends}"
    
    get_interval_on_ref_in_region(a=alignments[0], start=44536, end=44546)
    
    # for the whole (39674, 41744) interval on the padded consensus, create a bed file to print the traced back reference positions
    path_bed = DATA_DIR / "consensus_class" / "INV.7.consensus_interval.bed"
    # chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
    interval = cons7.consensus_padding.consensus_interval_on_sequence_with_padding
    colors = {0: "255,0,0", 1: "0,255,0", 2: "0,0,255"}
    with open(path_bed, "w") as f:
        for j in range(3):
            strand = "+" if not alignments[j].is_reverse else "-"
            isize = 10
            for i in range(interval[0], interval[1], isize):
                start = i
                end = min(i + isize, interval[1])
                ref_start, ref_end = get_interval_on_ref_in_region(a=alignments[j], start=start, end=end)
                print(f"6\t{ref_start}\t{ref_end}\t{start}-{end}\t0\t{strand}\t{ref_start+2}\t{ref_end-2}\t{colors[j]}", file=f)
    
    # do a more precise test which checks all the individual positions in an interval on the alignments
    from svirlpool.util.util import Direction, get_read_pitx_on_ref
    path_bed_singles = DATA_DIR / "consensus_class" / "INV.7.consensus_interval_singles.bed"
    with open(path_bed_singles, "w") as f:
        for j in range(3):
            strand = "+" if not alignments[j].is_reverse else "-"
            for i in range(*interval, 5):
                ref_pos = get_read_pitx_on_ref(alignment=alignments[j], position=i, direction=Direction.NONE)[0]
                print(f"6\t{ref_pos}\t{ref_pos+1}\t{i}\t0\t{strand}\t{ref_pos-1}\t{ref_pos+1}\t{colors[j]}", file=f)
    
    #assert expected_d == result_d, f"got {result_d}, expected {expected_d}"


def test_get_read_position_on_ref_simulated():
    # generate synthetic reads from a reference sequence
    # add insertions and deletions and test if the correct positions are traced back
    from svirlpool.util.util import (create_alignments_to_reference,
                                     delete_interval, duplication,
                                     generate_sequence, inversion,
                                     reverse_complement)
    ref = generate_sequence(1000, seed=1)
    # create a simple read on the + strand from pos 100 to 800 and trace back position 100 on the read to get 200 on the reference
    read_seq = ref[100:800]
    alignments = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq],
    )
    result = get_read_position_on_ref(alignment=alignments[0], position=100, direction=Direction.RIGHT)
    expected = 200
    assert result == expected, f"got {result}, expected {expected}"
    result_left = get_read_position_on_ref(alignment=alignments[0], position=0, direction=Direction.RIGHT)
    expected_left = 100
    assert result_left == expected_left, f"got {result_left}, expected {expected_left}"
    result_right = get_read_position_on_ref(alignment=alignments[0], position=766, direction=Direction.LEFT)
    expected_right = 800
    assert result_right == expected_right, f"got {result_right}, expected {expected_right}"
    
    # now test a reverse complemented read
    read_seq_rc = reverse_complement(ref[100:800])
    alignments_rc = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_rc],
    )
    result_rc = get_read_position_on_ref(alignment=alignments_rc[0], position=100, direction=Direction.LEFT)
    expected_rc = 700
    assert result_rc == expected_rc, f"got {result_rc}, expected {expected_rc}"
    result_rc_left = get_read_position_on_ref(alignment=alignments_rc[0], position=0, direction=Direction.RIGHT)
    expected_rc_left = 800
    assert result_rc_left == expected_rc_left, f"got {result_rc_left}, expected {expected_rc_left}"
    result_rc_right = get_read_position_on_ref(alignment=alignments_rc[0], position=700, direction=Direction.LEFT)
    expected_rc_right = 100
    assert result_rc_right == expected_rc_right, f"got {result_rc_right}, expected {expected_rc_right}"
    
    # test to trace back a position in a deletion
    read_seq_del = delete_interval(ref, a=480, b=500)
    alignments_del = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_del],
    )
    # test position 480 nd 481 on the read should both map to 500 on the reference
    result_del_480 = get_read_position_on_ref(alignment=alignments_del[0], position=480, direction=Direction.LEFT)
    expected_del_480 = 480
    assert result_del_480 == expected_del_480, f"got {result_del_480}, expected {expected_del_480}"
    result_del_481 = get_read_position_on_ref(alignment=alignments_del[0], position=481, direction=Direction.RIGHT)
    expected_del_481 = 501 # or should this be 500?
    assert result_del_481 == expected_del_481, f"got {result_del_481}, expected {expected_del_481}"
    # test with Direction.NONE
    result_del_none_480 = get_read_position_on_ref(alignment=alignments_del[0], position=490, direction=Direction.NONE)
    expected_del_none_480 = 490
    assert result_del_none_480 == expected_del_none_480, f"got {result_del_none_480}, expected {expected_del_none_480}"
    result_del_none_481 = get_read_position_on_ref(alignment=alignments_del[0], position=491, direction=Direction.NONE)
    expected_del_none_481 = 491
    assert result_del_none_481 == expected_del_none_481, f"got {result_del_none_481}, expected {expected_del_none_481}"
    
    # test the deletion with a reverse complemented read
    read_seq_del_rc = reverse_complement(delete_interval(ref, a=480, b=500))
    alignments_del_rc = create_alignments_to_reference(
        reference=ref,
        reads=[read_seq_del_rc],
    )
    # test position 500 and 501 on the read should map to 500, 480 on the reference
    result_del_rc_500 = get_read_position_on_ref(alignment=alignments_del_rc[0], position=499, direction=Direction.LEFT)
    expected_del_rc_500 = 501
    assert expected_del_rc_500 == result_del_rc_500, f"got {result_del_rc_500}, expected {expected_del_rc_500}"
    result_del_rc_501 = get_read_position_on_ref(alignment=alignments_del_rc[0], position=500, direction=Direction.RIGHT)
    expected_del_rc_501 = 480
    assert expected_del_rc_501 == result_del_rc_501, f"got {result_del_rc_501}, expected {expected_del_rc_501}"
    # test with Direction.NONE
    result_del_rc_none_510 = get_read_position_on_ref(alignment=alignments_del_rc[0], position=510, direction=Direction.NONE)
    expected_del_rc_none_510 = 510
    assert expected_del_rc_none_510 == result_del_rc_none_510, f"got {result_del_rc_none_510}, expected {expected_del_rc_none_510}"
    result_del_rc_none_511 = get_read_position_on_ref(alignment=alignments_del_rc[0], position=511, direction=Direction.NONE)
    expected_del_rc_none_511 = 509
    assert expected_del_rc_none_511 == result_del_rc_none_511, f"got {result_del_rc_none_511}, expected {expected_del_rc_none_511}"
    


def test_get_consensus_core_alignment_interval_on_reference_inv():
    # # test 11.0 and 11.1
    # alignments_path = DATA_DIR / "consensus_class" / "INV.11.alignments.json.gz"

    # alignments = [aln.to_pysam() for aln in load_alignments(
    #     alignments_path
    # )]

    # alignments_110 = [aln for aln in alignments if aln.query_name == "11.0"]

    # cons1 = load_test_data("INV.110.consensus.json.gz")
    # # for each alignment, find the core interval on reference
    # results = {
    #     i: consensus_class.get_consensus_core_alignment_interval_on_reference(
    #         consensus=cons1, alignment=alignments_110[i]
    #     )
    #     for i in range(len(alignments_110))
    # }

    # expected = {
    #     0: ("6", 169093103, 169093632),
    #     1: ("6", 169093632, 169094647),
    #     2: ("6", 169094647, 169095172),
    # }

    # for i, res in results.items():
    #     assert res == expected[i], f"Alignment {i}: got {res}, expected {expected[i]}"

    # alignments_111 = [aln for aln in alignments if aln.query_name == "11.1"]

    # cons2 = load_test_data("INV.111.consensus.json.gz")
    # results2 = {
    #     i: consensus_class.get_consensus_core_alignment_interval_on_reference(
    #         consensus=cons2, alignment=alignments_111[i]
    #     )
    #     for i in range(len(alignments_111))
    # }

    # expected2 = {
    #     0: ("6", 169093103, 169093621),
    #     1: ("6", 169093626, 169094647),
    #     2: ("6", 169094647, 169094656),
    # }

    # for i, res in results2.items():
    #     assert res == expected2[i], f"Alignment {i}: got {res}, expected {expected2[i]}"
    
    # test 7.0
    alignments_path = DATA_DIR / "consensus_class" / "INV.7.alignments.json.gz"
    alignments_70 = [aln.to_pysam() for aln in load_alignments(
        alignments_path
    )]
    cons3 = load_test_data("INV.7.consensus.json.gz")
    results3 = {
        i: consensus_class.get_consensus_core_alignment_interval_on_reference(
            consensus=cons3, alignment=alignments_70[i]
        )
        for i in range(len(alignments_70))
    }
    
    print(results3)
    
    expected3 = {
        0: ('6', 107171181, 107170880),
        1: ('6', 107169206, 107170880),
        2: ('6', 107167422, 107169206)}
    for i, res in results3.items():
        assert res == expected3[i], f"Alignment {i}: got {res}, expected {expected3[i]}"

test_get_read_position_on_ref_simulated()
# test_get_interval_on_ref_in_region_realdata()
#test_get_consensus_core_alignment_interval_on_reference_inv()

