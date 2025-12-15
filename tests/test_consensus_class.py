import json
from gzip import open as gzip_open
from gzip import open as gzopen
from pathlib import Path

import cattrs
from pysam import AlignmentFile

from svirlpool.localassembly import consensus_class
from svirlpool.util.datatypes import Alignment
from svirlpool.util.util import get_interval_on_ref_in_region

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


def test_get_consensus_core_alignment_interval_on_reference_inv():
    # test 11.0 and 11.1
    alignments_path = DATA_DIR / "consensus_class" / "inv.11.bam"
    try:
        alignments = list(AlignmentFile(str(alignments_path)).fetch())
    except Exception as e:
        raise RuntimeError(f"Error loading alignments from {alignments_path}: {e}")
    # test 11.0
    alignments_110 = [aln for aln in alignments if aln.query_name == "11.0"]

    cons1 = load_test_data("INV.110.consensus.json.gz")
    # for each alignment, find the core interval on reference
    results = {
        i: consensus_class.get_consensus_core_alignment_interval_on_reference(
            consensus=cons1, alignment=alignments_110[i]
        )
        for i in range(len(alignments_110))
    }

    expected = {
        0: ("6", 169093103, 169093632),
        1: ("6", 169093632, 169094647),
        2: ("6", 169094647, 169095172),
    }

    for i, res in results.items():
        assert res == expected[i], f"Alignment {i}: got {res}, expected {expected[i]}"

    alignments_111 = [aln for aln in alignments if aln.query_name == "11.1"]

    cons2 = load_test_data("INV.111.consensus.json.gz")
    results2 = {
        i: consensus_class.get_consensus_core_alignment_interval_on_reference(
            consensus=cons2, alignment=alignments_111[i]
        )
        for i in range(len(alignments_111))
    }

    expected2 = {
        0: ("6", 169093103, 169093621),
        1: ("6", 169093626, 169094647),
        2: ("6", 169094647, 169094656),
    }

    for i, res in results2.items():
        assert res == expected2[i], f"Alignment {i}: got {res}, expected {expected2[i]}"
