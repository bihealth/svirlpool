import json
from gzip import open as gzip_open
from gzip import open as gzopen
from pathlib import Path

import cattrs

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
#     alignment_path = "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.1.aln.json"
#     with open(alignment_path, "w") as debug_f:
#         json.dump(alignment.unstructure(), debug_f, indent=4)

# if consensus_obj.ID == "11.1":
#     cons = consensus_obj.unstructure()
#     alignment.unstructure()
#     consensus_path = "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.2.consensus.json"
#     with open(consensus_path, "w") as debug_f:
#         json.dump(cons, debug_f, indent=4)
#     alignment_path = "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.2.aln.json"
#     with open(alignment_path, "w") as debug_f:
#         json.dump(alignment.unstructure(), debug_f, indent=4)


def load_test_data(_name: str) -> consensus_class.Consensus:
    consensus_path = DATA_DIR / "consensus_class" / _name # INV.110.consensus.json.gz or INV.111.consensus.json.gz
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
        a=alignments[0].to_pysam(),
        start=1000,
        end=3000
    )
    expected_a = (1000, 2000)
    assert expected_a == result_a, f"got {result_a}, expected {expected_a}"
    # 2000H2000M2000H - reverse
    result_b = get_interval_on_ref_in_region(
        a=alignments[1].to_pysam(),
        start=1000,
        end=3000
    )
    expected_b = (3000, 4000)
    assert expected_b == result_b, f"got {result_b}, expected {expected_b}"
    # 4000S2000M - forward
    result_c = get_interval_on_ref_in_region(
        a=alignments[2].to_pysam(),
        start=3000,
        end=5000
    )
    expected_c = (4000, 5000)
    assert expected_c == result_c, f"got {result_c}, expected {expected_c}"

    
def test_get_consensus_core_alignment_interval_on_reference_inv():
    # test 11.0 and 11.1
    alignments_path = DATA_DIR / "consensus_class" / "inv.11.bam"
    alignments = list(pysam.AlignmentFile(alignments_path).fetch())
    # test 11.0
    alignments_110 = [aln for aln in alignments if aln.query_name == "11.0"]
    
    cons1 = load_test_data("INV.110.consensus.json.gz")
    # for each alignment, find the core interval on reference
    result = consensus_class.get_consensus_core_alignment_interval_on_reference(
        consensus=cons1,
        alignment=aln1.to_pysam(),
    )
    expected = ("6", 169093104, 169093633)
    assert result is not None
    chrom, start, end = result
    assert expected[0] == chrom
    assert abs(start - expected[1]) <= 100 and abs(end - expected[2]) <= 100, (
        f"start {start} vs expected {expected[1]}; end {end} vs expected {expected[2]}"
    )
    # test 11.1
    cons2, aln2 = load_test_data("INV.2")
    result2 = consensus_class.get_consensus_core_alignment_interval_on_reference(
        consensus=cons2,
        alignment=aln2.to_pysam(),
    )
    expected2 = ("6", 169093104, 169093633)
    assert result2 is not None
    chrom2, start2, end2 = result2
    assert expected2[0] == chrom2
    assert abs(start2 - expected2[1]) <= 100 and abs(end2 - expected2[2]) <= 100, (
        f"start {start2} vs expected {expected2[1]}; end {end2} vs expected {expected2[2]}"
    )

test_get_consensus_core_alignment_interval_on_reference_inv()
