import json
from gzip import open as gzip_open
from pathlib import Path

import cattrs

from svirlpool.localassembly import consensus_class
from svirlpool.util.datatypes import Alignment

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


def load_test_data(prefix: str) -> tuple[consensus_class.Consensus, Alignment]:
    # prefix is like "INV.1"
    filename_consensus = f"{prefix}.consensus.json.gz"
    filename_alignment = f"{prefix}.aln.json.gz"
    consensus_path = DATA_DIR / "consensus_class" / filename_consensus
    # load the gzipped json dumped structured objects
    with gzip_open(consensus_path, "rt") as f:
        consensus_data = f.read()
        consensus: consensus_class.Consensus = cattrs.structure(
            json.loads(consensus_data), consensus_class.Consensus
        )
    # load the corresponding alignment
    alignment_path = DATA_DIR / "consensus_class" / filename_alignment
    with gzip_open(alignment_path, "rt") as f:
        alignment_data = f.read()
        alignment: Alignment = cattrs.structure(json.loads(alignment_data), Alignment)
    return consensus, alignment


# consensus_class.get_consensus_core_alignment_interval_on_reference


def test_get_consensus_core_alignment_interval_on_reference_inv():
    # test 11.0
    cons1, aln1 = load_test_data("INV.1")
    result = consensus_class.get_consensus_core_alignment_interval_on_reference(
        consensus=cons1,
        alignment=aln1.to_pysam(),
    )
    expected = ("6", 169093172, 169095124)
    assert result is not None
    chrom, start, end = result
    assert chrom == expected[0]
    assert abs(start - expected[1]) <= 100, (
        f"start {start} vs expected {expected[1]}; consensus padding: {cons1.consensus_padding}"
    )
    assert abs(end - expected[2]) <= 100, (
        f"end {end} vs expected {expected[2]}; consensus padding: {cons1.consensus_padding}"
    )
    # test 11.1
    cons2, aln2 = load_test_data("INV.2")
    result2 = consensus_class.get_consensus_core_alignment_interval_on_reference(
        consensus=cons2,
        alignment=aln2.to_pysam(),
    )
    expected2 = ("6", 169093172, 169095124)
    assert result2 is not None
    chrom2, start2, end2 = result2
    assert chrom2 == expected2[0]
    assert abs(start2 - expected2[1]) <= 100, (
        f"start {start2} vs expected {expected2[1]}; consensus padding: {cons2.consensus_padding}"
    )
    assert abs(end2 - expected2[2]) <= 100, (
        f"end {end2} vs expected {expected2[2]}; consensus padding: {cons2.consensus_padding}"
    )
