import json
from gzip import open as gzip_open
from pathlib import Path

import cattrs
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus, consensus_class
from svirlpool.util.util import dict_to_seqrecord

DATA_DIR = Path(__file__).parent / "data" / "consensus"

# insterted into create_padding_for_consensus
# This is how I saved the data, then gzipped it:
# DEBUG: write input to a json file to test it.
# consensus_object can easily be serialized with .unsutructure() and then json dumped
# SeqRecord objects are serialized with util.seqrecord_to_dict()
# if consensus_object.ID in ["11.1", "11.0"]:
#     # remove letter annotations from sequences in cutreads and read_records for better readability
#     for cutread in cutreads.values():
#         cutread.letter_annotations = {}
#     for read in read_records.values():
#         read.letter_annotations = {}
#     iddict = {"11.1":2, "11.0":1}
#     output_path = f"/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus/consensus_padding.{iddict[consensus_object.ID]}.json"
#     with open(output_path, "w") as f:
#         json.dump(
#             {
#                 "consensus_object": consensus_object.unstructure(),
#                 "cutreads": {
#                     name: util.seqrecord_to_dict(cutread)
#                     for name, cutread in cutreads.items()
#                 },
#                 "read_records": {
#                     name: util.seqrecord_to_dict(read)
#                     for name, read in read_records.items()
#                 },
#             },
#             f,
#             indent=4,
#         )


def test_create_padding_for_consensus() -> None:
    # load data for files consensus_padding.1.json.gz, consensus_padding.2.json.gz
    with gzip_open(DATA_DIR / "consensus_padding.1.json.gz", "rt") as f:
        data1 = json.load(f)
        consensus_object1: consensus_class.Consensus = cattrs.structure(
            data1["consensus_object"], consensus_class.Consensus
        )
        cutreads1: dict[str, SeqRecord] = {
            name: dict_to_seqrecord(rec) for name, rec in data1["cutreads"].items()
        }
        read_records1: dict[str, SeqRecord] = {
            name: dict_to_seqrecord(rec) for name, rec in data1["read_records"].items()
        }
    with gzip_open(DATA_DIR / "consensus_padding.2.json.gz", "rt") as f:
        data2 = json.load(f)
        consensus_object2: consensus_class.Consensus = cattrs.structure(
            data2["consensus_object"], consensus_class.Consensus
        )
        cutreads2: dict[str, SeqRecord] = {
            name: dict_to_seqrecord(rec) for name, rec in data2["cutreads"].items()
        }
        read_records2: dict[str, SeqRecord] = {
            name: dict_to_seqrecord(rec) for name, rec in data2["read_records"].items()
        }
    # generate results
    result1 = consensus.create_padding_for_consensus(
        consensus_object=consensus_object1,
        cutreads=cutreads1,
        read_records=read_records1,
    )
    expected_1 = consensus_class.ConsensusPadding(
        sequence="",
        readname_left="f7ada504-4c28-4be8-8bb5-1d17f4888286",
        readname_right="a4c08354-027f-4cf0-8534-777467668ca1",
        padding_size_left=31691,
        padding_size_right=24725,
        consensus_interval_on_sequence_with_padding=(31691, 42777),
    )
    assert result1 == expected_1, f"Expected {expected_1}, got {result1}"

    result2 = consensus.create_padding_for_consensus(
        consensus_object=consensus_object2,
        cutreads=cutreads2,
        read_records=read_records2,
    )
    expected_2 = consensus_class.ConsensusPadding(
        sequence="",
        readname_left="b4423873-b948-436f-80e8-213f79b01b2a",
        readname_right="6707faed-82b1-4351-b153-657b4c12bbfe",
        padding_size_left=35204,
        padding_size_right=20184,
        consensus_interval_on_sequence_with_padding=(35204, 46326),
    )
    assert result2 == expected_2, f"Expected {expected_2}, got {result2}"
