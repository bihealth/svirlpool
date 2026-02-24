# %%
import json
from gzip import open as gzip_open
from pathlib import Path

import cattrs
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus, consensus_class
from svirlpool.localassembly.consensus import (
    get_max_extents_of_read_alignments_on_cr,
    get_read_alignment_intervals_in_region,
)
from svirlpool.util.datatypes import Alignment
from svirlpool.util.util import (
    align_reads_with_minimap,
    dict_to_seqrecord,
    generate_sequence,
    reverse_complement,
    write_sequences_to_fasta,
)

# %%
DATA_DIR = Path(__file__).parent / "data" / "consensus"


def generate_test_data_for_read_trimming_tests() -> None:
    R: list[str] = generate_sequence(2000, seed=0)  # 2 kb random DNA seq
    R0 = R[100:500]
    R1 = R[-500:-100]
    x = reverse_complement(R[850:950])
    r0 = generate_sequence(500, seed=1)
    r1 = generate_sequence(100, seed=2)
    # assemble the read sequence:
    read0 = R0 + r0 + x + r1 + R1
    read1 = reverse_complement(read0)
    path_test_ref = DATA_DIR / "test_trimming_ref.fasta"
    path_test_reads = DATA_DIR / "test_trimming_reads.fasta"
    # What is expected from this data:
    # reads are 1500 bp long
    # cutting on the reference happens at 100 and -100 (1900)
    # which leaves the reads at a total length of 1500 bp
    # and the first 400 bp and the last 400 bp should be aligned
    # break ends should be at 500 and 1500 on the reference, and at 400 and -400 on the reads
    # additional break ends should be at 850 and 950 on the reference, and at 900 and 1000 on the forward read,
    # and at 600 and 700 on the reverse read.
    write_sequences_to_fasta(
        seqs=[R], chrnames=True, path=path_test_ref, prefix="trimming_ref"
    )
    write_sequences_to_fasta(
        seqs=[read0, read1],
        chrnames=False,
        path=path_test_reads,
        prefix="trimming_reads",
    )
    align_reads_with_minimap(
        reference=str(path_test_ref),
        reads=str(path_test_reads),
        bamout=DATA_DIR / "test_trimming_alignments.bam",
        aln_args=" -z100,100 -r100,100",
    )


def generate_test_trimming_alignments_json() -> None:
    """Convert test_trimming_alignments.bam to a serialisable JSON.gz file.

    Reads every alignment from the BAM, wraps it in a datatypes.Alignment
    object (which is serialisable via cattrs), and writes the result to
    ``DATA_DIR / "test_trimming_alignments.json.gz"``.
    Call this once after regenerating the BAM with
    ``generate_test_data_for_read_trimming_tests``.
    """
    bam_path = DATA_DIR / "test_trimming_alignments.bam"
    json_path = DATA_DIR / "test_trimming_alignments.json.gz"
    alignments: list[Alignment] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as f:
        for aln in f.fetch(until_eof=True):
            alignments.append(Alignment.from_pysam(aln))
    with gzip_open(str(json_path), "wt") as f:
        json.dump({"alignments": [aln.unstructure() for aln in alignments]}, f)


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
    with gzip_open(DATA_DIR / "consensus_padding.11.0.json.gz", "rt") as f:
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
        result1 = consensus.create_padding_for_consensus(
            consensus_object=consensus_object1,
            cutreads=cutreads1,
            read_records=read_records1,
        )
    # print all but not the sequences
    expected_1 = consensus_class.ConsensusPadding(
        sequence="",
        readname_left="b4423873-b948-436f-80e8-213f79b01b2a",
        readname_right="f7ada504-4c28-4be8-8bb5-1d17f4888286",
        padding_size_left=39674,
        padding_size_right=36191,
        consensus_interval_on_sequence_with_padding=(39674, 41744),
    )

    assert abs(result1.padding_size_left - expected_1.padding_size_left) < 5, (
        f"Expected left padding {expected_1.padding_size_left}, got {result1.padding_size_left}"
    )
    assert abs(result1.padding_size_right - expected_1.padding_size_right) < 5, (
        f"Expected right padding {expected_1.padding_size_right}, got {result1.padding_size_right}"
    )
    assert (
        abs(
            result1.consensus_interval_on_sequence_with_padding[0]
            - expected_1.consensus_interval_on_sequence_with_padding[0]
        )
        < 5
    ), (
        f"Expected consensus start {expected_1.consensus_interval_on_sequence_with_padding[0]}, got {result1.consensus_interval_on_sequence_with_padding[0]}"
    )
    assert (
        abs(
            result1.consensus_interval_on_sequence_with_padding[1]
            - expected_1.consensus_interval_on_sequence_with_padding[1]
        )
        < 5
    ), (
        f"Expected consensus end {expected_1.consensus_interval_on_sequence_with_padding[1]}, got {result1.consensus_interval_on_sequence_with_padding[1]}"
    )

    with gzip_open(DATA_DIR / "consensus_padding.11.1.json.gz", "rt") as f:
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

    result2 = consensus.create_padding_for_consensus(
        consensus_object=consensus_object2,
        cutreads=cutreads2,
        read_records=read_records2,
    )
    expected_2 = consensus_class.ConsensusPadding(
        sequence="",
        readname_left="800ed304-2134-49fd-9a24-b93f26713816",
        readname_right="4c08354-027f-4cf0-8534-777467668ca1",
        padding_size_left=28150,
        padding_size_right=29199,
        consensus_interval_on_sequence_with_padding=(28150, 30228),
    )
    # print(result2)
    # counter-lengths are
    # assert a rough equivalence of the results
    assert abs(result2.padding_size_left - expected_2.padding_size_left) < 5, (
        f"Expected left padding {expected_2.padding_size_left}, got {result2.padding_size_left}"
    )
    assert abs(result2.padding_size_right - expected_2.padding_size_right) < 5, (
        f"Expected right padding {expected_2.padding_size_right}, got {result2.padding_size_right}"
    )
    assert (
        abs(
            result2.consensus_interval_on_sequence_with_padding[0]
            - expected_2.consensus_interval_on_sequence_with_padding[0]
        )
        < 5
    ), (
        f"Expected consensus start {expected_2.consensus_interval_on_sequence_with_padding[0]}, got {result2.consensus_interval_on_sequence_with_padding[0]}"
    )
    assert (
        abs(
            result2.consensus_interval_on_sequence_with_padding[1]
            - expected_2.consensus_interval_on_sequence_with_padding[1]
        )
        < 5
    ), (
        f"Expected consensus end {expected_2.consensus_interval_on_sequence_with_padding[1]}, got {result2.consensus_interval_on_sequence_with_padding[1]}"
    )


def test_create_padding_for_consensus_forward_breakends() -> None:
    with gzip_open(DATA_DIR / "consensus_padding.forward_breakends.json.gz", "rt") as f:
        data = json.load(f)
        consensus_object: consensus_class.Consensus = cattrs.structure(
            data["consensus_object"], consensus_class.Consensus
        )
        cutreads: dict[str, SeqRecord] = {
            name: dict_to_seqrecord(rec) for name, rec in data["cutreads"].items()
        }
        read_records: dict[str, SeqRecord] = {
            name: dict_to_seqrecord(rec) for name, rec in data["read_records"].items()
        }
        result = consensus.create_padding_for_consensus(
            consensus_object=consensus_object,
            cutreads=cutreads,
            read_records=read_records,
        )

    expected = consensus_class.ConsensusPadding(
        sequence="",
        readname_left="read_left",
        readname_right="read_right",
        padding_size_left=432,
        padding_size_right=916,
        consensus_interval_on_sequence_with_padding=(432, 1708),
    )
    assert abs(result.padding_size_left - expected.padding_size_left) < 5, (
        f"Expected left padding {expected.padding_size_left}, got {result.padding_size_left}"
    )
    assert abs(result.padding_size_right - expected.padding_size_right) < 5, (
        f"Expected right padding {expected.padding_size_right}, got {result.padding_size_right}"
    )
    assert (
        abs(
            result.consensus_interval_on_sequence_with_padding[0]
            - expected.consensus_interval_on_sequence_with_padding[0]
        )
        < 5
    ), (
        f"Expected consensus start {expected.consensus_interval_on_sequence_with_padding[0]}, got {result.consensus_interval_on_sequence_with_padding[0]}"
    )
    assert (
        abs(
            result.consensus_interval_on_sequence_with_padding[1]
            - expected.consensus_interval_on_sequence_with_padding[1]
        )
        < 5
    ), (
        f"Expected consensus end {expected.consensus_interval_on_sequence_with_padding[1]}, got {result.consensus_interval_on_sequence_with_padding[1]}"
    )


# %%


def test_trim_reads_INVDEL() -> None:
    # test consensus.trim_reads
    # Test data: 2 kb reference, two reads (forward and reverse complement).
    # Each read has 3 alignments:
    #   - ref 100-501 (left anchor)
    #   - ref 849-955 (inverted segment)
    #   - ref 1500-1900 (right anchor)
    # Expected: both reads are trimmed to full length (1500 bp, start=0, end=1500)
    # since the outermost reference extent spans ref 100-1900.
    json_path = DATA_DIR / "test_trimming_alignments.json.gz"
    reads_path = DATA_DIR / "test_trimming_reads.fasta"

    # Load all alignments keyed by crID=0 from the serialised JSON.gz file
    # (avoids pysam BAM parsing issues in automated testing)
    with gzip_open(str(json_path), "rt") as f:
        _data = json.load(f)
    _loaded: list[Alignment] = cattrs.structure(_data["alignments"], list[Alignment])
    dict_alignments: dict[int, list[pysam.AlignedSegment]] = {
        0: [aln.to_pysam() for aln in _loaded]
    }

    # Load full read sequences from FASTA (original orientation)
    read_records: dict[str, SeqRecord] = {
        rec.id: rec for rec in SeqIO.parse(reads_path, "fasta")
    }

    # Compute per-read alignment intervals over the full reference (length 2000)
    all_alns = [aln for alns in dict_alignments.values() for aln in alns]
    dict_all_intervals = get_read_alignment_intervals_in_region(
        region_start=0,
        regions_end=2000,
        alignments=all_alns,
        buffer_clipped_length=0,
    )
    # Reduce to max extents: (read_start, read_end, ref_start_chr, ref_start, ref_end_chr, ref_end)
    intervals = get_max_extents_of_read_alignments_on_cr(dict_all_intervals)

    result = consensus.trim_reads(
        dict_alignments=dict_alignments,
        intervals=intervals,
        read_records=read_records,
    )
    forward_read_name = "trimming_reads518f25de"
    reverse_read_name = "trimming_reads6f83015b"

    expected_description_a = "crID=0,start=0,end=1500,ref_start_chr=trimming_ref0,ref_start=100,ref_end_chr=trimming_ref0,ref_end=1900"
    expected_description_b = "crID=0,start=0,end=1500,ref_start_chr=trimming_ref0,ref_start=100,ref_end_chr=trimming_ref0,ref_end=1900"

    assert result[forward_read_name].description == expected_description_a, (
        f"Expected description for forward read: {expected_description_a}, got: {result[forward_read_name].description}"
    )
    assert result[reverse_read_name].description == expected_description_b, (
        f"Expected description for reverse read: {expected_description_b}, got: {result[reverse_read_name].description}"
    )


# %%
