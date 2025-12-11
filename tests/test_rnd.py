import json
from gzip import open as gzip_open
from pathlib import Path

import cattrs
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from svirlpool.localassembly import consensus
from svirlpool.localassembly.consensus import consensus_class
from svirlpool.util.util import (Direction, dict_to_seqrecord,
                                 get_alignments_to_reference,
                                 get_read_position_on_ref)

DATA_DIR = Path("/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/")

# def load_test_data(prefix: str) -> tuple[consensus_class.Consensus, Alignment]:
#     # prefix is like "INV.1"
#     filename_consensus = f"{prefix}.consensus.json.gz"
#     filename_alignment = f"{prefix}.aln.json.gz"
#     consensus_path = DATA_DIR / "consensus_class" / filename_consensus
#     # load the gzipped json dumped structured objects
#     with gzip_open(consensus_path, "rt") as f:
#         consensus_data = f.read()
#         consensus: consensus_class.Consensus = cattrs.structure(
#             json.loads(consensus_data), consensus_class.Consensus
#         )
#     # load the corresponding alignment
#     alignment_path = DATA_DIR / "consensus_class" / filename_alignment
#     with gzip_open(alignment_path, "rt") as f:
#         alignment_data = f.read()
#         alignment: Alignment = cattrs.structure(json.loads(alignment_data), Alignment)
#     return consensus, alignment

def debug_stuff() -> None:
    # load data for files consensus_padding.1.json.gz, consensus_padding.2.json.gz
    with gzip_open(DATA_DIR / "consensus/consensus_padding.1.json.gz", "rt") as f:
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
    with gzip_open(DATA_DIR / "consensus/consensus_padding.2.json.gz", "rt") as f:
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


    read_paddings_for_consensus = consensus._get_all_read_padding_intervals(
        consensus_object=consensus_object1, cutreads=cutreads1, read_records=read_records1
    )
    print("read_paddings_for_consensus:")
    print(read_paddings_for_consensus)
    expected_read_paddings_for_consensus = {
        'a4c08354-027f-4cf0-8534-777467668ca1': (0, 6221, 6221, 30946, True),
        '8782cefd-4a2c-4908-bf28-5d946f3163e7': (0, 7723, 7723, 29683, True),
        'ef0d336d-ea8f-4489-aa8a-c7da86133a36': (10586, 20817, 10231, 20817, False),
        'd4f7ea49-305d-4fa0-9221-428922d60591': (5097, 14342, 9245, 14342, False),
        'aa6212dc-3ebe-4964-a522-8f987863c000': (0, 6609, 6609, 10070, True),
        'af009bf3-8fd5-47ea-8f5e-d71595f72b7d': (0, 8820, 8820, 8820, True),
        'f7ada504-4c28-4be8-8bb5-1d17f4888286': (31691, 41755, 10064, 41755, True),
        'c719dee5-9997-4f73-b933-cc56d1043a43': (0, 1356, 1356, 25377, False),
        '800ed304-2134-49fd-9a24-b93f26713816': (28150, 29051, 901, 29051, True),
        '1fbfd7ac-d2da-478e-9941-eb39a3303870': (12979, 13663, 684, 13663, True)}
    print('\n===============================================\n')

    padding_sizes_per_read = consensus._get_padding_sizes_per_read(
        read_paddings_for_consensus=read_paddings_for_consensus
    )
    print("padding_sizes_per_read:")
    print(padding_sizes_per_read)
    expected_padding_sizes_per_read = {
        'a4c08354-027f-4cf0-8534-777467668ca1': (0, 24725),
        '8782cefd-4a2c-4908-bf28-5d946f3163e7': (0, 21960),
        'ef0d336d-ea8f-4489-aa8a-c7da86133a36': (0, 10586),
        'd4f7ea49-305d-4fa0-9221-428922d60591': (0, 5097),
        'aa6212dc-3ebe-4964-a522-8f987863c000': (0, 3461),
        'af009bf3-8fd5-47ea-8f5e-d71595f72b7d': (0, 0),
        'f7ada504-4c28-4be8-8bb5-1d17f4888286': (31691, 0),
        'c719dee5-9997-4f73-b933-cc56d1043a43': (24021, 0),
        '800ed304-2134-49fd-9a24-b93f26713816': (28150, 0), 
        '1fbfd7ac-d2da-478e-9941-eb39a3303870': (12979, 0)}
    print('\n===============================================\n')

    padding_reads = consensus._get_padding_read_names_of_consensus(
        padding_sizes_per_read=padding_sizes_per_read
    )
    print("padding_reads:")
    print(padding_reads)
    """('f7ada504-4c28-4be8-8bb5-1d17f4888286', 'a4c08354-027f-4cf0-8534-777467668ca1')"""
    print('\n===============================================\n')

    padding_intervals = consensus._get_read_padding_intervals(
        read_paddings_for_consensus=read_paddings_for_consensus,
        padding_read_names_of_consensus=padding_reads,
    )
    print("padding_intervals:")
    print(padding_intervals)
    """((0, 31691), (6221, 30946))"""
    print('\n===============================================\n')
    
    padding = consensus._create_padding_object(
        cons=consensus_object1,
        read_paddings_for_consensus=read_paddings_for_consensus,
        padding_sizes_per_read=padding_sizes_per_read,
        padding_reads=padding_reads,
        padding_intervals=padding_intervals,
        read_records=read_records1,
    )
    print("padding:")
    print(padding)
    expected_padding = consensus_class.ConsensusPadding(
        sequence='...',
        readname_left='f7ada504-4c28-4be8-8bb5-1d17f4888286',
        readname_right='a4c08354-027f-4cf0-8534-777467668ca1',
        padding_size_left=31691,
        padding_size_right=24725,
        consensus_interval_on_sequence_with_padding=(31691, 42777)) # (start , size of consenus - right padding)
    print('\n===============================================\n')

    # to prove the correctness of the padding sizes of the chosen reads,
    # align them to the consenusus sequence and print their cigar strings

    # create a SeqRecord of the consensus sequence
    consensus_seqrecord = SeqRecord(
        seq=Seq(consensus_object1.consensus_sequence),
        id="consensus_with_padding",
        description="consensus sequence with padding"
    )
    reads = [read_records1[padding_reads[0]], read_records1[padding_reads[1]]]

    alignments = get_alignments_to_reference(
        reads = reads,
        reference = consensus_seqrecord,
    )
    for aln in alignments:
        print(f"Read {aln.query_name}, flag: {aln.flag} aligned to consensus with CIGAR: {aln.cigarstring}")
    
    """
    Read f7ada504-4c28-4be8-8bb5-1d17f4888286, flag: 0 aligned to consensus with CIGAR: 31691S ... 10S ### seems correct
    Read f7ada504-4c28-4be8-8bb5-1d17f4888286, flag: 2048 aligned to consensus with CIGAR: 30501H ... 10983H
    Read f7ada504-4c28-4be8-8bb5-1d17f4888286, flag: 2064 aligned to consensus with CIGAR: 37658H ... 3824H

    Read a4c08354-027f-4cf0-8534-777467668ca1, flag: 272 aligned to consensus with CIGAR: 6689S ... 23967S
    Read a4c08354-027f-4cf0-8534-777467668ca1, flag: 2048 aligned to consensus with CIGAR: 23955H ... 6678H
    Read a4c08354-027f-4cf0-8534-777467668ca1, flag: 0 aligned to consensus with CIGAR: 1S ... 24720S ### seems correct
    Read a4c08354-027f-4cf0-8534-777467668ca1, flag: 2064 aligned to consensus with CIGAR: 4964H ... 25686H
    """
    
    # test 11.0
    # cons1, aln1 = load_test_data("INV.1")
    # result = consensus_class.get_consensus_core_alignment_interval_on_reference(
    #     consensus=cons1,
    #     alignment=aln1.to_pysam(),
    # )
    # expected = ("6", 169093104, 169093633)

def test_get_read_position_on_ref_inversion11() -> None:
    alns = DATA_DIR / "util" / "inv.11.bam"
    pysam_alignments = list(pysam.AlignmentFile(alns).fetch())
    #for each alignment, query regular positions on the read
    query_positions = list(range(0,70_000, 500))
    point_result_pairs:dict[int, list[tuple[int,int]]] = {}
    # for each alignment, get the reference position for each query position
    for i,aln in enumerate(pysam_alignments):
        point_result_pairs[i] = [] # collect results for this alignment
        for pos in query_positions:
            pos_ref = get_read_position_on_ref(
                alignment=aln,
                position=pos,
                direction=Direction.NONE
            )
            point_result_pairs[i].append((pos, pos_ref))
    # now compare the results. Print each position and the resulting reference positions
    for i in range(len(query_positions)):
        pos = query_positions[i]
        refs = [point_result_pairs[aln_idx][i][1] for aln_idx in range(len(pysam_alignments))]
        print(f"{pos}: {refs}")
    # additional output:
    # for each alignment, create a bed file of the reference position (+1) and the pos (key) value
    # for that, sort the point_result_pairs by reference position for each alignment
    for i, aln in enumerate(pysam_alignments):
        output_base = DATA_DIR / "util" / f"inv.11.aln{i}.readpos.bed"
        with open(output_base, "w") as bedf:
            sorted_positions = sorted(point_result_pairs[i], key=lambda x: x[1] if x[1] is not None else -1)
            for read_pos, ref_pos in sorted_positions:
                if ref_pos is not None:
                    bedf.write(f"{aln.reference_name}\t{ref_pos}\t{ref_pos+1}\t{read_pos}\n")
    # for each alignment, print the inferred read length
    for i, aln in enumerate(pysam_alignments):
        print(f"Alignment {i} inferred read length: {aln.infer_read_length()}")
    
test_get_read_position_on_ref_inversion11()
