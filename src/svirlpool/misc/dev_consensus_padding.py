# %% experimental script to test consensus sequence elongation and trimming

from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..scripts import (
    consensus,
    consensus_class,
    datatypes,
    util,
)

# %% functions

# def get_reads_of_consensus(consensus:consensus_class.Consensus) -> list[str]:
#     return [readname for _,__,readname in consensus.intervals_cutread_alignments]

# def aln_overlaps_candidate_regions(aln:pysam.AlignedSegment, candidate_regions:list[datatypes.CandidateRegion]) -> bool:
#     # check if the alignment overlaps with any of the candidate regions
#     it_aln = Interval(begin=aln.reference_start, end=aln.reference_end)
#     its_crs = {cr.crID:Interval(begin=cr.referenceStart, end=cr.referenceEnd) for cr in candidate_regions}
#     for cr in candidate_regions:
#         if aln.reference_name != cr.chr:
#             return False
#         if it_aln.overlaps(its_crs[cr.crID]):
#             return True
#     return False

# def get_padding_alignments_of_consensus(
#         consensus:consensus_class.Consensus,
#         alignments:list[pysam.AlignedSegment],
#         candidate_regions:list[datatypes.CandidateRegion]) -> tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
#     # find all candidate regions whose crIDs are in the consensus crIDs
#     selected_crs = [cr for cr in candidate_regions if cr.crID in consensus.crIDs]
#     # filter alignments for the reads that are in the candidate regions and overlap with at least one of the candidate regions
#     reads_of_consensus = get_reads_of_consensus(consensus)
#     filtered_alns = [aln for aln in alignments if aln.query_name in reads_of_consensus and aln_overlaps_candidate_regions(aln, selected_crs)]
#     # find the leftmost alignment
#     leftmost_aln = min(filtered_alns, key=lambda aln: aln.reference_start)
#     # find the rightmost alignment
#     rightmost_aln = max(filtered_alns, key=lambda aln: aln.reference_end)
#     return leftmost_aln, rightmost_aln # not perfect given muli-breaks consensus seqs, but good for now.


def get_read_paddings_for_consensus(
    consensus_object: consensus_class.Consensus, read_records: dict[str, SeqRecord]
) -> dict[str, tuple[int, int, int, int, bool]]:
    """Generates a dict[readname, (min start, max end, cutread length, read length, is_forward)] for all reads in the consensus object."""
    directions: dict[str, bool] = {
        rc.alignment.readname: rc.alignment.to_pysam().is_forward
        for rc in consensus_object.reconstructible_reads
    }

    data: dict[str, tuple[int, int]] = {}
    for rc in consensus_object.reconstructible_reads:
        cutread_description_dict = consensus_class.parse_description(rc.description)
        if rc.alignment.readname not in data:
            data[rc.alignment.readname] = (
                cutread_description_dict["start"],
                cutread_description_dict["end"],
            )
        else:
            data[rc.alignment.readname] = (
                min(data[rc.alignment.readname][0], cutread_description_dict["start"]),
                max(data[rc.alignment.readname][1], cutread_description_dict["end"]),
            )
    result: dict[str, tuple[int, int, int, int]] = {}
    # now add the read length to each entry
    for readname, (start, end) in data.items():
        if readname in read_records:
            read_length = len(read_records[readname].seq)
            result[readname] = (
                start,
                end,
                end - start,
                read_length,
                directions[readname],
            )
    return result


def get_padding_sizes_per_read(
    read_paddings_for_consensus: dict[str, tuple[int, int, int, int]],
) -> dict[str, tuple[int, int]]:
    """Generates a dict[readname, (left padding size, right padding size)] for all reads in the consensus object."""
    """A padding is the sequence that overshoots the cutread alignment on the left and right side (in the orientation of the consensus)."""
    result: dict[str, tuple[int, int]] = {}
    for readname, (
        start,
        end,
        cutread_length,
        read_length,
        forward,
    ) in read_paddings_for_consensus.items():
        alpha: int = start if forward else read_length - end
        beta: int = read_length - end if forward else start
        result[readname] = (alpha, beta)
    return result


def get_padding_read_names_of_consensus(
    padding_sizes_per_read: dict[str, tuple[int, int]],
) -> tuple[str, str]:
    """Returns the read names of the left and right padding reads."""
    return (
        max(padding_sizes_per_read.items(), key=lambda x: x[1][0])[0],
        max(padding_sizes_per_read.items(), key=lambda x: x[1][1])[0],
    )


def get_read_padding_intervals(
    read_paddings_for_consensus: dict[str, tuple[int, int]],
    padding_read_names_of_consensus: tuple[str, str],
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Returns the left and right intervals on the two padding reads, ignorant toward clipped bases."""
    left_read, right_read = padding_read_names_of_consensus
    start, end, cutread_length, read_length, forward = read_paddings_for_consensus[
        left_read
    ]
    left_padding = (0, start) if forward else (end, read_length)
    start, end, cutread_length, read_length, forward = read_paddings_for_consensus[
        right_read
    ]
    right_padding = (end, read_length) if forward else (0, start)
    return left_padding, right_padding


def create_padding(
    cons: consensus_class.Consensus,
    read_paddings_for_consensus: dict[str, tuple[int, int, int, int]],
    padding_sizes_per_read: dict[str, tuple[int, int]],
    padding_reads: tuple[str, str],
    padding_intervals: tuple[tuple[int, int], tuple[int, int]],
    read_records: dict[str, SeqRecord],
) -> datatypes.ConsensusPadding:
    """Creates a new ConsensusPadding object with the padding sequence and the read names."""

    left_padding: Seq = read_records[padding_reads[0]].seq[
        padding_intervals[0][0] : padding_intervals[0][1]
    ]
    # needs to be reverse complimented if the cutread alignment is reverse
    if read_paddings_for_consensus[padding_reads[0]][4] == False:
        left_padding = left_padding.reverse_complement()

    right_padding: Seq = read_records[padding_reads[1]].seq[
        padding_intervals[1][0] : padding_intervals[1][1]
    ]
    # needs to be reverse complimented if the cutread alignment is reverse
    if read_paddings_for_consensus[padding_reads[1]][4] == False:
        right_padding = right_padding.reverse_complement()
    # add padding to the consensus sequence
    padded_consensus_sequence = Seq(
        left_padding + cons.consensus_sequence + right_padding
    )

    new_padding = datatypes.ConsensusPadding(
        sequence=str(padded_consensus_sequence),
        readname_left=padding_reads[0],
        readname_right=padding_reads[1],
        padding_size_left=padding_sizes_per_read[padding_reads[0]][0],
        padding_size_right=padding_sizes_per_read[padding_reads[1]][1],
        consensus_interval_on_sequence_with_padding=(
            padding_sizes_per_read[padding_reads[0]][0],
            padding_sizes_per_read[padding_reads[0]][0] + len(cons.consensus_sequence),
        ),
    )
    return new_padding


def add_padding_to_consensus(
    consensus_object: consensus_class.Consensus, read_records: dict[str, SeqRecord]
) -> datatypes.ConsensusPadding:
    """Creates a ConsensusPadding object with the padding sequence, the consensus interval on the padded squence, and the read names."""

    read_paddings_for_consensus = get_read_paddings_for_consensus(
        consensus_object=consensus_object, read_records=read_records
    )

    padding_sizes_per_read = get_padding_sizes_per_read(
        read_paddings_for_consensus=read_paddings_for_consensus
    )

    padding_reads = get_padding_read_names_of_consensus(
        padding_sizes_per_read=padding_sizes_per_read
    )

    padding_intervals = get_read_padding_intervals(
        read_paddings_for_consensus=read_paddings_for_consensus,
        padding_read_names_of_consensus=padding_reads,
    )
    padding = create_padding(
        cons=consensus_object,
        read_paddings_for_consensus=read_paddings_for_consensus,
        padding_sizes_per_read=padding_sizes_per_read,
        padding_reads=padding_reads,
        padding_intervals=padding_intervals,
        read_records=read_records,
    )
    return padding


# %%

path_alignments = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/platinum/test/region.bam"
)
output_cut_reads = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/platinum/test/region_cutreads.fasta"
)
output_cutread_alns = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/platinum/test/region_cutreads.bam"
)
reference_index = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/chm13v2/chm13v2.0.mmi"
)
consensus_result = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/validation/platinum/svirlpool/chm13/NA12877.30x/consensus/20/consensus.4084.txt"
)

# %% create artificial candidate region

crID = 4084

# get all readnames from the alignments and create a dummy list of sv signals for the candidate region so that the reads get parsed
alignments = list(pysam.AlignmentFile(path_alignments).fetch())
sv_signals = [
    datatypes.ExtendedSVsignal(
        readname=aln.query_name,
        chr="chr3",
        chrID=2,
        coverage=29,
        samplename="dummy",
        forward=aln.is_forward,
        repeatID=-1,
        strength=1.0,
        ref_start=aln.reference_start,
        ref_end=aln.reference_end,
        read_start=aln.query_alignment_start,
        read_end=aln.query_alignment_end,
        size=100,
        sv_type=0,
    )
    for aln in alignments
    if aln.is_secondary == False and aln.is_supplementary == False
]

crs_dict = {
    0: datatypes.CandidateRegion(
        crID=crID,
        chr="chr3",
        referenceID="chm13v2_0",
        referenceStart=90515304,
        referenceEnd=90544881,
        sv_signals=sv_signals,
    )
}

# %%
alns, alns_wt = consensus.get_read_alignments_for_crs(
    crs=list(crs_dict.values()), alignments=path_alignments
)
dict_all_intervals = consensus.get_read_alignment_intervals_in_cr(
    crs=list(crs_dict.values()), dict_alignments=alns, buffer_clipped_length=4000
)
max_intervals = consensus.get_max_extents_of_read_alignments_on_cr(
    dict_all_intervals=dict_all_intervals
)
read_records: dict[str, SeqRecord] = consensus.get_full_read_sequences_of_alignments(
    dict_alignments=alns, path_alignments=path_alignments
)
cutreads = consensus.trim_reads(
    dict_alignments=alns, intervals=max_intervals, read_records=read_records
)
# %% write cutreads to output_cut_reads and align to reference index
# with open(output_cut_reads, "w") as output_handle:
#     for readname, cutread in cutreads.items():
#         SeqIO.write(cutread, output_handle, "fasta")
# util.align_reads_with_minimap(reference=reference_index,reads=output_cut_reads,bamout=output_cutread_alns,threads=8)
# %% load consensus sequence
# consensus_results = next(consensuses_to_svprimitives.parse_crs_container_results(consensus_result))
# consensus_IDs = list(consensus_results.consensus_dicts.keys())
# %% select reads to use fro padding

# each consensus_ID is a consensus sequence. it has reads assigned to it. Access them with
# e.g. consensus_results.consensus_dicts['4084.0'].intervals_cutread_alignments[0][2]

# cons = consensus_results.consensus_dicts[f'{str(crID)}.0']
# alignments_flat = [aln for crID,A in alns.items() for aln in A]

# aln_l, aln_r = get_padding_alignments_of_consensus(
#         consensus=cons,
#         alignments=alignments_flat,
#         candidate_regions=list(crs_dict.values()))

# %%
# test - find aln of bf39db65-4001-4538-80c8-eafe4cea0113
# aln_l = next((aln for aln in alignments_flat if aln.query_name == "37014ece-55a4-4df7-9b26-68b57c337497"), None) # has 5 clipped bases right?

# %%

# # get the cutread alignments of the same reads as in aln_l and aln_r
# # find them in consensus.reconstructible_reads (a list of ReconstructibleSequence objects)
# # the query name is present in the alignment.readname member
# cutread_l = next((cutread for cutread in cons.reconstructible_reads if cutread.alignment.readname == aln_l.query_name), None)
# cutread_r = next((cutread for cutread in cons.reconstructible_reads if cutread.alignment.readname == aln_r.query_name), None)
# # parse the string to a dict. value,key pair are separated by comma. a pair is separated by an equal sign
# cutread_l_description_dict = consensus.parse_description(cutread_l.description)
# cutread_r_description_dict = consensus.parse_description(cutread_r.description)
# cutread_l_pysam = cutread_l.alignment.to_pysam()
# cutread_r_pysam = cutread_r.alignment.to_pysam()
# cutread_l_length = cutread_l_pysam.infer_read_length()
# cutread_r_length = cutread_r_pysam.infer_read_length()

# util.query_start_end_on_read(cutread_l_pysam)
# util.query_start_end_on_read(cutread_r_pysam)

# # translate the cut-read start and end positions on the consensus to the read start and end positions
# cutread_l_description_dict['end'] - cutread_l_description_dict['start'] # cut positions on read
# # trace back cutread_l_description_dict['start'] - min(util.query_start_end_on_read(cutread_l_pysam))
# cutread_l_pysam.cigartuples[0] # are now at the start


# %%
# to test, align the padded consensus sequence to the reference
path_test_consensus_padded = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/platinum/test/padded_consensus.fasta"
)
path_test_alignment = Path(
    "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/platinum/test/padded_consensus.bam"
)
with open(path_test_consensus_padded, "w") as output_handle:
    SeqIO.write(
        SeqRecord(seq=padded_consensus_sequence, id=f"padded.{crID}"),
        output_handle,
        "fasta",
    )
util.align_reads_with_minimap(
    reference=reference_index,
    reads=path_test_consensus_padded,
    bamout=path_test_alignment,
    threads=4,
)

# %%

# combine all to a function that takes arguments: consensus_object,
