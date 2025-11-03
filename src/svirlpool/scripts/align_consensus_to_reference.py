# script to align consensus sequences to reference genome and
# additionally align the clipped tials to the reference genome
# %%
import argparse
import subprocess
import tempfile
import typing
from pathlib import Path
from shlex import split

import numpy as np
import pysam
from logzero import logger as log
from tqdm import tqdm

from . import util

# %%
# # load test alignments
# alns = Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d01/HG002.consensus.bam")
# # load alignments to list
# alignments = list(pysam.AlignmentFile(alns, "rb"))
# a=alignments[44]
# reference_positions = [7588766,7589107]


# %%
# =============================================================================
# cut alignments at given reference positions
# returns a list of alignedSegments with properly cut cigar strings and other data.
# 1) find cutting positions on read and reference
# 2) cut cigar strings, sequences and qualities
# 3) create new alignedSegments
# 4) return tuple(list of final reference_positions, list of new alignedSegments)
def cut_alignment(
    a: pysam.AlignedSegment, reference_positions: typing.List[int]
) -> typing.Tuple[typing.List[int], typing.List[pysam.AlignedSegment]]:
    # find cutting positions on read and reference
    read_positions = [
        util.get_ref_position_on_read(
            alignment=a, direction=util.Direction.LEFT, position=p
        )
        for p in reference_positions
    ]
    # then transfer back the read positions to reference positions to ensure full transitiveness
    ref_positions = [
        util.get_read_position_on_ref(
            alignment=a, direction=util.Direction.LEFT, position=p
        )
        for p in list(zip(*read_positions))
    ]
    read_positions = [
        util.get_ref_position_on_read(
            alignment=a, direction=util.Direction.RIGHT, position=p
        )
        for p in list(zip(*ref_positions))
    ]
    # now cut the cigar strings, sequences and qualities
    # keep all other information from the original alignment in each segment
    segments = []
    # get the first segment
    new_segment: pysam.AlignedSegment = a.copy()
    new_segment.query_sequence = a.query_sequence[: read_positions[0]]
    new_segment.query_qualities = a.query_qualities[: read_positions[0]]
    new_segment.cigartuples = a.cigartuples[: ref_positions[0]]
    new_segment.reference_end = a.reference_start + ref_positions[0]
    new_segment.query_name = f"{a.query_name}_0"
    segments.append(new_segment)
    # get the middle segments
    for i in range(1, len(ref_positions)):
        new_segment: pysam.AlignedSegment = a.copy()
        new_segment.query_sequence = a.query_sequence[
            read_positions[i - 1] : read_positions[i]
        ]
        new_segment.query_qualities = a.query_qualities[
            read_positions[i - 1] : read_positions[i]
        ]
        new_segment.cigartuples = a.cigartuples[ref_positions[i - 1] : ref_positions[i]]


# get_reference_sequence_from_read_alignment
# %%


def split_bam(path_consensus_to_ref_alignments: Path) -> None:
    # check if unaligned reads exist
    last_alignment = list(pysam.AlignmentFile(path_consensus_to_ref_alignments, "rb"))[
        -1
    ]
    if not last_alignment.is_unmapped:
        # can return if all alignments are mapped
        return None
    # split bam file by copying all unaligned reads to a new bam file and removing the unaligned from path_consensus_to_ref_alignments
    cmd_unaligned = f"samtools view -b -f 4 {path_consensus_to_ref_alignments}"
    path_unaligned = str(path_consensus_to_ref_alignments.with_suffix(".unaligned.bam"))
    subprocess.check_call(split(cmd_unaligned), stdout=open(path_unaligned, "wb"))
    tmp_aligned = tempfile.NamedTemporaryFile(suffix=".bam", delete=True)
    cmd_aligned = f"samtools view -b -F 4 {path_consensus_to_ref_alignments}"
    subprocess.check_call(split(cmd_aligned), stdout=open(tmp_aligned.name, "wb"))
    cmd_mv = f"cp {tmp_aligned.name} {path_consensus_to_ref_alignments}"
    subprocess.check_call(split(cmd_mv))
    cmd_index_aligned = f"samtools index {path_consensus_to_ref_alignments}"
    subprocess.check_call(split(cmd_index_aligned))
    cmd_index_unaligned = f"samtools index {path_unaligned}"
    subprocess.check_call(split(cmd_index_unaligned))
    log.info(
        f"split {path_consensus_to_ref_alignments} into {path_consensus_to_ref_alignments} and {path_unaligned}"
    )
    return path_unaligned


def align_consensus_sequences_to_reference(
    path_reference: Path,
    path_consensus: Path,
    bamout: Path,
    threads: int,
    aln_args: str,
    tech: str,
) -> None:
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_bamout = Path(tmp_dir.name) / bamout.name
    log.info("aligning consensus sequences to reference with minimap2")
    util.align_reads_with_minimap(
        reference=path_reference,
        bamout=tmp_bamout,
        reads=path_consensus,
        tech=tech,
        threads=threads,
        aln_args=aln_args,
    )
    # refine alignments by aligning the clipped tails of the consensus sequences to the reference
    util.add_clipped_tails_to_alignments(
        alignments=tmp_bamout,
        reference=path_reference,
        threads=threads,
        tech=tech,
        output=bamout,
        aln_args=aln_args,
        min_clipped_length=300,
    )
    tmp_dir.cleanup()


def align_consensus_to_ref_new(
    path_reference: Path,
    path_consensus: Path,
    path_consensus_to_ref_alignments: Path,
    aln_args: str,
    threads: int,
) -> None:
    # first, create 5 different parameter sets for minimap2.
    # the idea is to increase the gap open costs and reduce the mismatch penalty incrementally.
    # asm20 starts with -O6,26
    # first, create a temp directory to store the minimap2 index and the alignment files
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_dir_path = Path(tmp_dir.name)
    # if path_reference is a fasta file, create a minimap2 index
    if path_reference.suffix in [".fa", ".fasta"]:
        log.info(f"creating minimap2 index for {path_reference}..")
        minimap2_index = tmp_dir_path / "reference.mmi"
        cmd_index = f"minimap2 -d {minimap2_index} -k19 -w19 {path_reference}"
        subprocess.check_call(split(cmd_index))
        path_reference = minimap2_index
    # args k,w,U already set in index creation
    minimap2_params = [
        "-U50,500 --rmq -r1k,100k -g10k -A1 -B3 -O16,41 -E2,1 -s200 -z200 -N50",
        "-U50,500 --rmq -r1k,100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50",
        "-U50,500 --rmq -r1k,100k -g10k -A1 -B9 -O16,41 -E2,1 -s200 -z200 -N50",
        "-U50,500 --rmq -r1k,100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50",
    ]
    for i, params in enumerate(minimap2_params):
        path_tmp_aln = tmp_dir_path / f"aln_{i}.bam"
        align_consensus_sequences_to_reference(
            path_reference=path_reference,
            path_consensus=path_consensus,
            bamout=path_tmp_aln,
            aln_args=f"{str(aln_args)} {params} --secondary=no",
            tech="map-ont",
            threads=threads,
        )
        split_bam(path_consensus_to_ref_alignments=path_tmp_aln)
    # Now, it has to be decided which alignments are to be picked.
    # all alignments can be loaded and then given some heuristic, it can be decided which alignments are to be picked.
    # the heuristics should consider the number of gaps, the number of mismatches, the aligned length and the alignment score.
    tmp_unsorted_sam = tmp_dir_path / "unsorted.bam"
    log.info(
        f"pick best alignments from {len(minimap2_params)} alignment files and write to {tmp_unsorted_sam}"
    )
    pick_best_alignments(
        alignment_files=[
            tmp_dir_path / f"aln_{i}.bam" for i in range(len(minimap2_params))
        ],
        bam_output=tmp_unsorted_sam,
    )
    # sort and index the final bam file and write to path_consensus_to_ref_alignments
    log.info(
        f"sort and index {tmp_unsorted_sam} and write to {path_consensus_to_ref_alignments}"
    )
    cmd_sort = f"samtools sort -@ {threads} -O BAM -o {path_consensus_to_ref_alignments} {tmp_unsorted_sam}"
    cmd_index = f"samtools index {path_consensus_to_ref_alignments}"
    subprocess.check_call(split(cmd_sort))
    subprocess.check_call(split(cmd_index))


def align_consensus_to_ref(
    path_reference: Path,
    path_consensus: Path,
    path_consensus_to_ref_alignments: Path,
    aln_args: str,
    threads: int,
) -> None:
    align_consensus_sequences_to_reference(
        path_reference=path_reference,
        path_consensus=path_consensus,
        bamout=path_consensus_to_ref_alignments,
        aln_args=str(aln_args) + " --secondary=no",
        tech="asm10",
        threads=threads,
    )
    log.info("split unaligned reads from consensus to ref alignments")
    split_bam(path_consensus_to_ref_alignments=path_consensus_to_ref_alignments)


# re-write this. align_consensus_to_ref should apply five different parameter sets.
# a minimap2 index should be created for the reference genome and used in all 5 stepts to reduce redundancy.
# map the first parameter set with minimap2 and dump the minimap2 index into a file.
# then use this index for the next 4 mappings.


def pick_best_alignments(alignment_files: typing.List[Path], bam_output: Path) -> None:
    # load alignments from all files in alignment_files to a dict {i:{query_name:[alignments]}}
    # pick for each query_name the best alignments from the given alignment files
    # write each picked alignment to the bam_output
    log.info(f"load and sort alignments from {len(alignment_files)} alignment files")
    all_alignments = dict()
    for i, alns in enumerate(alignment_files):
        dict_alns = dict()
        for a in tqdm(pysam.AlignmentFile(alns, "rb")):
            if a.query_name in dict_alns:
                dict_alns[a.query_name].append(a)
            else:
                dict_alns[a.query_name] = [a]
        all_alignments[i] = dict_alns
    # create set of all query_names in all alignments
    query_names = set(
        [a for alns_dict in all_alignments.values() for a in alns_dict.keys()]
    )
    # first write the header to the bam_output. The header can be taken from the first alignment file
    with pysam.AlignmentFile(alignment_files[0], "rb") as f:
        header = f.header
    with pysam.AlignmentFile(bam_output, "wb", header=header) as f:
        for qn in query_names:
            # compute all scores of all alignments that are found for the query_name
            # and pick the alignments that have the highest score.
            # write these alignments to the bam_output
            alns = [
                aln_dict[qn] for aln_dict in all_alignments.values() if qn in aln_dict
            ]
            scores = [score_of_consensus_alignments(a) for a in alns]
            best_alns = alns[np.argmax(scores)]
            for a in best_alns:
                f.write(a)


# score_of_consensus_alignments computes a score based on the number of gaps, the number of mismatches,
# and the sum of aligned regions and the number of fragments.
def score_of_consensus_alignments(
    alns: typing.List[pysam.AlignedSegment], min_mapq: int = 40
) -> float:
    # for each alignment, get the stats (matches,mismatches,n_gaps)
    stats = np.array(
        [
            [get_bp_matches(aln), get_bp_mismatches(aln), sum(get_n_gaps(aln))]
            for aln in alns
            if aln.mapping_quality >= min_mapq
        ],
        dtype=int,
    )
    if len(stats) == 0:
        return 0
    n_fragments = len(stats)
    sum_matches = np.sum(stats[:, 0])
    sum_mismatches = np.sum(stats[:, 1])
    sum_gaps = np.sum(stats[:, 2])
    # score very bad if mismatches are more than 2% of the aligned bases
    #   the function to get a score on the mismatches is a very steep function that starts at 2% of the aligned bases.
    #   Any value below is a score of 1.
    score_mismatches = (
        1
        if sum_mismatches / sum_matches < 0.02
        else 1 / (1 + np.exp(100 * (sum_mismatches / sum_matches - 0.02)))
    )
    # score is very good if the number of gaps is very low
    score_gaps = 1 / 1 if sum_gaps <= 1 else 1 / (2 ** (sum_gaps - 1))
    # score is very good if the number of fragments is very low
    score_fragments = 1 / 1 if n_fragments <= 1 else 1 / (2 ** (n_fragments - 1))
    return score_mismatches * score_gaps * score_fragments * sum_matches


def get_bp_matches(aln: pysam.AlignedSegment) -> int:
    return aln.get_cigar_stats()[0][0]


def get_bp_mismatches(aln: pysam.AlignedSegment) -> int:
    # mismatches + deletions + insertions = aln.get_cigar_stats()[0][10]
    bp_deletions = aln.get_cigar_stats()[0][2]
    bp_insertions = aln.get_cigar_stats()[0][1]
    return aln.get_cigar_stats()[0][10] - bp_deletions - bp_insertions


def get_n_gaps(aln: pysam.AlignedSegment, min_gap_size=8) -> typing.Tuple[int, int]:
    t, x = (np.array(list(x)) for x in zip(*aln.cigartuples))
    # the number of gaps is the number of 1s or 2s in t, where the corresponding element in x is larger than min_gap_size
    # first get the number of insertions that are larger than min_gap_size bp
    n_ins = np.sum((t == 1) & (x >= min_gap_size))
    # then get the number of deletions that are larger than min_gap_size bp
    n_del = np.sum((t == 2) & (x >= min_gap_size))
    return n_ins, n_del


def run(args, **kwargs):
    align_consensus_to_ref(
        path_reference=args.path_reference,
        path_consensus=args.path_consensus,
        path_consensus_to_ref_alignments=args.output,
        aln_args=args.aln_args,
        threads=args.threads,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="this script aligns consensus sequences to the reference genome and additionally aligns the clipped tails of the consensus sequences to the reference genome"
    )
    parser.add_argument(
        "-r",
        "--path-reference",
        type=Path,
        required=True,
        help="Input path to reference fasta file.",
    )
    parser.add_argument(
        "-i",
        "--path-consensus",
        type=Path,
        required=True,
        help="Input path to fasta file of concatenated consensus sequences.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output path to bam file of consensus-to-reference alignments.",
    )
    parser.add_argument(
        "-x",
        "--aln-args",
        type=str,
        required=False,
        default="",
        help="Arguments for minimap2 of the final consensus-to-reference alignments. E.g. -O4,42.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Number of threads to use.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
