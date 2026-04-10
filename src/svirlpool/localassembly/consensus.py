"""
Consensus Assembly Module

This script creates consensus assemblies from sets of sequencing reads.
It provides functionality for alignment processing, read clustering,
consensus generation, and scoring.
"""

# =============================================================================
# IMPORTS AND CONFIGURATION
# =============================================================================

import argparse
import hashlib
import json
import logging
import os
import random
import shlex
import shutil
import sqlite3
import subprocess
import tempfile
from pathlib import Path

import attrs
import cattrs
import matplotlib
import numpy as np
import pysam
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree
from sklearn.cluster import KMeans, SpectralClustering

from ..signalprocessing import alignments_to_rafs, copynumber_tracks
from ..util import datatypes, util
from ..util.signal_loss_logger import get_signal_loss_logger
from . import consensus_class, consensus_lib

matplotlib.use("Agg")  # Use non-interactive backend
logging.getLogger("matplotlib").setLevel(logging.WARNING)
log = logging.getLogger(__name__)


# =============================================================================
# FILE I/O AND BAM/SAM PROCESSING
# =============================================================================


def subsample_alignments(
    input_bamfile: Path, output_samfile: Path, number: int = 20, verbose: bool = False
) -> list[str]:
    """Filter BAM file to include only the longest reads up to 'number' and write to SAM format."""
    # Read all alignments from bam file using context manager
    with pysam.AlignmentFile(input_bamfile, "rb") as infile:
        alignments = list(infile)

        alignments_filtered = [
            a for a in alignments if a.query_name != a.reference_name
        ]
        alignments_by_length = sorted(
            alignments_filtered,
            key=lambda a: a.reference_end - a.reference_start,
            reverse=True,
        )
        alignments_sampled = alignments_by_length[
            : min(number, len(alignments_by_length))
        ]

        if len(alignments_by_length) == 0:
            # No valid alignments, create empty output file
            with pysam.AlignmentFile(output_samfile, "w", template=infile) as f:
                pass
            return []

        # Get readnames of selected alignments
        readnames = [a.query_name for a in alignments_sampled]

        # Write selected alignments to samfile
        with pysam.AlignmentFile(output_samfile, "w", template=infile) as f:
            for a in alignments_sampled:
                f.write(a)

    return readnames


def filter_sam_by_readnames(
    input_samfile: Path, output_samfile: Path, readnames: set[str]
) -> None:
    """Filter SAM file to include only alignments with readnames in the provided list."""
    with (
        pysam.AlignmentFile(input_samfile, "r") as infile,
        pysam.AlignmentFile(output_samfile, "w", template=infile) as outfile,
    ):
        for aln in infile:
            if aln.query_name in readnames:
                outfile.write(aln)


# =============================================================================
# READ SEQUENCE EXTRACTION AND PROCESSING
# =============================================================================


def get_all_supplementary_positions_with_full_read_sequence(
    aln: pysam.AlignedSegment,
) -> list[tuple[str, int]]:
    """
    Get the next supplementary alignment with a full read sequence.
    Only apply this to alignments with hard clippings!
    Returns None if no such alignment is found.
    """
    if not aln.has_tag("SA"):
        raise ValueError(
            f"Alignment {aln.query_name} in {aln.reference_name}:{aln.reference_start}-{aln.reference_end} \
has no supplementary alignments. This should not happen. For each aligned read should exist at least one alignment \
that brings the whole read sequence with it - identified by absence of 'H' or hard clippings in the cigar string. \
This is only tested with minimap2 alignments so far. Consider reporting to the author if you have used a different \
read aligner so we can fix it."
        )
    result = []
    lines = aln.get_tag("SA").split(";")
    for line in lines:
        fields = line.split(",")
        if len(fields) < 6:
            continue
        # check if H not in cigar string, if so, return chrom,start
        result.append((fields[0], int(fields[1])))
    return result


def get_full_read_sequences_of_alignments(
    dict_alignments: dict[int, list[pysam.AlignedSegment]], path_alignments: Path
) -> dict[str, SeqRecord]:
    """Retrieves all DNA sequences of all given read alignments. The read DNA is in original orientation, as it was given in the original fasta file."""
    # iterate all alignments of a sampleID across all crIDs
    crIDs = dict_alignments.keys()
    dict_read_sequences: dict[
        str, SeqRecord
    ] = {}  # {(readname:(aln_reverse:bool,query_sequence:str))}
    dict_supplementary_positions: dict[str, list[tuple[str, int]]] = {}
    for crID in crIDs:
        for aln in dict_alignments[crID]:
            if (
                not aln.query_sequence
                or aln.infer_read_length() != len(aln.query_sequence)
                and aln.query_name not in dict_supplementary_positions
            ):
                dict_supplementary_positions[aln.query_name] = (
                    get_all_supplementary_positions_with_full_read_sequence(aln=aln)
                )
    # assort all supplementary_positions to a dict of the form {chromosome:[start,start+1]}
    dict_positions = {
        chrom: []
        for chrom in {
            chr
            for positions in dict_supplementary_positions.values()
            for chr, pos in positions
        }
    }
    for _readname, positions in dict_supplementary_positions.items():
        for chr_, pos in positions:
            dict_positions[chr_].append(pos)
    # then sort the list of start positions in each chromosome
    for chrom in dict_positions.keys():
        dict_positions[chrom] = sorted(dict_positions[chrom])
    # then merge any that are within 50kb of the next one
    # to do so, build an intervaltree, where 50_000 is subtracted from start and added to end
    for chrom in dict_positions.keys():
        it = IntervalTree()
        for start in dict_positions[chrom]:
            it[max(0, start - 50_000) : start + 50_000] = 1
        # then merge the intervaltree
        it.merge_overlaps(data_reducer=lambda a, b: a + b)
        # then extract the merged intervals and fill dict_positions with (start,end) tuples
        dict_positions[chrom] = [(interval.begin, interval.end) for interval in it]
    # the regions are then used to fetch the alignments from the bam file
    regions = [
        (chrom, start, end)
        for chrom in dict_positions.keys()
        for start, end in dict_positions[chrom]
    ]
    with pysam.AlignmentFile(path_alignments, "rb") as f:
        for region in regions:
            for aln in f.fetch(region[0], region[1], region[2]):
                if aln.query_name in dict_supplementary_positions.keys():
                    if not aln.query_sequence or not aln.infer_read_length() == len(
                        aln.query_sequence
                    ):
                        continue
                    seq = (
                        Seq(aln.query_sequence).reverse_complement()
                        if aln.is_reverse
                        else Seq(aln.query_sequence)
                    )
                    qualities = {"phred_quality": None}
                    if aln.query_qualities:
                        qualities = (
                            {"phred_quality": aln.query_qualities[::-1]}
                            if aln.is_reverse
                            else {"phred_quality": aln.query_qualities}
                        )
                    dict_read_sequences[aln.query_name] = SeqRecord(
                        seq=seq,
                        letter_annotations=qualities,
                        id=aln.query_name,
                        name=aln.query_name,
                    )
            # check if all keys of dict_read_sequences have a SeqRecord
            # if so, break this loop
            if all(
                readname in dict_read_sequences
                for readname in dict_supplementary_positions.keys()
            ):
                break
    # check if for each readname in dict_supplementary_positions, there is a SeqRecord in dict_read_sequences
    for readname in dict_supplementary_positions.keys():
        if readname not in dict_read_sequences:
            raise ValueError(
                f"readname {readname} not found in dict_read_sequences. Make sure that at least one alignment of the read brings the full read sequence with it."
            )
    # for all other alignments, just add the sequence and reverse flag
    for crID in crIDs:
        for aln in dict_alignments[crID]:
            if aln.query_name not in dict_read_sequences:
                seq = (
                    Seq(aln.query_sequence).reverse_complement()
                    if aln.is_reverse
                    else Seq(aln.query_sequence)
                )
                qualities = (
                    {"phred_quality": aln.query_qualities}
                    if aln.query_qualities
                    else None
                )
                dict_read_sequences[aln.query_name] = SeqRecord(
                    seq=seq,
                    letter_annotations=qualities,
                    id=aln.query_name,
                    name=aln.query_name,
                )
    return dict_read_sequences


# =============================================================================
# ALIGNMENT INTERVAL PROCESSING
# =============================================================================


def get_read_alignments_for_crs(
    crs: list[datatypes.CandidateRegion], alignments: Path
) -> tuple[
    dict[int, list[pysam.AlignedSegment]], dict[int, list[pysam.AlignedSegment]]
]:
    """Returns a dict of the form crID:{sampleID:[alignments]}"""
    dict_alignments: dict[int, list[pysam.AlignedSegment]] = {}
    dict_alignments_wt: dict[int, list[pysam.AlignedSegment]] = {}
    _readnames_in_signals = {signal.readname for cr in crs for signal in cr.sv_signals}
    # if no read alignments can be found for one sample,
    for cr in crs:
        # check if cr is of type datatypes.CandidateRegion
        assert isinstance(cr, datatypes.CandidateRegion), (
            f"cr {cr} is not of type datatypes.CandidateRegion"
        )
        # crate a dict of sampleID:aug_readname
        with pysam.AlignmentFile(alignments, "rb") as f:
            for aln in f.fetch(cr.chr, max(cr.referenceStart, 0), cr.referenceEnd):
                if aln.is_secondary:
                    continue
                if cr.crID not in dict_alignments:
                    dict_alignments[cr.crID] = []
                dict_alignments[cr.crID].append(aln)
                # if aln.query_name in readnames_in_signals:
                #     if cr.crID not in dict_alignments:
                #         dict_alignments[cr.crID] = []
                #     dict_alignments[cr.crID].append(aln)
                # else:  # add reads as wt if they don't appear in the signals
                #     if cr.crID not in dict_alignments_wt:
                #         dict_alignments_wt[cr.crID] = []
                #     dict_alignments_wt[cr.crID].append(aln)
    return dict_alignments, dict_alignments_wt


def get_read_alignment_intervals_in_region(
    region_start: int,
    regions_end: int,
    alignments: list[pysam.AlignedSegment],
    buffer_clipped_length: int,
) -> dict[str, list[tuple[int, int, str, int, int]]]:
    """Extract read alignment intervals within a specific region."""
    # check if all elements in alignments are of type pysam.AlignedSegment
    for aln in alignments:
        assert isinstance(aln, pysam.AlignedSegment), (
            f"aln {aln} is not a pysam.AlignedSegment"
        )
    dict_intervals = {}
    for aln in alignments:
        cr_start, cr_end = region_start, regions_end
        read_start, read_end = util.get_interval_on_read_in_region(
            a=aln,
            start=cr_start,
            end=cr_end,
            buffer_clipped_length=buffer_clipped_length,
        )
        ref_start, ref_end = util.get_interval_on_ref_in_region(
            a=aln, start=read_start, end=read_end
        )
        # check read_start and read_end if they are numerical and end-start is greater than 0
        if read_start is not None and read_end is not None and read_end > read_start:
            if aln.query_name not in dict_intervals:
                dict_intervals[aln.query_name] = []
            dict_intervals[aln.query_name].append((
                read_start,
                read_end,
                aln.reference_name,
                ref_start,
                ref_end,
            ))
            # dict_intervals[read_aug_name].append((min(read_start,read_end),max(read_start,read_end),min(ref_start,ref_end),max(ref_start,ref_end)))
    return dict_intervals


# crs are connected candidate regions!
def get_read_alignment_intervals_in_cr(
    crs: list[datatypes.CandidateRegion],
    buffer_clipped_length: int,
    dict_alignments: dict[int, list[pysam.AlignedSegment]],
) -> dict[str, list[tuple[int, int, str, int, int]]]:
    """
    Extract read alignment intervals in candidate regions.

    dict_alignments is a dict of the form {crID:[pysam.AlignedSegment]}.
    Returns a dict of the form {readname:[(read_start,read_end,ref_chr,ref_start,ref_end)]}.
    """
    # check form of dict_alignments. i.e. are all keys an integer?
    # are all values of dict_alignments a list of pysam.AlignedSegment?
    for crID in dict_alignments.keys():
        assert isinstance(crID, int), f"crID {crID} is not an integer"
        # check if value is a list
        assert isinstance(dict_alignments[crID], list), (
            f"value {dict_alignments[crID]} is not a list, but a {type(dict_alignments[crID])}"
        )
        for aln in dict_alignments[crID]:
            assert isinstance(aln, pysam.AlignedSegment), (
                f"aln {aln} is not a pysam.AlignedSegment"
            )
    dict_all_intervals: dict[str, list[tuple[int, int, str, int, int]]] = {}
    cr_extents = {cr.crID: (cr.referenceStart, cr.referenceEnd) for cr in crs}
    # get the maximum insertion size of all original alignments in the candidate regions
    max_insertion_size = max(
        [sv.size for cr in crs for sv in cr.sv_signals if sv.sv_type == 0], default=0
    )
    used_buffer_clipped_length = max(buffer_clipped_length, max_insertion_size)
    for crID in dict_alignments.keys():
        # error: update is not correct. it should be appending.
        intervals = get_read_alignment_intervals_in_region(
            alignments=dict_alignments[crID],
            buffer_clipped_length=used_buffer_clipped_length,
            region_start=cr_extents[crID][0],
            regions_end=cr_extents[crID][1],
        )
        # to each readname (key) append the according interval
        for readname in intervals.keys():
            if readname not in dict_all_intervals:
                dict_all_intervals[readname] = []
            dict_all_intervals[readname].extend(intervals[readname])
    return dict_all_intervals


def get_max_extents_of_read_alignments_on_cr(
    dict_all_intervals: dict[str, list[tuple[int, int, str, int, int]]],
) -> dict[str, tuple[int, int, str, int, str, int]]:
    """Returns a dict of the form {readname:(read_start,read_end,ref_start,ref_end)}"""
    dict_max_extents = {}
    for readname, intervals in dict_all_intervals.items():
        min_read_start = min(
            start for (start, end, ref_chr, ref_start, ref_end) in intervals
        )
        max_read_end = max(
            end for (start, end, ref_chr, ref_start, ref_end) in intervals
        )
        min_ref = min(
            (ref_start, ref_chr)
            for (start, end, ref_chr, ref_start, ref_end) in intervals
        )
        max_ref = max(
            (ref_end, ref_chr)
            for (start, end, ref_chr, ref_start, ref_end) in intervals
        )
        dict_max_extents[readname] = (
            min_read_start,
            max_read_end,
            min_ref[1],
            min_ref[0],
            max_ref[1],
            max_ref[0],
        )
    return dict_max_extents


def trim_reads(
    dict_alignments: dict[int, list[pysam.AlignedSegment]],
    intervals: dict[str, tuple[int, int, str, int, str, int]],
    read_records: dict[str, SeqRecord],
) -> dict[str, SeqRecord]:
    """
    Trim reads based on alignment intervals.

    Returns a dict of the form {aug_readname:SeqRecord}.
    The description holds the cut positions on both read and reference.
    """
    # dict_alignments is of the form: dict_alignments[cr.crID]=[pysam.AlignedSegment]
    # returns a list of read sequence records
    # check if dict_alignments is of the correct form
    for crID in dict_alignments.keys():
        assert isinstance(crID, int), f"crID {crID} is not an integer"
        for aln in dict_alignments[crID]:
            assert isinstance(aln, pysam.AlignedSegment), (
                f"aln {aln} is not a pysam.AlignedSegment"
            )
    # check if intervals is of the correct form
    for readname in intervals.keys():
        assert isinstance(readname, str), f"readname {readname} is not a string"
        assert isinstance(intervals[readname], tuple), (
            f"intervals[readname] {intervals[readname]} is not a tuple"
        )
        assert len(intervals[readname]) == 6, (
            f"intervals[readname] {intervals[readname]} is not of length 6"
        )
        # check if each tuple is of types: int,int,str,int,str,int
        for _key, (
            read_start,
            read_end,
            ref_start_chr,
            ref_start,
            ref_end_chr,
            ref_end,
        ) in intervals.items():
            assert isinstance(read_start, int), (
                f"read_start {read_start} is not an integer"
            )
            assert isinstance(read_end, int), f"read_end {read_end} is not an integer"
            assert isinstance(ref_start_chr, str), (
                f"ref_start_chr {ref_start_chr} is not a string"
            )
            assert isinstance(ref_start, int), (
                f"ref_start {ref_start} is not an integer"
            )
            assert isinstance(ref_end_chr, str), (
                f"ref_end_chr {ref_end_chr} is not a string"
            )
            assert isinstance(ref_end, int), f"ref_end {ref_end} is not an integer"

    cut_reads: dict[str, SeqRecord] = {}
    for crID in dict_alignments.keys():
        for aln in dict_alignments[crID]:
            if aln.query_name in intervals:
                # fixed this :)
                # if aln.cigartuples[-1][0] ==5 or aln.cigartuples[-1][0] == 5:
                #     raise(f"hard clipped read {read_aug_name} found in alignments. The complete DNA sequence of the read cannot be extracted.")
                if aln.query_name in cut_reads:
                    continue
                if aln.infer_read_length() != len(read_records[aln.query_name]):
                    raise ValueError(
                        f"read length of {aln.query_name} does not match the length of the read record. read length from alignment: {aln.infer_read_length()}, read length from read record: {len(read_records[aln.query_name])}"
                    )
                start, end, ref_start_chr, ref_start, ref_end_chr, ref_end = intervals[
                    aln.query_name
                ]
                record = read_records[aln.query_name][start:end]
                record.name = aln.query_name
                record.id = aln.query_name
                chosen_ref_start = (
                    (ref_start_chr, ref_start)
                    if ref_start < ref_end
                    else (ref_end_chr, ref_end)
                )
                chosen_ref_end = (
                    (ref_start_chr, ref_start)
                    if ref_start >= ref_end
                    else (ref_end_chr, ref_end)
                )
                # log.info(f"{read_aug_name} cut from {start} to {end} on ref: {ref_start} to {ref_end}")
                record.description = f"crID={crID},start={start},end={end},ref_start_chr={chosen_ref_start[0]},ref_start={chosen_ref_start[1]},ref_end_chr={chosen_ref_end[0]},ref_end={chosen_ref_end[1]}"
                cut_reads[aln.query_name] = record
    return cut_reads


# =============================================================================
# CONSENSUS GENERATION WITH RACON
# =============================================================================


def find_representative_read(
    similarity_matrix: np.ndarray,
    sim_read_names: list[str],
    cluster_read_names: list[str],
) -> str:
    """Find the centroid read in a cluster using the similarity matrix.

    The centroid is the read that maximises the mean similarity to all other
    reads in the cluster (equivalently, minimises the mean distance).

    Args:
        similarity_matrix: Symmetric (n, n) similarity matrix with values in [0, 1].
        sim_read_names: Ordered read names corresponding to matrix rows/columns.
        cluster_read_names: Read names belonging to the cluster of interest.

    Returns:
        The name of the representative (centroid) read.
    """
    if len(cluster_read_names) == 1:
        return cluster_read_names[0]

    name_to_idx = {name: i for i, name in enumerate(sim_read_names)}
    cluster_indices = np.array([name_to_idx[r] for r in cluster_read_names])

    # Sub-matrix for the cluster
    sub_sim = similarity_matrix[np.ix_(cluster_indices, cluster_indices)]
    # Mean similarity of each read to all other reads in the cluster
    mean_similarities = sub_sim.mean(axis=1)
    best_local_idx = int(np.argmax(mean_similarities))
    chosen_read = cluster_read_names[best_local_idx]
    log.debug(
        f"Chosen representative read: {chosen_read} with mean similarity {mean_similarities[best_local_idx]:.4f} in cluster of size {len(cluster_read_names)}"
    )
    return chosen_read


def filter_ava_alignments_for_racon(
    ava_bam: Path,
    representative_read: str,
    output_sam: Path,
) -> int:
    """Filter an all-vs-all BAM to produce a SAM suitable for racon polishing.

    From the all-vs-all BAM, keeps only alignments where the reference (target)
    is ``representative_read``.  For each query read, only the longest alignment
    (by number of query bases consumed) is retained, guaranteeing at most one
    alignment per query name and avoiding racon's duplicate-sequence error.

    Returns the number of alignments written.
    """
    # Pass 1: collect the longest alignment per query read that maps to the representative
    best: dict[str, pysam.AlignedSegment] = {}
    with pysam.AlignmentFile(str(ava_bam), mode="rb") as bam:
        ref_names = bam.references
        if representative_read not in ref_names:
            log.warning(
                f"Representative read '{representative_read}' not found as a reference "
                f"in AVA BAM '{ava_bam}'."
            )
            return 0
        for aln in bam.fetch(representative_read):
            if aln.is_unmapped or aln.query_name is None:
                continue
            # Skip the self-alignment: the representative read aligned to itself
            # would appear in both the SAM and the reference FASTA, causing racon to
            # complain about a duplicate sequence with unequal data.
            if aln.query_name == representative_read:
                continue
            query_bases = aln.query_alignment_length or 0
            prev = best.get(aln.query_name)
            if prev is None or query_bases > (prev.query_alignment_length or 0):
                best[aln.query_name] = aln

    # Pass 2: write filtered alignments to SAM
    with pysam.AlignmentFile(str(ava_bam), mode="rb") as bam:
        header = bam.header
    with pysam.AlignmentFile(str(output_sam), mode="w", header=header) as out_sam:
        for aln in best.values():
            out_sam.write(aln)

    return len(best)


def make_consensus_with_racon_subsampled(
    reference_fasta: Path,
    sam_alignments: Path,
    reads_fasta: Path,
    name: str,
    threads: int,
    consensus_fasta_path: Path,
    max_alns: int = 20,
    verbose: bool = False,
) -> str | None:
    tmp_sam = tempfile.NamedTemporaryFile(
        prefix="filtered.", suffix=".sam", delete=False
    )
    subsampled_readnames: list[str] = subsample_alignments(
        input_bamfile=sam_alignments,
        output_samfile=Path(tmp_sam.name),
        number=max_alns,
        verbose=verbose,
    )
    if len(subsampled_readnames) == 0:
        log.warning(
            f"no reads could be subsampled for consensus of {name}. Returning the reference as consensus."
        )
        return str(next(SeqIO.parse(reference_fasta, "fasta")).seq)
    contig = make_consensus_with_racon(
        reference_fasta=reference_fasta,
        sam_alignments=Path(tmp_sam.name),
        reads_fasta=reads_fasta,
        name=name,
        threads=threads,
        consensus_fasta_path=consensus_fasta_path,
    )
    return contig


def make_consensus_with_racon(
    reference_fasta: Path,
    sam_alignments: Path,
    reads_fasta: Path,
    name: str,
    threads: int,
    consensus_fasta_path: Path,
) -> str | None:
    """
    Create a consensus sequence from aligned reads using racon.

    Tries to create a consensus sequence from a set of reads that were aligned
    to a reference or representant. If the produced consensus is empty, it returns None.
    """
    # check if reference_fasta really contains a sequence
    if not any(SeqIO.parse(reference_fasta, "fasta")):
        raise ValueError(
            f"reference_fasta {reference_fasta} does not contain any sequences"
        )
    cmd_racon = f"racon -t {threads} {str(reads_fasta)} {str(sam_alignments)} {str(reference_fasta)}"
    try:
        with open(consensus_fasta_path, "w") as f:
            subprocess.check_call(shlex.split(cmd_racon), stdout=f)
        return str(next(SeqIO.parse(consensus_fasta_path, "fasta")).seq)
    except subprocess.CalledProcessError as e:
        log.warning(
            f"racon failed for {name} with command\n{cmd_racon}\nand error: {e}"
        )
        return None
    except Exception as e:
        log.warning(
            f"racon failed to produce a consensus for {name} with command\n{cmd_racon}\nand error: {e}"
        )
        return None


def align_reads_to_record(
    reference: Path,
    reads: Path,
    path_alignments: Path,
    threads: int,
    aln_args: str = " --secondary=no -Y --sam-hit-only -U25,75",
) -> None:
    """
    Align reads to a reference sequence.

    Takes a SeqRecord and a path to a fastq file of reads and aligns the reads
    to the reference. The alignments are written to a sam file at path_alignments.
    """
    # index fasta file with samtools index
    subprocess.run(shlex.split(f"samtools faidx {reference}"), check=True)
    log.info("aligning reads to representative")
    cmd_align = shlex.split(
        f"minimap2 -a -x map-ont -t {threads} {aln_args} {reference} {reads}"
    )
    with open(path_alignments, "w") as f:
        subprocess.check_call(cmd_align, stdout=f)


# =============================================================================
# CONSENSUS GENERATION WITH LAMASSEMBLE
# =============================================================================

# cmd_lamassembly = "lamassemble --name {consensus_dict['ID']} --all -P {threads} -f fa -s 0 {strlamassemble_mat)} {tmp_reads.name}" # shoul dbe written to the opened tmp consensus file


def make_consensus_with_lamassemble(
    lamassemble_mat: Path,
    reads_file: Path,
    output: Path,
    consensus_name: str,
    threads: int,
    timeout: int,
    verbose: bool = False,
) -> str | None:
    """Assemble reads with lamassemble and return the consensus sequence as a string."""
    # try to run lamassemble. if it fails or the output is empty, return None.
    # lamassemble writes its output to the command line, so the output should be caught from there.
    cmd_lamassemble = f"lamassemble --name {consensus_name} --all -P {threads} -f fa -s 0 {str(lamassemble_mat)} {str(reads_file)}"
    log.info(
        f"Running lamassemble with command:\n{cmd_lamassemble}\nwith timeout of {timeout} seconds"
    )
    result: str = ""
    try:
        with open(output, "w") as f:
            subprocess.check_call(
                shlex.split(cmd_lamassemble), stdout=f, timeout=timeout
            )
        # read the output file and return the sequence as a string
        try:
            consensus_sequence = next(SeqIO.parse(output, "fasta"))
            if len(consensus_sequence.seq) == 0:
                if verbose:
                    log.warning(
                        f"lamassemble produced an empty consensus for {consensus_name}"
                    )
                return None
            result = str(consensus_sequence.seq)
        except StopIteration:
            if verbose:
                log.warning(
                    f"lamassemble failed to produce a consensus for {consensus_name}"
                )
            return None
    except subprocess.TimeoutExpired:
        if verbose:
            log.warning(
                f"lamassemble timed out for {consensus_name} after {timeout} seconds"
            )
        return None
    except subprocess.CalledProcessError as e:
        if verbose:
            log.warning(f"lamassemble failed for {consensus_name} with error: {e}")
        return None
    return result


# =============================================================================
# MAIN CONSENSUS PROCESSING ALGORITHMS
# =============================================================================


def final_consensus(
    samplename: str,
    min_indel_size: int,
    min_bnd_size: int,
    reads_fasta: Path,
    consensus_fasta_path: Path,
    consensus_sequence: str,
    ID: str,
    crIDs: list[int],
    original_regions: list[tuple[str, int, int]],
    threads: int,
    verbose: bool,
    tmp_dir_path: Path | None = None,
) -> consensus_class.Consensus | None:
    # check if the consensus_fasta_path or reads_fasta are empty. if so, raise an exception
    if os.path.getsize(consensus_fasta_path) == 0:
        raise ValueError(
            f"consensus fasta {consensus_fasta_path} is empty. Cannot proceed."
        )
    if os.path.getsize(reads_fasta) == 0:
        raise ValueError(f"reads fasta {reads_fasta} is empty. Cannot proceed.")

    min_cr = min(crIDs)
    tmp_alignments = tempfile.NamedTemporaryFile(
        dir=tmp_dir_path,
        prefix=f"consensus.{str(min_cr)}.final_tmp_alignments.",
        suffix=".bam",
        delete=False if tmp_dir_path else True,
    )

    util.align_reads_with_minimap(
        bamout=tmp_alignments.name,
        reads=reads_fasta,
        reference=consensus_fasta_path,
        aln_args=" -Y --sam-hit-only --secondary=no  -U 10,25 -H",
        threads=threads,
    )

    # Parse alignments using context manager
    with pysam.AlignmentFile(tmp_alignments.name, mode="r") as aln_file:
        cut_read_alns: list[pysam.AlignedFragment] = list(aln_file)
    if len(cut_read_alns) == 0:
        log.warning(
            f"No alignments found for consensus {ID} in {tmp_alignments.name}. Returning None."
        )
        return None

    # Extract SV signals from alignments
    cut_read_alignment_signals = [
        parse_ReadAlignmentSignals_from_alignment(
            alignment=aln,
            min_signal_size=min_indel_size,
            min_bnd_size=min_bnd_size,
            samplename=samplename,
        )
        for aln in cut_read_alns
    ]

    # Extract alignment intervals
    intervals_cutread_alignments = [
        (aln.reference_start, aln.reference_end, aln.query_name, aln.is_forward)
        for aln in cut_read_alns
    ]
    # if they are empty, then raise an exception
    if len(intervals_cutread_alignments) == 0:
        raise ValueError(
            f"No alignment intervals found for consensus {ID} in {tmp_alignments.name}. Cannot proceed."
        )
    if len(cut_read_alignment_signals) == 0:
        log.warning(
            f"No alignment signals found for consensus {ID} -> no distortions available?. Returning None."
        )

    # Create final Consensus object
    result = consensus_class.Consensus(
        ID=ID,
        crIDs=crIDs,
        original_regions=original_regions,
        consensus_sequence=consensus_sequence,
        cut_read_alignment_signals=cut_read_alignment_signals,
        intervals_cutread_alignments=intervals_cutread_alignments,
    )
    log.debug(
        f"Consensus object created. ID:{result.ID}, \
            crIDs:{result.crIDs}, \
            original_regions:{result.original_regions}, \
            consensus_sequence length:{len(result.consensus_sequence) if result.consensus_sequence else 'None'}, \
            cut read intervals:{len(result.intervals_cutread_alignments)}"
    )

    if verbose:
        alignments_to_rafs.display_ascii_alignments(
            alignments=cut_read_alns, terminal_width=125
        )

    return result


def add_cutread_alignments_to_consensus_inplace(
    alns: list[pysam.AlignedSegment],
    consensus: consensus_class.Consensus,
    min_indel_size: int,
    min_bnd_size: int,
    samplename: str,
) -> None:
    cut_read_alignment_signals = [
        parse_ReadAlignmentSignals_from_alignment(
            alignment=aln,
            min_signal_size=min_indel_size,
            min_bnd_size=min_bnd_size,
            samplename=samplename,
        )
        for aln in alns
    ]
    intervals_cutread_alignments = [
        (aln.reference_start, aln.reference_end, aln.query_name, aln.is_forward)
        for aln in alns
    ]
    consensus.cut_read_alignment_signals.extend(cut_read_alignment_signals)
    consensus.intervals_cutread_alignments.extend(intervals_cutread_alignments)
    log.debug(
        f"Added {len(alns)} alignments to consensus {consensus.ID}. Now has {len(consensus.cut_read_alignment_signals)} signals and {len(consensus.intervals_cutread_alignments)} intervals."
    )


def deterministic_hash(s: str) -> int:
    """Create a deterministic hash from a string."""
    return int(hashlib.md5(s.encode()).hexdigest()[:8], 16)


# def estimate_adaptive_separation_threshold(
#         alignment_scores: dict[str, float],
#         expected_clusters:int,
#         default_threshold: float = 5.0,
#         outlier_fraction=0.1,
#         min_scores_required: int = 4,
#         min_result=3.0,
#         verbose: bool = False) -> float:
#     expected_clusters = 2 #max(1, expected_clusters) # TODO: number of haplotypes should be tied to the copy number at this point.
#     if len(alignment_scores) <= min_scores_required - 1:
#         # Not enough scores to calculate meaningful pairwise differences
#         return default_threshold
#     if outlier_fraction < 0.0 or outlier_fraction > 1.0:
#         raise ValueError("outlier_fraction must be between 0 and 1")
#     if not min_result > 0.0:
#         raise ValueError("min_result must be greater than 0.0")
#     # Sort scores and calculate pairwise differences
#     sorted_scores:list[float] = sorted(alignment_scores.values())
#     pairwise_diffs:list[float] = [abs(sorted_scores[i+1] - sorted_scores[i])
#                       for i in range(len(sorted_scores) - 1)]

#     if len(pairwise_diffs) < min_scores_required - 1:
#         return default_threshold

#     # remove the expected outliers
#     pairwise_diffs: list[float] = sorted(pairwise_diffs)[int(round((len(pairwise_diffs) * (1.0 -outlier_fraction)))):]

#     # expected_clusters cannot be larger than the number of pairwise differences left after outlier removal
#     # and must be at least 1
#     expected_clusters = min(max(1, expected_clusters), len(pairwise_diffs))

#     # Ensure we don't try to access beyond the list bounds
#     if len(pairwise_diffs) == 0:
#         return default_threshold

#     # Calculate the percentile of pairwise differences
#     result = pairwise_diffs[-expected_clusters]

#     if verbose:
#         log.info(f"Adaptive separation threshold base value: {result:.3f} with given expected_clusters: {expected_clusters}, outlier_fraction: {outlier_fraction}, min_scores_required: {min_scores_required}, total scores: {len(alignment_scores)}")
#     result = max(result, min_result)
#     return result


# =============================================================================
# ITERATIVE CONSENSUS CLUSTERING WITH RACON OR LAMASSEMBLE
# =============================================================================


@attrs.define
class ConsensusClusteringObject:
    reference_fasta: Path
    reads_fasta: Path
    sam_alignments: Path
    indexed_bam: Path | None
    consensus_fasta: Path | None
    consensus_name: str | None
    chosen_readnames: set[str] = set()

    def __repr__(self) -> str:
        return f"ConsensusClusteringObject(reference_fasta={self.reference_fasta}, reads_fasta={self.reads_fasta}, sam_alignments={self.sam_alignments}, indexed_bam={self.indexed_bam}, consensus_fasta={self.consensus_fasta}, consensus_name={self.consensus_name}, chosen_readnames={self.chosen_readnames})"


def filter_consensus_paths_object_by_readnames(
    cpo: ConsensusClusteringObject, chosen_reads: list[str]
) -> None:
    """Filter the reads and alignments in a ConsensusClusteringObject by a list of chosen read names."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        reads_fasta = cpo.reads_fasta  # they are really fastq files
        sam_alignments = cpo.sam_alignments
        # at first, copy the reads_fasta and sam_alignments to the tmp_dir
        tmp_reads_fasta = tempfile.NamedTemporaryFile(
            prefix="reads.", suffix=".fasta", dir=tmp_dir, delete=False
        )
        shutil.copyfile(reads_fasta, tmp_reads_fasta.name)
        tmp_sam_alignments = tempfile.NamedTemporaryFile(
            prefix="alignments.", suffix=".sam", dir=tmp_dir, delete=False
        )
        shutil.copyfile(sam_alignments, tmp_sam_alignments.name)
        # now, filter the reads and alignments
        filter_reads_fastq_by_readnames(
            input_fastq=Path(tmp_reads_fasta.name),
            output_fastq=reads_fasta,
            readnames=set(chosen_reads),
        )
        filter_sam_by_readnames(
            input_samfile=Path(tmp_sam_alignments.name),
            output_samfile=sam_alignments,
            readnames=set(chosen_reads),
        )


def align_cutreads_to_representative(
    representative: str,
    readnames: list[str],
    cutreads: dict[str, datatypes.SeqRecord],
    tmp_dir: str | Path,
    compress: bool = True,
    show_ascii_alignments: bool = False,
    threads: int = 1,
) -> ConsensusClusteringObject | None:
    """Write the representative read and cut reads to separate files.
    Args:
        representative (str): DNA sequence of the representative read
        readnames (list[str]): list of read names to write to the reads fasta
        cutreads (dict[str, datatypes.SeqRecord]): dict of cut reads (name:SeqRecord)
        tmp_dir (str | Path): path to the temporary directory

    Returns:
        tuple[Path, Path, Path, Path]: paths to the (indexed) representative fasta,
        the reads fasta, and the alignments paths (sam and indexed bam).
    """
    # write representative to tmp fasta
    tmp_ref = tempfile.NamedTemporaryFile(
        prefix="representative.", suffix=".fasta", dir=tmp_dir, delete=False
    )
    with open(tmp_ref.name, "w") as f:
        SeqIO.write(cutreads[representative], f, "fasta")
    cmd_index = f"samtools faidx {tmp_ref.name}"
    subprocess.check_call(shlex.split(cmd_index))
    # write reads to tmp fasta
    tmp_reads = tempfile.NamedTemporaryFile(
        prefix="reads.", suffix=".fasta", dir=tmp_dir, delete=False
    )
    with open(tmp_reads.name, "w") as f:
        SeqIO.write([cutreads[readname] for readname in readnames], f, "fasta")
    # align reads to reference
    tmp_alignments_sam = tempfile.NamedTemporaryFile(
        prefix="alignments.", suffix=".sam", dir=tmp_dir, delete=False
    )
    tmp_alignments_bam = tempfile.NamedTemporaryFile(
        prefix="alignments.", suffix=".bam", dir=tmp_dir, delete=False
    )
    align_reads_to_record(
        reference=tmp_ref.name,
        reads=tmp_reads.name,
        path_alignments=tmp_alignments_sam.name,
        threads=threads,
    )

    if show_ascii_alignments:
        with pysam.AlignmentFile(tmp_alignments_sam.name, mode="r") as aln_file:
            alignments_to_rafs.display_ascii_alignments(alignments=list(aln_file))

    if compress:
        # convert sam to bam and index
        cmd_view = f"samtools view -b {tmp_alignments_sam.name}"
        with open(tmp_alignments_bam.name, "wb") as f:
            subprocess.check_call(shlex.split(cmd_view), stdout=f)
        cmd_index_bam = f"samtools index {tmp_alignments_bam.name}"
        subprocess.check_call(shlex.split(cmd_index_bam))
    return ConsensusClusteringObject(
        reference_fasta=Path(tmp_ref.name),
        reads_fasta=Path(tmp_reads.name),
        sam_alignments=Path(tmp_alignments_sam.name),
        indexed_bam=Path(tmp_alignments_bam.name) if compress else None,
    )


def filter_reads_fastq_by_readnames(
    input_fastq: Path, output_fastq: Path, readnames: set[str]
) -> None:
    """Filter a fasta file by read names."""
    with open(output_fastq, "w") as out_f:
        for record in SeqIO.parse(input_fastq, "fasta"):
            if record.id in readnames:
                SeqIO.write(record, out_f, "fasta")


def prepare_and_run_racon(
    reads_fasta: Path,
    consensus_fasta_path: Path,
    name: str,
    threads: int,
    ava_bam: Path,
    representative_read: str,
    tmp_dir_path: Path | None = None,
    verbose: bool = False,
) -> str | None:
    """Extract alignments to the representative from the AVA BAM and polish with racon.

    Args:
        reads_fasta: Path to a FASTQ file containing the cluster reads.
        consensus_fasta_path: Path where the consensus FASTA will be written.
        name: Name for the consensus sequence.
        threads: Number of threads for racon.
        ava_bam: Path to the all-vs-all alignments BAM file.
        representative_read: Name of the chosen representative (centroid) read.
        tmp_dir_path: Optional directory for temporary files (kept on disk if set).
        verbose: Enable verbose logging.

    Returns:
        The consensus sequence string, or None on failure.
    """
    log.info(f"prepare and run consensus assembly with racon for '{name}'")
    records = {r.id: r for r in SeqIO.parse(reads_fasta, "fasta")}
    if len(records) == 0:
        log.warning(f"No reads in {reads_fasta} for racon consensus '{name}'.")
        return None
    if len(records) == 1:
        rec = next(iter(records.values()))
        with open(consensus_fasta_path, "w") as f:
            SeqIO.write(rec, f, "fasta")
        return str(rec.seq)

    if representative_read not in records:
        log.warning(
            f"Representative read '{representative_read}' not found in reads_fasta for '{name}'. Falling back to longest read."
        )
        representative_read = max(records, key=lambda r: len(records[r].seq))

    if verbose:
        log.info(
            f"racon: using '{representative_read}' "
            f"(length {len(records[representative_read].seq)}) as reference for '{name}'"
        )

    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tmp_dir:
        # 1. Write the representative read to a temporary FASTA file.
        tmp_ref = tempfile.NamedTemporaryFile(
            prefix="racon_ref.",
            suffix=".fasta",
            dir=tmp_dir,
            delete=False if tmp_dir_path else True,
        )
        with open(tmp_ref.name, "w") as f:
            SeqIO.write(records[representative_read], f, "fasta")

        # 2. Extract AVA alignments to the representative, one per query (longest kept).
        tmp_sam = tempfile.NamedTemporaryFile(
            prefix="racon_aln.",
            suffix=".sam",
            dir=tmp_dir,
            delete=False if tmp_dir_path else True,
        )
        n_alns = filter_ava_alignments_for_racon(
            ava_bam=ava_bam,
            representative_read=representative_read,
            output_sam=Path(tmp_sam.name),
        )
        if n_alns == 0:
            log.warning(
                f"No alignments found to representative '{representative_read}' "
                f"in AVA BAM '{ava_bam}' for '{name}'."
            )
            return None

        # 3. Polish with racon (reads_fasta is already FASTA — no quality strings).
        consensus_sequence = make_consensus_with_racon(
            reference_fasta=Path(tmp_ref.name),
            sam_alignments=Path(tmp_sam.name),
            reads_fasta=reads_fasta,
            name=name,
            threads=threads,
            consensus_fasta_path=consensus_fasta_path,
        )
    return consensus_sequence


# This function is run after cut reads of one cluster have been aligned to their representative. This means that there are
# cut reads of one cluster written to a file (use for lamassemble or racon), and alignments in sam format to be used with
# racon, and the representative read in fasta format - to be used with racon.
def assemble_consensus(
    lamassemble_mat: Path | None | str,
    name: str,
    reads_fasta: Path,
    consensus_fasta_path: Path,
    threads: int,
    timeout: int,
    method: str,
    ava_bam: Path | None = None,
    representative_read: str | None = None,
    tmp_dir_path: Path | None = None,
    verbose: bool = False,
) -> str | None:
    if method not in ("lamassemble", "racon"):
        raise ValueError(
            f"Unknown consensus method '{method}'. Choose 'lamassemble' or 'racon'."
        )

    with tempfile.TemporaryDirectory():
        if method == "racon":
            if ava_bam is None or representative_read is None:
                log.warning(
                    "ava_bam and representative_read are required when method='racon'. "
                    "Returning None."
                )
                return None
            consensus_sequence: str | None = prepare_and_run_racon(
                reads_fasta=reads_fasta,
                consensus_fasta_path=consensus_fasta_path,
                name=name,
                threads=threads,
                ava_bam=ava_bam,
                representative_read=representative_read,
                tmp_dir_path=tmp_dir_path,
                verbose=verbose,
            )
            if consensus_sequence is not None:
                return consensus_sequence
        else:
            consensus_sequence = make_consensus_with_lamassemble(
                lamassemble_mat=lamassemble_mat,
                reads_file=reads_fasta,
                output=consensus_fasta_path,
                consensus_name=name,
                timeout=timeout,
                threads=threads,
                verbose=verbose,
            )
            if consensus_sequence is not None:
                return consensus_sequence
    return None


def partition_reads_spectral(
    similarity_matrix: np.ndarray, read_names: list[str], n_clusters: int
) -> dict[str, int]:
    """Cluster reads via spectral clustering on a precomputed similarity matrix.

    Args:
        similarity_matrix: Symmetric (n, n) similarity matrix with values in [0, 1].
        read_names: Ordered list of read names corresponding to matrix rows/columns.
        n_clusters: Number of clusters to produce.

    Returns:
        Dict mapping each read name to its cluster label (0-indexed).
    """
    clustering = SpectralClustering(
        n_clusters=n_clusters,
        affinity="precomputed",
        assign_labels="kmeans",
        random_state=42,
    )

    labels = clustering.fit_predict(similarity_matrix)
    return {read: int(label) for read, label in zip(read_names, labels, strict=True)}


def consensus_while_clustering(
    samplename: str,
    lamassemble_mat: Path | str | None,
    pool: dict[str, SeqRecord],
    candidate_regions: dict[int, datatypes.CandidateRegion],
    partitions: int,
    timeout: int,
    consensus_method: str,
    threads: int = 1,
    tmp_dir_path: Path | None = None,
    figures_dir: Path | None = None,
    verbose: bool = False,
    densities_weight: float = 1.0,
    max_intra_distance: float = -1.0,
) -> dict[str, consensus_class.Consensus] | None:
    # all-vs-all alignments with minimap2
    # calculate penalties (called score_ras_from_alignments, but its really not a score, but a penalty)
    # construct a graph of reads with edges between reads given by the alignment penalties. If no penalty is found, there is no edge.
    # find a partition of the graph given N clusters that minimizes the sum of edge weights between clusters and contains all reads (min coverage set)
    crIDs = [cr.crID for cr in candidate_regions.values()]
    result: dict[str, consensus_class.Consensus] | None = None
    max_ava_attempts = 4
    ava_pool = dict(pool)  # mutable copy for subsampling
    excluded_reads: list[str] = []  # reads dropped by subsampling
    try:
        with tempfile.TemporaryDirectory(
            dir=tmp_dir_path, delete=False if tmp_dir_path else True
        ) as tmp_dir:
            # Retry all-vs-all alignment with subsampling on failure
            ava_succeeded = False
            for ava_attempt in range(max_ava_attempts):
                if len(ava_pool) < 2:
                    log.warning(
                        f"Only {len(ava_pool)} read(s) remaining after subsampling. "
                        "Cannot perform all-vs-all alignment."
                    )
                    break

                if ava_attempt > 0:
                    # Subsample 50% of the current pool
                    read_names_list = sorted(ava_pool.keys())
                    n_keep = max(2, len(read_names_list) // 2)
                    rng = random.Random(42 + ava_attempt)
                    kept = rng.sample(read_names_list, n_keep)
                    dropped = [r for r in read_names_list if r not in set(kept)]
                    excluded_reads.extend(dropped)
                    ava_pool = {r: ava_pool[r] for r in kept}
                    log.info(
                        f"AVA attempt {ava_attempt + 1}/{max_ava_attempts}: "
                        f"subsampled to {len(ava_pool)} reads "
                        f"({len(excluded_reads)} reads excluded so far)."
                    )

                # 1) write all read sequences to a temporary fasta file
                tmp_all_reads = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="all_reads.",
                    suffix=".fasta",
                    delete=False if tmp_dir_path else True,
                )
                with open(tmp_all_reads.name, "w") as f:
                    SeqIO.write(ava_pool.values(), f, "fasta")
                # 2) align all reads to each other with minimap2
                tmp_all_vs_all_sam = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="all_vs_all.",
                    suffix=".bam",
                    delete=False if tmp_dir_path else True,
                )
                try:
                    util.align_reads_with_minimap(
                        timeout=timeout,
                        bamout=tmp_all_vs_all_sam.name,
                        reads=tmp_all_reads.name,
                        reference=tmp_all_reads.name,
                        aln_args=" --sam-hit-only --secondary=yes -U 25,35 -H",
                        tech="ava-ont",
                        threads=threads,
                    )
                except (
                    TimeoutError,
                    subprocess.TimeoutExpired,
                    subprocess.CalledProcessError,
                ) as e:
                    log.warning(
                        f"AVA alignment failed on attempt {ava_attempt + 1}/{max_ava_attempts} "
                        f"with {len(ava_pool)} reads: {e}"
                    )
                    continue

                # 3) parse alignments using context manager
                with pysam.AlignmentFile(tmp_all_vs_all_sam.name, mode="r") as aln_file:
                    all_vs_all_alignments = list(aln_file)
                if len(all_vs_all_alignments) == 0:
                    log.warning(
                        f"No alignments found in AVA on attempt {ava_attempt + 1}/{max_ava_attempts}. "
                        "Retrying with fewer reads."
                    )
                    continue

                ava_succeeded = True
                break

            if not ava_succeeded:
                log.error(
                    f"All {max_ava_attempts} AVA alignment attempts failed. Returning None."
                )
                return None

            if len(excluded_reads) > 0:
                log.info(
                    f"AVA alignment succeeded with {len(ava_pool)} of {len(pool)} reads. "
                    f"{len(excluded_reads)} excluded reads will be added to unused reads."
                )

            # 4) Build similarity matrix via consensus_lib pipeline
            ava_path = Path(tmp_all_vs_all_sam.name)
            all_reads = list(ava_pool.keys())
            read_lengths = {name: len(rec.seq) for name, rec in ava_pool.items()}
            if verbose:
                # print all read lengths
                for name, length in read_lengths.items():
                    log.info(f"Read {name} has length {length}")

            if figures_dir is not None:
                figures_dir.mkdir(parents=True, exist_ok=True)
                consensus_lib.visualize_ava_alignments(
                    ava_path, figures_dir / "ava_alignment_lengths.png"
                )
                consensus_lib.visualize_alignment_presence(
                    ava_path, figures_dir / "alignment_presence.png"
                )

            ava_signals = consensus_lib.parse_sv_signals_from_ava_alignments(
                ava_alignments=ava_path,
                min_signal_size=12,
                min_bnd_size=100,
            )
            size_densities = consensus_lib.sv_size_densities_from_ava_signals(
                ava_signals, figures_dir=figures_dir
            )
            densities = consensus_lib.importance_densities_from_ava_signals(
                ava_signals=ava_signals,
                size_densities=size_densities,
            )
            directed_signals = consensus_lib.parse_directed_signals_from_ava(
                ava_alignments=ava_path,
                min_signal_size=12,
                min_bnd_size=100,
            )

            similarity_matrix, sim_read_names = (
                consensus_lib.pairwise_similarity_matrix(
                    directed_signals=directed_signals,
                    densities=densities,
                    read_lengths=read_lengths,
                    all_read_names=all_reads,
                    densities_weight=densities_weight,
                )
            )

            if figures_dir is not None:
                consensus_lib.visualize_importance_densities(
                    densities=densities,
                    output=figures_dir / "importance_densities.png",
                )
                consensus_lib.visualize_size_similarity(
                    read_lengths=read_lengths,
                    output=figures_dir / "size_similarity.png",
                )
                consensus_lib.visualize_raw_signal_matrix(
                    directed_signals=directed_signals,
                    all_read_names=sim_read_names,
                    output=figures_dir / "raw_signal_matrix.png",
                )
                consensus_lib.visualize_normalized_similarity_matrix(
                    similarity_matrix=similarity_matrix,
                    read_names=sim_read_names,
                    output=figures_dir / "normalized_similarity_matrix.png",
                )

            # Detect outlier reads
            outliers = consensus_lib.detect_outlier_reads(
                similarity=similarity_matrix,
                read_names=sim_read_names,
                read_lengths=read_lengths,
                n_clusters=partitions,
            )
            isolated = list(outliers)
            # Add reads that were excluded during AVA subsampling
            isolated.extend(excluded_reads)
            outlier_set = set(outliers)
            well_connected = [r for r in sim_read_names if r not in outlier_set]
            log.debug(
                f"Well-connected reads: {len(well_connected)}, Outlier reads: {len(isolated)}"
            )

            if len(well_connected) < partitions:
                effective_partitions = max(1, len(well_connected))
            else:
                effective_partitions = partitions
            log.debug(
                f"Effective number of clusters for spectral clustering: {effective_partitions}. Original requested partitions: {partitions}"
            )

            # Cluster well-connected reads using their sub-matrix
            if len(well_connected) > 1:
                wc_indices = [
                    i
                    for i, name in enumerate(sim_read_names)
                    if name not in outlier_set
                ]
                wc_similarity = similarity_matrix[np.ix_(wc_indices, wc_indices)]
                clustering_result = partition_reads_spectral(
                    similarity_matrix=wc_similarity,
                    read_names=well_connected,
                    n_clusters=effective_partitions,
                )
            else:
                clustering_result = {well_connected[0]: 0} if well_connected else {}

            if figures_dir is not None:
                _cluster_palette = [
                    "#e41a1c",
                    "#377eb8",
                    "#4daf4a",
                    "#984ea3",
                    "#ff7f00",
                    "#a65628",
                    "#f781bf",
                ]
                node_colors: dict[str, str] = {
                    name: _cluster_palette[cid % len(_cluster_palette)]
                    for name, cid in clustering_result.items()
                }
                for name in isolated:
                    node_colors[name] = "#aaaaaa"
                consensus_lib.visualize_similarity_graph(
                    similarity_matrix=similarity_matrix,
                    read_names=sim_read_names,
                    output=figures_dir / "similarity_graph.png",
                    node_colors=node_colors,
                )

            # 3. Proceed with consensus generation for each cluster
            for cluster_id in set(clustering_result.values()):
                chosen_reads = [
                    read for read, cid in clustering_result.items() if cid == cluster_id
                ]
                consensus_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="consensus.",
                    suffix=".fasta",
                    delete=False if tmp_dir else True,
                )
                reads_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="reads.",
                    suffix=".fasta",
                    delete=False if tmp_dir else True,
                )
                with open(reads_fasta.name, "w") as f:
                    SeqIO.write(
                        [pool[readname] for readname in chosen_reads], f, "fasta"
                    )
                consensus_name = f"{min(crIDs)}.{cluster_id}"

                # Find centroid read for racon
                _representative = (
                    find_representative_read(
                        similarity_matrix=similarity_matrix,
                        sim_read_names=sim_read_names,
                        cluster_read_names=chosen_reads,
                    )
                    if consensus_method == "racon"
                    else None
                )

                consensus_sequence: str | None = assemble_consensus(
                    lamassemble_mat=lamassemble_mat,
                    name=consensus_name,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    threads=threads,
                    timeout=timeout,
                    method=consensus_method,
                    ava_bam=ava_path if consensus_method == "racon" else None,
                    representative_read=_representative,
                    tmp_dir_path=tmp_dir_path,
                    verbose=verbose,
                )

                if consensus_sequence is None:
                    log.warning(
                        f"consensus assembly failed for cluster {cluster_id} with {len(chosen_reads)} reads. Skipping this cluster."
                    )
                    # add the chosen reads to the isolated reads
                    isolated.extend(chosen_reads)
                    continue

                if result is None:
                    result = {}
                consensus = final_consensus(
                    samplename=samplename,
                    min_indel_size=8,
                    min_bnd_size=100,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    consensus_sequence=consensus_sequence,
                    ID=consensus_name,
                    crIDs=crIDs,
                    original_regions=[
                        (cr.chr, cr.referenceStart, cr.referenceEnd)
                        for cr in candidate_regions.values()
                    ],
                    threads=threads,
                    verbose=verbose,
                    tmp_dir_path=Path(tmp_dir),
                )
                if consensus is not None:
                    result[consensus_name] = consensus
                else:
                    log.warning(
                        f"final consensus generation failed for cluster {cluster_id} with {len(chosen_reads)} reads. Skipping this cluster."
                    )
                    # add the chosen reads to the isolated reads
                    isolated.extend(chosen_reads)
                    continue
            # if there is no cluster, but only isolated reads, try to rescue them by trying to create a consensus from them
            # write all isolated reads to one fasta file
            # create a tmp consensus fasta file
            if len(isolated) > 0 and result is None or len(result) == 0:
                log.info(
                    f"Trying to rescue {len(isolated)} isolated reads by creating one consensus from them."
                )
                consensus_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="consensus.isolated.",
                    suffix=".fasta",
                    delete=False if tmp_dir else True,
                )
                reads_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="reads.isolated.",
                    suffix=".fasta",
                    delete=False if tmp_dir else True,
                )
                with open(reads_fasta.name, "w") as f:
                    SeqIO.write(
                        [pool[readname] for readname in isolated if readname in pool],
                        f,
                        "fasta",
                    )
                consensus_name = f"{min(crIDs)}.isolated"

                # Find centroid read for racon (isolated reads)
                _isolated_in_pool = [r for r in isolated if r in pool]
                _representative = (
                    find_representative_read(
                        similarity_matrix=similarity_matrix,
                        sim_read_names=sim_read_names,
                        cluster_read_names=_isolated_in_pool,
                    )
                    if consensus_method == "racon" and len(_isolated_in_pool) > 0
                    else None
                )

                consensus_sequence: str | None = assemble_consensus(
                    lamassemble_mat=lamassemble_mat,
                    name=consensus_name,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    threads=threads,
                    timeout=timeout,
                    method=consensus_method,
                    ava_bam=ava_path if consensus_method == "racon" else None,
                    representative_read=_representative,
                    tmp_dir_path=tmp_dir_path,
                    verbose=verbose,
                )

                if consensus_sequence is not None:
                    if result is None:
                        result = {}
                    consensus = final_consensus(
                        samplename=samplename,
                        min_indel_size=8,
                        min_bnd_size=100,
                        reads_fasta=Path(reads_fasta.name),
                        consensus_fasta_path=Path(consensus_fasta.name),
                        consensus_sequence=consensus_sequence,
                        ID=consensus_name,
                        crIDs=crIDs,
                        original_regions=[
                            (cr.chr, cr.referenceStart, cr.referenceEnd)
                            for cr in candidate_regions.values()
                        ],
                        threads=threads,
                        verbose=verbose,
                        tmp_dir_path=Path(tmp_dir),
                    )
                    if consensus is not None:
                        result[consensus_name] = consensus
                    else:
                        log.warning(
                            f"final consensus generation failed for isolated reads with {len(isolated)} reads. Returning None."
                        )
                else:
                    log.warning(
                        f"consensus assembly failed for isolated reads with {len(isolated)} reads. Returning None."
                    )

    except Exception as e:
        log.warning(
            f"consensus while clustering failed with exception: {e}. Returning None."
        )
    return result


def consensus_while_clustering_with_kmeans(
    samplename: str,
    dict_summed_indels: dict[str, list[int]],
    lamassemble_mat: Path | str | None,
    pool: dict[str, SeqRecord],
    candidate_regions: dict[int, datatypes.CandidateRegion],
    max_k: int,
    variance_threshold: float,
    distance_threshold: float,
    consensus_method: str,
    threads: int = 1,
    tmp_dir_path: Path | None = None,
    timeout: int = 120,
    verbose: bool = False,
) -> dict[str, consensus_class.Consensus] | None:
    """Cluster reads by their summed indel distribution using KMeans and assemble consensus per cluster.

    If there is a very clear separation in dict_summed_indels (2 dimensions: sum insertions,
    sum deletions), the all-vs-all alignment step can be skipped and reads can be directly
    assembled per cluster.

    The variance within a cluster should be low (below variance_threshold), and the distance
    between the cluster centroids should be high (above distance_threshold).

    Returns None if no good clustering is found (caller should fall back to
    consensus_while_clustering).
    """
    crIDs = [cr.crID for cr in candidate_regions.values()]

    # Need at least 2 reads for meaningful clustering
    if len(dict_summed_indels) < 2:
        log.debug("Not enough reads for KMeans clustering. Returning None.")
        return None

    # Filter out reads whose sequence length is an outlier.
    # Uses the same gap-based approach as consensus_lib._outlier_reads_by_length:
    # reads are sorted by length, gaps exceeding `length_factor` start a new
    # cluster, and clusters smaller than a threshold are flagged as outliers.
    read_lengths = {rn: len(pool[rn].seq) for rn in dict_summed_indels if rn in pool}
    size_outliers: set[str] = set()
    if len(read_lengths) >= 3:
        length_factor = 1.2
        sorted_by_len = sorted(read_lengths, key=lambda r: read_lengths[r])
        clusters: list[list[str]] = [[sorted_by_len[0]]]
        for i in range(1, len(sorted_by_len)):
            prev_len = read_lengths[sorted_by_len[i - 1]]
            curr_len = read_lengths[sorted_by_len[i]]
            if prev_len > 0 and curr_len / prev_len > length_factor:
                clusters.append([])
            clusters[-1].append(sorted_by_len[i])
        min_cluster_size = max(2, int(np.sqrt(len(read_lengths) / max(max_k, 1))))
        for cluster in clusters:
            if len(cluster) < min_cluster_size:
                size_outliers.update(cluster)
        if size_outliers:
            log.info(
                f"KMeans: filtered {len(size_outliers)} size-outlier read(s): "
                f"{sorted(size_outliers)}"
            )

    # Build the feature matrix: each read has [sum_insertions, sum_deletions]
    readnames = sorted(
        rn for rn in dict_summed_indels.keys() if rn not in size_outliers
    )

    if len(readnames) < 2:
        log.debug(
            "Not enough reads after size-outlier filtering for KMeans clustering. Returning None."
        )
        return None

    X = np.array([dict_summed_indels[rn] for rn in readnames], dtype=np.float64)

    # 1) Try KMeans clustering with k = 1 .. max_k
    #    Greedily pick the first k where intra-cluster variance is low
    #    and inter-cluster distance is high.
    chosen_k: int | None = None
    chosen_labels: np.ndarray | None = None

    for k in range(1, max_k + 1):
        if k > len(readnames):
            break

        kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
        labels = kmeans.fit_predict(X)
        centroids = kmeans.cluster_centers_

        # Compute max intra-cluster variance (Euclidean distance from points to centroid)
        max_intra_variance = 0.0
        for cluster_id in range(k):
            mask = labels == cluster_id
            if mask.sum() == 0:
                continue
            cluster_points = X[mask]
            distances = np.linalg.norm(cluster_points - centroids[cluster_id], axis=1)
            cluster_variance = np.mean(distances)
            max_intra_variance = max(max_intra_variance, cluster_variance)

        # Compute min inter-cluster distance between centroids
        min_inter_distance = float("inf")
        if k > 1:
            for i in range(k):
                for j in range(i + 1, k):
                    d = np.linalg.norm(centroids[i] - centroids[j])
                    min_inter_distance = min(min_inter_distance, d)
        else:
            min_inter_distance = 0.0

        if verbose:
            log.info(
                f"KMeans k={k}: max_intra_variance={max_intra_variance:.2f}, "
                f"min_inter_distance={min_inter_distance:.2f}"
            )

        # For k=1, accept if variance is low (homogeneous pool)
        # Apply a stricter threshold (half) so that a single cluster is only
        # accepted when the reads are truly homogeneous.
        if k == 1:
            if max_intra_variance <= variance_threshold / 2.0:
                chosen_k = 1
                chosen_labels = labels
                break
            # variance too high with k=1 -> try splitting
            continue

        # For k>1, require low variance AND high inter-cluster separation
        if (
            max_intra_variance <= variance_threshold
            and min_inter_distance >= distance_threshold
        ):
            chosen_k = k
            chosen_labels = labels
            break

    if chosen_k is None or chosen_labels is None:
        log.debug(
            f"No suitable KMeans clustering found for k=1..{max_k}. Returning None."
        )
        return None

    if verbose:
        log.info(f"Chosen KMeans clustering with k={chosen_k}")
        for i, rn in enumerate(readnames):
            log.info(
                f"  {rn}: cluster={chosen_labels[i]}, "
                f"ins={dict_summed_indels[rn][0]}, del={dict_summed_indels[rn][1]}"
            )

    # 2) Assemble consensus for each cluster
    result: dict[str, consensus_class.Consensus] | None = None
    failed_reads: list[str] = []  # list(size_outliers)

    try:
        with tempfile.TemporaryDirectory(
            dir=tmp_dir_path, delete=False if tmp_dir_path else True
        ) as tmp_dir:
            for cluster_id in range(chosen_k):
                chosen_reads = [
                    rn
                    for rn, label in zip(readnames, chosen_labels, strict=True)
                    if label == cluster_id
                ]
                if len(chosen_reads) == 0:
                    continue

                # Filter to reads that are actually in the pool
                chosen_reads = [rn for rn in chosen_reads if rn in pool]
                if len(chosen_reads) == 0:
                    continue

                consensus_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix=f"consensus.kmeans.{cluster_id}.",
                    suffix=".fasta",
                    delete=False if tmp_dir_path else True,
                )
                reads_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix=f"reads.kmeans.{cluster_id}.",
                    suffix=".fasta",
                    delete=False if tmp_dir_path else True,
                )
                with open(reads_fasta.name, "w") as f:
                    SeqIO.write(
                        [pool[readname] for readname in chosen_reads], f, "fasta"
                    )
                consensus_name = f"{min(crIDs)}.{cluster_id}"

                consensus_sequence: str | None = assemble_consensus(
                    lamassemble_mat=lamassemble_mat,
                    name=consensus_name,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    threads=threads,
                    timeout=timeout,
                    method=consensus_method,
                    ava_bam=None,
                    representative_read=None,
                    verbose=verbose,
                )

                if consensus_sequence is None:
                    log.warning(
                        f"consensus assembly failed for KMeans cluster {cluster_id} "
                        f"with {len(chosen_reads)} reads. Skipping this cluster."
                    )
                    failed_reads.extend(chosen_reads)
                    continue

                if result is None:
                    result = {}
                consensus = final_consensus(
                    samplename=samplename,
                    min_indel_size=8,
                    min_bnd_size=100,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    consensus_sequence=consensus_sequence,
                    ID=consensus_name,
                    crIDs=crIDs,
                    original_regions=[
                        (cr.chr, cr.referenceStart, cr.referenceEnd)
                        for cr in candidate_regions.values()
                    ],
                    threads=threads,
                    verbose=verbose,
                    tmp_dir_path=Path(tmp_dir),
                )
                if consensus is not None:
                    result[consensus_name] = consensus
                else:
                    log.warning(
                        f"final consensus generation failed for KMeans cluster {cluster_id} "
                        f"with {len(chosen_reads)} reads. Skipping this cluster."
                    )
                    failed_reads.extend(chosen_reads)
                    continue

            # Try to rescue failed reads by assembling them into one consensus
            if len(failed_reads) > 0 and (result is None or len(result) == 0):
                log.info(
                    f"Trying to rescue {len(failed_reads)} failed reads from KMeans "
                    f"clustering by creating one consensus from them."
                )
                consensus_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="consensus.kmeans.rescue.",
                    suffix=".fasta",
                    delete=False if tmp_dir_path else True,
                )
                reads_fasta = tempfile.NamedTemporaryFile(
                    dir=tmp_dir,
                    prefix="reads.kmeans.rescue.",
                    suffix=".fasta",
                    delete=False if tmp_dir_path else True,
                )
                rescue_reads = [rn for rn in failed_reads if rn in pool]
                with open(reads_fasta.name, "w") as f:
                    SeqIO.write(
                        [pool[readname] for readname in rescue_reads], f, "fasta"
                    )
                consensus_name = f"{min(crIDs)}.rescue"

                consensus_sequence = assemble_consensus(
                    lamassemble_mat=lamassemble_mat,
                    name=consensus_name,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    threads=threads,
                    timeout=timeout,
                    method=consensus_method,
                    ava_bam=None,
                    representative_read=None,
                    verbose=verbose,
                )

                if consensus_sequence is not None:
                    if result is None:
                        result = {}
                    consensus = final_consensus(
                        samplename=samplename,
                        min_indel_size=8,
                        min_bnd_size=100,
                        reads_fasta=Path(reads_fasta.name),
                        consensus_fasta_path=Path(consensus_fasta.name),
                        consensus_sequence=consensus_sequence,
                        ID=consensus_name,
                        crIDs=crIDs,
                        original_regions=[
                            (cr.chr, cr.referenceStart, cr.referenceEnd)
                            for cr in candidate_regions.values()
                        ],
                        threads=threads,
                        verbose=verbose,
                        tmp_dir_path=Path(tmp_dir),
                    )
                    if consensus is not None:
                        result[consensus_name] = consensus
                    else:
                        log.warning(
                            f"final consensus generation failed for rescue reads "
                            f"with {len(rescue_reads)} reads. Returning None."
                        )
                else:
                    log.warning(
                        f"consensus assembly failed for rescue reads "
                        f"with {len(rescue_reads)} reads. Returning None."
                    )

    except Exception as e:
        log.warning(
            f"consensus_while_clustering_with_kmeans failed with exception: {e}. Returning None."
        )
        return None

    return result


def add_unaligned_reads_to_consensuses_inplace(
    samplename: str,
    consensus_objects: dict[str, consensus_class.Consensus],
    pool: dict[str, SeqRecord],
    min_indel_size: int,
    min_bnd_size: int,
) -> None:
    """Inplace: Add all reads that are not used in any consensus object to the consensus object that they align best to."""
    """A read is not added to a consensus if its SV signals are significantly more than all already added reads."""
    assert len(consensus_objects) > 0, (
        "consensus_sequences must contain at least one consensus object"
    )
    if len(pool) == 0:
        return consensus_objects
    # write consensus sequences to a fasta file
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_consensus = tempfile.NamedTemporaryFile(
            prefix="tmp_consensus.", suffix=".fasta", dir=tmp_dir, delete=True
        )
        tmp_reads = tempfile.NamedTemporaryFile(
            prefix="tmp_reads.", suffix=".fasta", dir=tmp_dir, delete=True
        )
        with open(tmp_consensus.name, "w") as f:
            for consensus in consensus_objects.values():
                SeqIO.write(
                    SeqRecord(
                        name=consensus.ID,
                        id=consensus.ID,
                        seq=Seq(consensus.consensus_sequence),
                    ),
                    f,
                    "fasta",
                )
        with open(tmp_reads.name, "w") as f:
            SeqIO.write(pool.values(), f, "fasta")
        # align all reads to all consensus sequences
        tmp_alignments = tempfile.NamedTemporaryFile(
            prefix="tmp_alignments.", suffix=".bam", dir=tmp_dir, delete=True
        )
        util.align_reads_with_minimap(
            reference=tmp_consensus.name,
            reads=tmp_reads.name,
            bamout=tmp_alignments.name,
            threads=1,
            tech="map-ont",
            aln_args="  --secondary=no --sam-hit-only -U 25,75 -H",
        )

        # for each alignment, choose the primary alignment's reference name to re-distribute the read to the consensus object
        redistribution: dict[str, str] = {}  # {readname:consensusID}
        with pysam.AlignmentFile(tmp_alignments.name, mode="r") as aln_file:
            alns = list(aln_file)
        for aln in alns:
            if aln.is_unmapped:
                continue
            if aln.is_supplementary:
                continue
            readname = aln.query_name
            consensusID = aln.reference_name
            if readname in redistribution:
                log.warning(
                    f"read {readname} already assigned to consensus {redistribution[readname]}. Skipping additional alignment to {consensusID}."
                )
                continue
            redistribution[readname] = consensusID

        for consensusID, consensus in consensus_objects.items():
            assigned_alignments = [
                aln
                for aln in alns
                if aln.query_name in redistribution
                and redistribution[aln.query_name] == consensusID
            ]
            if len(assigned_alignments) == 0:
                continue
            add_cutread_alignments_to_consensus_inplace(
                alns=assigned_alignments,
                consensus=consensus,
                min_indel_size=min_indel_size,
                min_bnd_size=min_bnd_size,
                samplename=samplename,
            )


# =============================================================================
# SIGNAL EXTRACTION AND SCORING ALGORITHMS
# =============================================================================


def extract_signals(rass: list[datatypes.ReadAlignmentSignals]) -> np.ndarray:
    """Create a numpy array with columns: ref_pos, ref_end, size, sv_type, k, readname_index"""
    readnames = {ras.read_name: i for i, ras in enumerate(rass)}
    signals = np.array(
        [
            np.array(
                [
                    sv.ref_start,
                    sv.ref_end,
                    sv.size,
                    sv.sv_type,
                    k,
                    readnames[ras.read_name],
                ],
                dtype=np.int32,
            )
            for k, ras in enumerate(rass)
            for sv in ras.SV_signals
        ],
        dtype=np.int32,
    )
    if signals.size == 0:
        return np.array([])
    signals = signals[signals[:, 0].argsort()]
    return signals


def closeness_function(x, k: int = 10):
    """Calculate closeness function for SV signal proximity scoring."""
    v = abs(x)
    maxrange = 1000
    cutoff = 0.05
    if v <= k:
        return 1.0
    if v > maxrange:
        return 0.0
    else:
        if np.log2(v - k) - cutoff <= 0.0:
            return 1.0
        else:
            return 1 / np.log2(v - k) - cutoff


def probability_of_sv_presence_with_neighbors(
    svSignal: datatypes.SVsignal,
    neighbors: np.ndarray,
    N_samples: int,
    size_factor: float = 0.2,
    radius: int = 10,
) -> float:
    # a high prob is given, if many samples are represented in the neighbors
    # and the sv is very close to its neighbors
    # and the size of the signal is very similar to the neighbors (for insertions and deletions)
    # check if neighbors is empty. If so, return 0.0
    if neighbors.shape[0] == 0:
        return 0.0
    if svSignal.sv_type == 0:  # insertion
        margin = max(5, size_factor * svSignal.size)
        mask_similar_size = abs(neighbors[:, 2] - svSignal.size) <= margin
        # aplly distance function to score close signals higher than distant signals
        # pick one signal per sample that fits best
        neighbors_selected = neighbors[mask_similar_size]
        proximities = np.array([
            closeness_function(x, radius)
            for x in neighbors_selected[:, 0] - svSignal.ref_start
        ])
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(neighbors_selected[:, 5]):
            mask_sample = neighbors_selected[:, 5] == rn_ID
            sum_signal += np.max(proximities[mask_sample])
        return sum_signal / N_samples
    if svSignal.sv_type == 1:  # deletion
        # similar to insertion, but the signal is split into two parts
        # the left and the right side of the deletion
        margin = max(5, size_factor * svSignal.size)
        mask_similar_size = abs(neighbors[:, 2] - svSignal.size) <= margin
        # aplly distance function to score close signals higher than distant signals
        # pick one signal per sample that fits best
        neighbors_selected = neighbors[mask_similar_size]
        proximities_left = np.array([
            closeness_function(x, radius)
            for x in neighbors_selected[:, 0] - svSignal.ref_start
        ])
        proximities_right = np.array([
            closeness_function(x, radius)
            for x in neighbors_selected[:, 1] - (svSignal.ref_start + svSignal.size)
        ])
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(neighbors_selected[:, 5]):
            mask_sample = neighbors_selected[:, 5] == rn_ID
            sum_signal += np.max(proximities_left[mask_sample]) * 0.5
            sum_signal += np.max(proximities_right[mask_sample]) * 0.5
        return sum_signal / N_samples
    if svSignal.sv_type == 3 or svSignal.sv_type == 4:  # break end
        proximities = np.array([
            closeness_function(x, radius) for x in neighbors[:, 0] - svSignal.ref_start
        ])
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(neighbors[:, 5]):
            mask_sample = neighbors[:, 5] == rn_ID
            sum_signal += np.max(proximities[mask_sample])
        return sum_signal / N_samples
    raise ValueError(
        "sv must be of type AlignmentInsertion, AlignmentDeletion or AlignmentBreakend."
    )


def probability_of_sv_with_similar_svs(
    svSignal: datatypes.SVsignal,
    signals: np.ndarray,
    N_samples: int,
    radius: int = 10_000,
    size_tolerance: float = 0.05,
) -> float:
    # for each sample the best matching row in signals is found (matching by type and size)
    # the probability is the number of samples that have a similar signal in the vicinity.
    # It can't be calculated for breakends
    # works only fo SVs that are greater than a given size
    # check if signals is empty. If so, return 0.0
    if signals.shape[0] == 0:
        return 0.0
    if svSignal.sv_type == 0:  # insertion
        margin = max(5, size_tolerance * svSignal.size)
        # subset signals to signals of the same type and signals that max. differ in size by margin
        mask = (signals[:, 3] == 0) & (abs(signals[:, 2] - svSignal.size) <= margin)
        candidates = signals[mask]
        # apply closeness_function to the size difference of the candidates
        similarities = np.array([
            closeness_function(x, radius) for x in candidates[:, 2] - svSignal.size
        ])
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(candidates[:, 5]):
            mask_sample = candidates[:, 5] == rn_ID
            sum_signal += np.max(similarities[mask_sample])
        return sum_signal / N_samples
    if svSignal.sv_type == 1:  # deletion
        margin = max(5, size_tolerance * svSignal.size)
        # subset signals to signals of the same type and signals that max. differ in size by margin
        mask = (signals[:, 3] == 1) & (abs(signals[:, 2] - svSignal.size) <= margin)
        candidates = signals[mask]
        # apply closeness_function to the size difference of the candidates
        similarities = np.array([
            closeness_function(x, radius) for x in candidates[:, 2] - svSignal.size
        ])
        # pick the best matching signal for each sample
        sum_signal = 0.0
        for rn_ID in np.unique(candidates[:, 5]):
            mask_sample = candidates[:, 5] == rn_ID
            sum_signal += np.max(similarities[mask_sample])
        return sum_signal / N_samples
    raise ValueError(
        "sv must be of type SVsignal and must have sv_type 0 (insertion) or 1 (deletion)."
    )


def sv_neighbors(
    svSignal: datatypes.SVsignal, signals: np.ndarray, radius: int
) -> np.ndarray:
    # check if signals is empty. If so, return an empty array
    if signals.size == 0:
        return np.array([])
    if svSignal.sv_type == 0:  # insertion
        # check if there are signals of the same type in the vicinity
        start, end = svSignal.ref_start - radius, svSignal.ref_start + radius
        return signals[
            (signals[:, 3] == 0) & (signals[:, 0] >= start) & (signals[:, 0] <= end)
        ]
    if svSignal.sv_type == 1:  # deletion
        # check if there are signals of the same type in the vicinity
        start_left, end_left = svSignal.ref_start - radius, svSignal.ref_start + radius
        start_right, end_right = svSignal.ref_end - radius, svSignal.ref_end + radius
        mask_left = (
            (signals[:, 3] == 1)
            & (signals[:, 0] >= start_left)
            & (signals[:, 0] <= end_left)
        )
        mask_right = (
            (signals[:, 3] == 1)
            & (signals[:, 0] >= start_right)
            & (signals[:, 0] <= end_right)
        )
        return signals[mask_left | mask_right]
    if svSignal.sv_type == 3 or svSignal.sv_type == 4:  # break end
        # check if there are signals of the same type in the vicinity
        start, end = svSignal.ref_start - radius, svSignal.ref_start + radius
        return signals[
            ((signals[:, 3] == 3) | (signals[:, 3] == 4))
            & (signals[:, 0] >= start)
            & (signals[:, 0] <= end)
        ]
    raise ValueError("sv must be of type SVsignal and must have sv_type 0, 1, 3 or 4.")


# function to calculate both scores from probability_of_sv_presence_with_neighbors and probability_of_sv_with_similar_svs
# and picks the maximum of both
def probability_of_sv(
    svSignal: datatypes.SVsignal,
    signals: np.ndarray,
    N_samples: int,
    bias_for_locality: float = 0.5,
    size_factor: float = 0.2,
    size_tolerance: float = 0.05,
    radius_close: int = 30,
    radius_neighbors: int = 1000,
    radius_far: int = 10_000,
) -> float:
    """Computes the probability of the presence of a structural variant.
    Args:
        sv (_type_): structural variant of the type AlignmentInsertion, AlignmentDeletion or AlignmentBreakend
        signals (np.ndarray): array of signals with columns: [position_start,position_end,size,type,raf_index,read_index]
        N_samples (int): number of samples in all rafs
        size_factor (float, optional): The maximum difference in size of two SV signals when matching locally. Defaults to 0.05.
        radius_close (int, optional): radius to consider SV signals in locality sensitive approximation. Defaults to 10.
        radius_far (int, optional): radius to consider SV signals in similarity sensitive approximation. Defaults to 10_000.

    Returns:
        float: probability of the presence of the SV
    """
    # check input. is signals empty? is N_samples 0?
    if signals.size == 0:
        return 0.0
    if N_samples == 0:
        return 0.0
    neighbors = sv_neighbors(
        svSignal=svSignal, radius=radius_neighbors, signals=signals
    )
    score_locality = probability_of_sv_presence_with_neighbors(
        svSignal=svSignal,
        neighbors=neighbors,
        N_samples=N_samples,
        size_factor=size_factor,
        radius=radius_close,
    )
    # can compute similarity score only for insertions and deletions
    if svSignal.sv_type == 3 or svSignal.sv_type == 4:
        return score_locality
    score_similarity = probability_of_sv_with_similar_svs(
        svSignal=svSignal,
        signals=signals,
        N_samples=N_samples,
        size_tolerance=size_tolerance,
        radius=radius_far,
    )
    # weighted mean
    return (
        bias_for_locality * score_locality
        + (1.0 - bias_for_locality) * score_similarity
    )


def calc_ras_scores(
    rass: list[datatypes.ReadAlignmentSignals],
    reflen: int,
    size_factor: float = 0.05,
    radius_close: int = 15,
    radius_neighbors: int = 100,
    radius_far: int = 10_000,
    bias: float = 0.5,
) -> np.ndarray:
    # calc a score for each raf. The score of a raf is the sum of probabilities of each sv signal in the raf
    N_samples = len({ras.read_name for ras in rass})
    signals = extract_signals(rass=rass)
    # generate regions with repeats / low complexity
    # For each
    ras_scores = np.zeros(len(rass))
    for i, ras in enumerate(rass):
        for svSignal in ras.SV_signals:
            if svSignal.sv_type <= 2:
                weight = svSignal.size
            elif svSignal.sv_type == 3:
                # get clipped length
                clipped_length = svSignal.size
                clipped_length_on_ref = svSignal.ref_start
                weight = max(1000, min(clipped_length, clipped_length_on_ref))
            elif svSignal.sv_type == 4:
                # get clipped length
                clipped_length = svSignal.size
                clipped_length_on_ref = reflen - svSignal.ref_end
                weight = max(1000, min(clipped_length, clipped_length_on_ref))
            else:
                raise f"svSignal.sv_type {svSignal.sv_type} not recognized. expected one of 0,1,2,3,4."
            ras_scores[i] += (
                probability_of_sv(
                    svSignal=svSignal,
                    signals=signals,
                    N_samples=N_samples,
                    bias_for_locality=bias,
                    size_factor=size_factor,
                    radius_close=radius_close,
                    radius_neighbors=radius_neighbors,
                    radius_far=radius_far,
                )
                * weight
            )
    return ras_scores


def parse_ReadAlignmentSignals_from_alignment(
    samplename: str,
    alignment: pysam.AlignedSegment,
    min_signal_size: int,
    min_bnd_size: int,
) -> datatypes.ReadAlignmentSignals:
    ref_start, ref_end, read_start, read_end = alignments_to_rafs.get_start_end(
        alignment
    )
    sv_signals = alignments_to_rafs.parse_SVsignals_from_alignment(
        alignment=alignment,
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        min_signal_size=min_signal_size,
        min_bnd_size=min_bnd_size,
    )
    return datatypes.ReadAlignmentSignals(
        samplename=samplename,
        read_name=str(alignment.query_name),
        reference_name=str(alignment.reference_name),
        alignment_forward=bool(not alignment.is_reverse),
        SV_signals=sorted(sv_signals),
    )


def filter_excessive_QC_ras_inplace(
    rass: list[datatypes.ReadAlignmentSignals],
    read_sequences: dict[str, SeqRecord],
    reference_sequence: str,
) -> None:
    _loss_logger = get_signal_loss_logger()
    for ras in rass:
        filtered_signals = []
        for sv_signal in ras.SV_signals:
            if sv_signal.sv_type == 0:  # insertion (query sequence)
                dna_string = read_sequences[ras.read_name][
                    sv_signal.read_start : sv_signal.read_end
                ].seq
            elif sv_signal.sv_type == 1:  # deletion (reference sequence)
                dna_string = reference_sequence[sv_signal.ref_start : sv_signal.ref_end]
            else:
                filtered_signals.append(sv_signal)
                continue
            gc_frac = SeqUtils.gc_fraction(Seq(dna_string))
            if 0.1 < gc_frac <= 0.9:
                filtered_signals.append(sv_signal)
            else:
                _loss_logger.log_filtered(
                    stage="filter_excessive_QC_ras_inplace",
                    reason="extreme_GC_fraction",
                    details={
                        "read_name": ras.read_name,
                        "sv_type": sv_signal.sv_type,
                        "size": sv_signal.size,
                        "ref_start": sv_signal.ref_start,
                        "ref_end": sv_signal.ref_end,
                        "gc_fraction": round(gc_frac, 4),
                    },
                )
        _n_removed = len(ras.SV_signals) - len(filtered_signals)
        if _n_removed > 0:
            log.info(
                f"filter_excessive_QC_ras_inplace: read {ras.read_name}: "
                f"filtered {_n_removed}/{len(ras.SV_signals)} signals with extreme GC content"
            )
        ras.SV_signals = filtered_signals


def score_ras_from_alignments(
    samplename: str,
    alignments: list[pysam.AlignedSegment],
    read_sequences: dict[str, SeqRecord],
    reference_sequence: str,
    min_signal_size: int,
    min_bnd_size: int,
    reflen: int,
    verbose: bool = False,
) -> dict[str, float]:
    # score alignments
    # get rafs
    rass = [
        parse_ReadAlignmentSignals_from_alignment(
            samplename=samplename,
            alignment=a,
            min_bnd_size=min_bnd_size,
            min_signal_size=min_signal_size,
        )
        for a in alignments
    ]
    # filter rass sv signals for signals that have less or equal to 3 different 2mers
    filter_excessive_QC_ras_inplace(
        rass=rass, read_sequences=read_sequences, reference_sequence=reference_sequence
    )

    # score rafs (read_name might occur multiple times)
    # re-write raf_scores. It should first create statistics over all sv signal
    # then generate weights for each sv signal, based on a maximum score
    # of 1) a distance score and 2) a size score.

    # handicap = get_handicap(rass=rass,handicap_cutoff=0.9)

    # return dict of {read_name: score}
    # always pick the highest scoring raf for each read_name
    rass_scores = calc_ras_scores(rass=rass, reflen=reflen)
    rafs_scores_dict = {}
    for ras, score in zip(rass, rass_scores, strict=True):
        if ras.read_name in rafs_scores_dict:
            if rafs_scores_dict[ras.read_name] > score:
                rafs_scores_dict[ras.read_name] = score
        else:
            rafs_scores_dict[ras.read_name] = float(score)

    # log.debug(f"done scoring rafs with {len(set([ras.read_name for ras in rass]))} reads and {len(rass)} ReadAlignmentSignals objects.")
    # log.debug(f"Resulting: {rafs_scores_dict}.")
    # log.debug(f"handicap: {handicap}")
    if verbose:
        sorted_scores = sorted(rafs_scores_dict.items(), key=lambda x: x[1])
        print_str = "\n".join([f"{read}:\t{score}" for read, score in sorted_scores])
        log.info(f"scoring results:\n{print_str}")
    return rafs_scores_dict  # ,handicap


# =======================================================================================================================================================
#  consensus padding
# =======================================================================================================================================================


def parse_description(description: str) -> dict:
    result = {}
    for key_value in description.split(","):
        key, value = key_value.split("=")
        if "chr" in key:
            value = value
        else:
            result[key] = int(value)
    return result


def create_padding_for_consensus(
    consensus_object: consensus_class.Consensus,
    cutreads: dict[str, SeqRecord],
    read_records: dict[str, SeqRecord],
) -> consensus_class.ConsensusPadding:
    """Creates a ConsensusPadding object with the padding sequence, the consensus interval on the padded squence, and the read names."""

    read_paddings_for_consensus = _get_all_read_padding_intervals(
        consensus_object=consensus_object, cutreads=cutreads, read_records=read_records
    )

    padding_sizes_per_read = _get_padding_sizes_per_read(
        read_paddings_for_consensus=read_paddings_for_consensus
    )

    padding_reads = _get_padding_read_names_of_consensus(
        padding_sizes_per_read=padding_sizes_per_read
    )

    padding_intervals = _get_read_padding_intervals(
        read_paddings_for_consensus=read_paddings_for_consensus,
        padding_read_names_of_consensus=padding_reads,
    )
    padding = _create_padding_object(
        cons=consensus_object,
        read_paddings_for_consensus=read_paddings_for_consensus,
        padding_reads=padding_reads,
        padding_intervals=padding_intervals,
        read_records=read_records,
    )
    return padding


def _get_all_read_padding_intervals(
    consensus_object: consensus_class.Consensus,
    cutreads: dict[str, SeqRecord],
    read_records: dict[str, SeqRecord],
) -> dict[str, tuple[int, int, int, int, bool]]:
    """Generates a dict[readname, (min start, max end, cutread length, read length, is_forward)] for all reads in the consensus object."""
    directions: dict[str, bool] = {
        readname: forward
        for start, end, readname, forward in (
            consensus_object.intervals_cutread_alignments
        )
    }

    data: dict[str, tuple[int, int]] = {}
    for cutread in cutreads.values():
        if cutread.name not in directions.keys():
            continue
        cutread_description_dict = parse_description(cutread.description)
        if cutread.name not in data:
            data[cutread.name] = (
                cutread_description_dict["start"],
                cutread_description_dict["end"],
            )
        else:
            data[cutread.name] = (
                min(data[cutread.name][0], cutread_description_dict["start"]),
                max(data[cutread.name][1], cutread_description_dict["end"]),
            )
    result: dict[str, tuple[int, int, int, int, bool]] = {}
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


def _get_padding_sizes_per_read(
    read_paddings_for_consensus: dict[str, tuple[int, int, int, int, bool]],
) -> dict[str, tuple[int, int]]:
    """Generates a dict[readname, (left padding size, right padding size)] for all reads in the consensus object."""
    """A padding is the sequence that overshoots the cutread alignment on the left and right side (in the orientation of the consensus)."""
    result: dict[str, tuple[int, int]] = {}
    for readname, (
        start,
        end,
        _cutread_length,
        read_length,
        forward,
    ) in read_paddings_for_consensus.items():
        alpha: int = start if forward else read_length - end
        beta: int = read_length - end if forward else start
        result[readname] = (alpha, beta)
    return result


def _get_padding_read_names_of_consensus(
    padding_sizes_per_read: dict[str, tuple[int, int]],
) -> tuple[str, str]:
    """Returns the read names of the left and right padding reads."""
    return (
        max(padding_sizes_per_read.items(), key=lambda x: x[1][0])[0],
        max(padding_sizes_per_read.items(), key=lambda x: x[1][1])[0],
    )


def _get_read_padding_intervals(
    read_paddings_for_consensus: dict[str, tuple[int, int, int, int, bool]],
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


def _create_padding_object(
    cons: consensus_class.Consensus,
    read_paddings_for_consensus: dict[str, tuple[int, int, int, int, bool]],
    padding_reads: tuple[str, str],
    padding_intervals: tuple[tuple[int, int], tuple[int, int]],
    read_records: dict[str, SeqRecord],
) -> consensus_class.ConsensusPadding:
    """Creates a new ConsensusPadding object with the padding sequence and the read names."""

    left_padding: Seq = read_records[padding_reads[0]].seq[
        padding_intervals[0][0] : padding_intervals[0][1]
    ]
    # needs to be reverse complimented if the cutread alignment is reverse
    if not read_paddings_for_consensus[padding_reads[0]][4]:
        left_padding = left_padding.reverse_complement()

    right_padding: Seq = read_records[padding_reads[1]].seq[
        padding_intervals[1][0] : padding_intervals[1][1]
    ]
    # needs to be reverse complimented if the cutread alignment is reverse
    if not read_paddings_for_consensus[padding_reads[1]][4]:
        right_padding = right_padding.reverse_complement()
    # add padding to the consensus sequence
    padded_consensus_sequence = Seq(
        left_padding.lower() + cons.consensus_sequence.upper() + right_padding.lower()
    )

    new_padding = consensus_class.ConsensusPadding(
        sequence=str(padded_consensus_sequence),
        readname_left=padding_reads[0],
        readname_right=padding_reads[1],
        padding_size_left=len(left_padding),
        padding_size_right=len(right_padding),
        consensus_interval_on_sequence_with_padding=(
            len(left_padding),
            len(padded_consensus_sequence) - len(right_padding),
        ),
    )
    return new_padding


# =============================================================================
# DATABASE OPERATIONS AND SERIALIZATION
# =============================================================================


def deserialize_crs_container(crs_container: dict) -> dict:
    crs = [
        cattrs.structure(cr, datatypes.CandidateRegion) for cr in crs_container["crs"]
    ]
    connecting_reads = crs_container["connecting_reads"]
    return {"crs": crs, "connecting_reads": connecting_reads}


def load_crs_containers_from_db(
    path_db: Path, crIDs: list[int] | None = None
) -> dict[int, dict]:
    """Load candidate region containers from database."""
    # if crIDs is not None, then check if all elements are of type int
    # check if path_db is a valid path
    if not path_db.exists():
        raise FileNotFoundError(f"Database file {path_db} does not exist.")
    if not path_db.is_file():
        raise ValueError(f"{path_db} is not a file.")
    if crIDs is not None:
        assert type(crIDs) == list, "crIDs must be a list of integers."
        assert len(crIDs) > 0, (
            "crIDs must be a list of integers with at least one element."
        )
        if not all(isinstance(crID, int) for crID in crIDs):
            raise ValueError(f"crIDs must be a list of integers. Instead got {crIDs}")

    conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
    c = conn.cursor()
    if crIDs is not None:
        c.execute(
            """SELECT crID, data FROM containers WHERE crID IN ({0})""".format(
                ", ".join([str(crID) for crID in crIDs])
            )
        )
    else:
        c.execute("SELECT crID, data FROM containers")
    containers = {}
    for crID, data in c.fetchall():
        crs_container = json.loads(data)
        crs_container = deserialize_crs_container(crs_container)
        containers[crID] = crs_container
    c.close()
    conn.close()
    return containers


def summed_indel_distribution(
    alns: dict[int, list[pysam.AlignedSegment]],
    crs: dict[int, datatypes.CandidateRegion],
) -> dict[str, list[int]]:
    # Build an IntervalTree from the candidate regions so we can quickly check
    # whether a signal falls within any provided region.
    cr_tree = IntervalTree()
    for cr in crs.values():
        if cr.referenceStart < cr.referenceEnd:
            cr_tree.addi(cr.referenceStart, cr.referenceEnd)

    rafs: list[datatypes.ReadAlignmentFragment] = []
    for _crID, alnlist in alns.items():
        for aln in alnlist:
            rafs.append(
                alignments_to_rafs.parse_ReadAlignmentFragment_from_alignment(
                    alignment=aln,
                    samplename="None",
                    min_signal_size=12,
                    min_bnd_size=100,
                    filter_density_radius=500,
                    filter_density_min_bp=30,
                )
            )

    # create a dict of the form read:sum_signals ({str,int})
    read_sum_signals: dict[str, int] = {
        aln.query_name: [0, 0]
        for aln in [aln for alnlist in alns.values() for aln in alnlist]
    }
    # now add all summed insertions, if they exist
    # only keep signals whose reference position overlaps a candidate region interval
    for raf in rafs:
        insertions = [
            signal.size
            for signal in raf.SV_signals
            if signal.sv_type == 0
            and cr_tree.overlaps(
                signal.ref_start, max(signal.ref_start + 1, signal.ref_end)
            )
        ]
        deletions = [
            signal.size
            for signal in raf.SV_signals
            if signal.sv_type == 1
            and cr_tree.overlaps(
                signal.ref_start, max(signal.ref_start + 1, signal.ref_end)
            )
        ]
        if len(insertions) > 0:
            read_sum_signals[raf.read_name][0] += sum(insertions)
        if len(deletions) > 0:
            read_sum_signals[raf.read_name][1] += sum(deletions)
    return read_sum_signals


def print_indel_distribution(read_sum_signals: dict[str, list[int]]) -> None:
    for readname, (sum_inss, sum_dels) in sorted(
        read_sum_signals.items(), key=lambda x: x[1]
    ):
        print(f"{readname}: {sum_inss} insertions, {sum_dels} deletions")


def process_consensus_container(
    samplename: str,
    crs_dict: dict[int, datatypes.CandidateRegion],
    path_alignments: Path,
    copy_number_tracks: Path,
    lamassemble_mat: Path | str | None,
    timeout: int,
    buffer_clipped_length: int,
    consensus_method: str,
    threads: int = 1,
    tmp_dir_path: Path | str | None = None,
    figures_dir: Path | None = None,
    verbose: bool = False,
    densities_weight: float = 1.0,
    max_intra_distance: float = -1.0,
    cn_override: int | None = None,
) -> tuple[
    dict[str, consensus_class.Consensus], dict[int, list[datatypes.SequenceObject]]
]:
    # alns: dict[crID:list[AlignedSegment]]
    # compute max_copy_number from the intervals of the candidate regions and by querying the bgzipped and tabix indexed copy number tracks
    regions_for_cn_query = [
        (cr.chr, cr.referenceStart, cr.referenceEnd) for cr in crs_dict.values()
    ]
    if verbose:
        for cr in crs_dict.values():
            print(f"CR region on ref: {cr.chr}:{cr.referenceStart}-{cr.referenceEnd}")
    if cn_override is not None:
        max_copy_number = cn_override
        log.info(f"Using overridden copy number: {max_copy_number}")
    else:
        max_copy_number = max(
            2,
            copynumber_tracks.query_copynumber_from_regions(
                bgzip_bed=copy_number_tracks, regions=regions_for_cn_query
            ),
        )  # 2 minimum clusters - maybe this is really bad, idk.
        log.info(f"Maximum copy number for this container: {max_copy_number}")
    if max_copy_number > 4:
        # don't process this container. Too complex.
        log.warning(
            f"Maximum copy number {max_copy_number} exceeds threshold of 4. Skipping consensus building for this container."
        )
        # parse unused reads to datatypes.SequenceObject
        dict_unused_read_names: dict[int, set[str]] = {
            cr.crID: {signal.readname for signal in cr.sv_signals}
            for cr in crs_dict.values()
        }
        dict_unused_reads: dict[int, list[datatypes.SequenceObject]] = {
            crID: [] for crID in dict_unused_read_names
        }
        for crID, readnames in dict_unused_read_names.items():
            for readname in readnames:
                # get the read sequence from the alignment file
                with pysam.AlignmentFile(str(path_alignments), "rb") as samfile:
                    read_seqRecord = samfile.fetch(contig=readname)
                    for read in read_seqRecord:
                        read_sequence_object = datatypes.SequenceObject(
                            id=read.query_name,
                            name=read.query_name,
                            sequence=read.query_sequence,
                            description=read.query_name,
                            qualities=(
                                read.query_qualities
                                if read.query_qualities is not None
                                else None
                            ),
                        )
                        dict_unused_reads[crID].append(read_sequence_object)
        return {}, dict_unused_reads

    alns, _alns_wt = get_read_alignments_for_crs(
        crs=list(crs_dict.values()), alignments=path_alignments
    )
    if verbose:
        # print the qname and reference intervals of the alignments in alns for each alignment
        for crID, alnlist in alns.items():
            print(f"CR {crID} has {len(alnlist)} alignments:")
            for aln in alnlist:
                # either left or right break end
                clipped_left = (
                    aln.cigartuples[0][1] if aln.cigartuples[0][0] == 4 else 0
                )
                clipped_right = (
                    aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 4 else 0
                )
                print(
                    f"\t{aln.query_name}: {aln.reference_name}:{aln.reference_start}-{aln.reference_end}, clipped: left={clipped_left}, right={clipped_right}"
                )
    dict_all_intervals = get_read_alignment_intervals_in_cr(
        crs=list(crs_dict.values()),
        dict_alignments=alns,
        buffer_clipped_length=buffer_clipped_length,
    )
    if verbose:
        print(f"dict_all_intervals: {dict_all_intervals}")

    max_intervals = get_max_extents_of_read_alignments_on_cr(
        dict_all_intervals=dict_all_intervals
    )
    if verbose:
        print(f"max_intervals: {max_intervals}")
    read_records: dict[str, SeqRecord] = get_full_read_sequences_of_alignments(
        dict_alignments=alns, path_alignments=path_alignments
    )
    log.info("cutting reads from alignments")
    cutreads: dict[str, SeqRecord] = trim_reads(
        dict_alignments=alns, intervals=max_intervals, read_records=read_records
    )
    log.info(f"number of reads: {len(cutreads)}")
    log.info(
        f"summed trimmed reads bp: {sum(len(read.seq) for read in cutreads.values())}"
    )

    dict_summed_indels: dict[str, list[int]] = summed_indel_distribution(
        alns=alns, crs=crs_dict
    )
    if verbose:
        print_indel_distribution(dict_summed_indels)

    # =========================================== CONSENSUS BUILDING =========================================== #

    res: dict[str, consensus_class.Consensus] | None = None
    consensus_objects: dict[str, consensus_class.Consensus] = {}
    res = consensus_while_clustering_with_kmeans(
        samplename=samplename,
        dict_summed_indels=dict_summed_indels,
        lamassemble_mat=lamassemble_mat,
        pool=cutreads,
        candidate_regions=crs_dict,
        max_k=max_copy_number,
        variance_threshold=29.0,
        distance_threshold=29.0,
        threads=threads,
        tmp_dir_path=tmp_dir_path,
        timeout=timeout,
        verbose=verbose,
        consensus_method=consensus_method,
    )
    if not res:
        res = consensus_while_clustering(
            samplename=samplename,
            lamassemble_mat=lamassemble_mat,
            pool=cutreads,
            candidate_regions=crs_dict,
            partitions=max_copy_number,
            timeout=timeout,
            threads=threads,
            tmp_dir_path=tmp_dir_path,
            figures_dir=figures_dir,
            verbose=verbose,
            densities_weight=densities_weight,
            max_intra_distance=max_intra_distance,
            consensus_method=consensus_method,
        )
    if res is not None:
        consensus_objects.update(res)

    # TODO: check if the number of clusters is satisfying the expected number of alleles
    # via check_clustering_consensus_results
    # if not, then re-run consensus_while_clustering_with_racon with a higher separation_threshold (relaxation)
    # =========================================== CONSENSUS BUILDING END =========================================== #

    if len(consensus_objects) == 0:
        log.warning(
            "No consensus objects could be built from the candidate regions. Returning empty consensus dict and all reads as unused reads."
        )
        # parse unused reads to datatypes.SequenceObject
        dict_unused_read_names: dict[int, set[str]] = {
            cr.crID: {signal.readname for signal in cr.sv_signals}
            for cr in crs_dict.values()
        }
        dict_unused_reads: dict[int, list[datatypes.SequenceObject]] = {
            crID: [] for crID in dict_unused_read_names
        }
        for crID, readnames in dict_unused_read_names.items():
            for readname in readnames:
                read_seqRecord = cutreads[readname]
                read_sequence_object = datatypes.SequenceObject(
                    id=read_seqRecord.id,
                    name=read_seqRecord.name,
                    sequence=str(read_seqRecord.seq),
                    description=read_seqRecord.description,
                    qualities=(
                        read_seqRecord.letter_annotations["phred_quality"]
                        if "phred_quality" in read_seqRecord.letter_annotations
                        else None
                    ),
                )
                dict_unused_reads[crID].append(read_sequence_object)
        return {}, dict_unused_reads
        # TODO: implement assembly for regions without representatives

    # ================================ ADDING UNUSED READS TO CONSENSUS OBJECTS ================================ #
    # A last attempt is made to add unused reads to the consensus objects
    # identify all unused reads
    unused_reads: set[str] = set(cutreads.keys())
    for consensus in consensus_objects.values():
        unused_reads -= consensus.get_used_readnames()
    if verbose:
        log.info(
            f"Trying to add {len(unused_reads)} unused reads to the consensus objects."
        )
    pool_unused_reads = {readname: cutreads[readname] for readname in unused_reads}
    add_unaligned_reads_to_consensuses_inplace(
        pool=pool_unused_reads,
        samplename=samplename,
        consensus_objects=consensus_objects,
        min_bnd_size=100,
        min_indel_size=12,
    )
    if verbose:
        # again, find the unused reads
        unused_reads = set(cutreads.keys())
        for consensus in consensus_objects.values():
            unused_reads -= consensus.get_used_readnames()
        log.info(
            f"{len(unused_reads)} reads could not be added to the consensus objects:\n{'\n'.join(list(unused_reads))}"
        )

    # ================================ ADDING UNUSED READS TO CONSENSUS OBJECTS END ================================ #
    # add original SV singals (ExtendedSVsignal) to the consensus
    # each consensus has a list of reconstructible reads
    # each reconstructible read has a member alignment
    # this alignment has a member readname and function aug_name() that returns the augmented readname
    # use the augmented readname to get the signals from the cr

    # for consensus in consensus_objects.values():
    #     _readnames = [x.read_name for x in consensus.cut_read_alignment_signals]
    #     signals = [signal for cr in crs_dict.values() for signal in cr.sv_signals if signal.readname in _readnames]
    #     consensus.original_signals = signals
    set_unused_readnames = set(cutreads.keys())
    for consensus in consensus_objects.values():
        set_unused_readnames -= consensus.get_used_readnames()
    log.debug(
        f"Found {len(set_unused_readnames)} unused reads after consensus building."
    )

    unused_seqobjects: dict[str, datatypes.SequenceObject] = {}
    for readname in set_unused_readnames:
        read_seqRecord = cutreads[readname]
        unused_seqobjects[readname] = datatypes.SequenceObject(
            id=read_seqRecord.id,
            name=read_seqRecord.name,
            sequence=str(read_seqRecord.seq),
            description=read_seqRecord.description,
            qualities=(
                read_seqRecord.letter_annotations["phred_quality"]
                if "phred_quality" in read_seqRecord.letter_annotations
                else None
            ),
        )
    dict_unused_reads: dict[int, list[datatypes.SequenceObject]] = {}
    for crID, cr in crs_dict.items():
        # save all reads to dict_unused_reads that are in the cr and in set_unused_readnames
        readnames_in_cr = {signal.readname for signal in cr.sv_signals}
        dict_unused_reads[crID] = []
        for readname, seqobject in unused_seqobjects.items():
            if readname in readnames_in_cr:
                dict_unused_reads[crID].append(seqobject)

    # add padding to the consensus objects
    for consensus in consensus_objects.values():
        consensus.consensus_padding = create_padding_for_consensus(
            consensus_object=consensus, cutreads=cutreads, read_records=read_records
        )

    # for each consensus, create a dict
    return consensus_objects, dict_unused_reads


def load_crIDs_from_containers_db(path_db: Path) -> list[int]:
    conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT crID FROM containers")
    crIDs = [row[0] for row in c.fetchall()]
    c.close()
    conn.close()
    return crIDs


def crs_containers_to_consensus(
    samplename: str,
    input: Path,
    copy_number_tracks: Path,
    output: Path,
    lamassemble_mat: Path | str | None,
    path_alignments: Path,
    threads: int,
    buffer_clipped_sequence: int,
    timeout: int,
    consensus_method: str,
    crIDs: list[int] | None = None,
    tmp_dir_path: Path | str | None = None,
    figures_dir: Path | None = None,
    verbose: bool = False,
    densities_weight: float = 1.0,
    max_intra_distance: float = -1.0,
    fasta_debug_path: Path | None = None,
    cn_override: int | None = None,
) -> None:
    if lamassemble_mat is not None and not Path(lamassemble_mat).exists():
        raise FileNotFoundError(
            f"lamassemble matrix file {lamassemble_mat} does not exist."
        )
    if not input.exists():
        raise FileNotFoundError(f"Input file {input} does not exist.")
    if not input.is_file():
        raise FileNotFoundError(f"Input file {input} is not a file.")
    if not output.parent.exists():
        raise FileNotFoundError(f"Output directory {output.parent} does not exist.")
    if not output.parent.is_dir():
        raise NotADirectoryError(
            f"Output directory {output.parent} is not a directory."
        )
    if threads <= 0:
        raise ValueError("threads must be greater than 0.")
    if crIDs is not None:
        if not isinstance(crIDs, list):
            raise TypeError("crIDs must be a list of integers.")
        if len(crIDs) == 0:
            raise ValueError("crIDs must contain at least one element.")
        if not all(isinstance(crID, int) for crID in crIDs):
            raise TypeError("crIDs must be a list of integers.")
    log.info("load data")
    if crIDs:
        existing_crIDs = load_crIDs_from_containers_db(path_db=input)
        # log all crIDs that are missing in a warning
        if not all(crID in existing_crIDs for crID in crIDs):
            log.warning(
                f"crIDs {set(crIDs) - set(existing_crIDs)} are not in the database. They are skipped."
            )
        # if no real crIDs are left, return with warning
        if not any(crID in existing_crIDs for crID in crIDs):
            log.error("No crIDs to process. Exiting.")
            return
    log.info(f"Loading crs containers from database at {input}.")
    containers = load_crs_containers_from_db(path_db=input, crIDs=crIDs)
    log.info("loaded data")
    all_consensuses: dict[str, consensus_class.Consensus] = {}
    all_unused_reads: dict[int, list[datatypes.SequenceObject]] = {}
    for crID, container in containers.items():
        log.info(f"Processing crID {crID}.")
        crs_dict = {cr.crID: cr for cr in container["crs"]}
        # connecting_reads = container['connecting_reads']
        consensuses, unused_reads = process_consensus_container(
            samplename=samplename,
            crs_dict=crs_dict,
            tmp_dir_path=tmp_dir_path,
            path_alignments=path_alignments,
            copy_number_tracks=copy_number_tracks,
            threads=1,
            buffer_clipped_length=buffer_clipped_sequence,
            lamassemble_mat=lamassemble_mat,
            timeout=timeout,
            figures_dir=figures_dir,
            verbose=verbose,
            densities_weight=densities_weight,
            max_intra_distance=max_intra_distance,
            cn_override=cn_override,
            consensus_method=consensus_method,
        )
        all_unused_reads.update(unused_reads)
        if not consensuses:
            continue
        all_consensuses.update(consensuses)
    # do nothing and log a warning if no consensuses were generated
    result = consensus_class.CrsContainerResult(
        consensus_dicts=all_consensuses, unused_reads=all_unused_reads
    )

    # for debugging, check if all consensus objects in result.consensus_dicts have non empty original_regions
    for consensusID, consensus in result.consensus_dicts.items():
        if len(consensus.original_regions) == 0:
            raise ValueError(f"Consensus {consensusID} has no original regions.")

    if verbose:
        # list each consensus by name and the number of reads that were used to generate it
        log.info(
            f"A total of {len(result.consensus_dicts)} consensuses were generated."
        )
        for consensusID, consensus in result.consensus_dicts.items():
            log.info(
                f"Consensus {consensusID} was generated with {len(set(consensus.get_used_readnames()))} reads and max depth {consensus.get_max_cutread_coverage()}."
            )
        for crID, unused_reads in result.unused_reads.items():
            log.info(
                f"crID {crID} has {len(unused_reads)} unused reads: {[seqobj.id for seqobj in unused_reads]}"
            )

    if verbose and tmp_dir_path is not None:
        # write all consensus padded sequences to a fasta file in the tmp dir
        fasta_path = (
            Path(tmp_dir_path) / f"{samplename}_consensus_padded_sequences.fasta"
        )
        log.info(f"Writing padded consensus sequences to {fasta_path}.")
        with open(fasta_path, "w") as f:
            for consensusID, consensus in result.consensus_dicts.items():
                if consensus.consensus_padding is not None:
                    f.write(f">{consensusID}\n{consensus.consensus_padding.sequence}\n")

    if fasta_debug_path is not None:
        log.info(f"Writing padded consensus sequences to {fasta_debug_path}.")
        with open(fasta_debug_path, "w") as f:
            for consensusID, consensus in result.consensus_dicts.items():
                if consensus.consensus_padding is not None:
                    f.write(f">{consensusID}\n{consensus.consensus_padding.sequence}\n")

    log.info(f"Writing result to {output}.")
    with open(output, "w") as f:
        print(json.dumps(result.unstructure()), file=f)
    log.info("done")


def run_consensus_script(args, **kwargs):
    if args.consensus_method == "lamassemble" and args.lamassemble_mat is None:
        raise ValueError(
            "--lamassemble-mat is required when --consensus-method is 'lamassemble'."
        )
    crs_containers_to_consensus(
        samplename=args.samplename,
        input=args.input,
        path_alignments=args.alignments,
        copy_number_tracks=args.copy_number_tracks,
        output=args.output,
        lamassemble_mat=args.lamassemble_mat,
        threads=args.threads,
        crIDs=args.crIDs,
        buffer_clipped_sequence=args.buffer_clipped_sequence,
        timeout=args.timeout,
        tmp_dir_path=args.tmp_dir_path,
        figures_dir=Path(args.figures_dir) if args.figures_dir else None,
        verbose=args.verbose,
        densities_weight=args.densities_weight,
        max_intra_distance=args.max_intra_distance,
        fasta_debug_path=args.fasta_debug_path,
        cn_override=args.cn_override,
        consensus_method=args.consensus_method,
    )


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================


def get_consensus_parser(
    description="Creates Consensus objects from crs container objects.",
):
    """Create and configure the argument parser for the consensus script."""
    parser = argparse.ArgumentParser(
        description="Reads crs container objects and create Consensus objects that are written to the output database."
    )
    parser.add_argument(
        "-s", "--samplename", type=str, required=True, help="Name of the sample."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the input database that contains crs container objects.",
    )
    parser.add_argument(
        "-a",
        "--alignments",
        type=Path,
        required=True,
        help="Path to the alignments file.",
    )
    parser.add_argument(
        "-cn",
        "--copy-number-tracks",
        type=Path,
        required=True,
        help="Path to the copy number track file in bgzipped and indexed bed format.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output file that is an unstructured CrsContainerResult object.",
    )
    parser.add_argument(
        "--lamassemble-mat",
        type=Path,
        required=False,
        default=None,
        help="Path to the lamassemble mat file. Required when --consensus-method is 'lamassemble'.",
    )
    parser.add_argument(
        "--consensus-method",
        type=str,
        choices=["lamassemble", "racon"],
        default="lamassemble",
        help="Method used for consensus assembly. 'lamassemble' or 'racon' (default).",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads to use."
    )
    parser.add_argument(
        "-c",
        "--crIDs",
        type=int,
        required=False,
        default=None,
        nargs="+",
        help="List of crIDs to process. If not given, all crIDs are processed.",
    )
    parser.add_argument(
        "--cn-override",
        type=int,
        required=False,
        default=None,
        help="Instead of using the copy number from the estimation track, this fixed value is used.",
    )
    parser.add_argument(
        "--buffer-clipped-sequence",
        type=int,
        required=False,
        default=500,
        help="If parts of reads are clipped, then a maximum of this number of bases are kept to cut the read. This should never be larger than --min-cr-size in the 'svirlpool run' command.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=240,
        help="Timeout for all vs all alignment and consensus generation. Default is 240 (seconds)",
    )
    parser.add_argument(
        "--tmp-dir-path",
        required=False,
        default=None,
        help="Path to a temporary directory.",
    )
    parser.add_argument(
        "--logfile", default=None, required=False, help="Path to the logfile."
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        required=False,
        default=None,
        help="Directory to write diagnostic figures (AVA alignments, importance densities, size similarity, etc.) per candidate-region cluster. Created if it does not exist.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="prints the alignments in each clustering iteration to the terminal.",
    )
    parser.add_argument(
        "--densities-weight",
        type=float,
        default=0.0,
        help="Scaling factor for importance density weights in the pairwise similarity matrix. "
        "0.0: densities have no effect (plain signal counts); "
        ">0.0: density influence amplified. Default is 0.0.",
    )
    parser.add_argument(
        "--max-intra-distance",
        type=float,
        default=-1.0,
        help="Hard upper bound on how much intra-cluster distance is tolerated. "
        "When set to a value > 0, any pair of reads whose signal-only distance exceeds "
        "this threshold will not be in the same cluster. "
        "-1.0 disables the threshold. Default is -1.0.",
    )
    parser.add_argument(
        "--fasta-debug-path",
        type=Path,
        required=False,
        default=None,
        help="If provided, write the final padded consensus sequences to a FASTA file at this path.",
    )
    return parser


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    """Main entry point for the consensus script."""
    parser = get_consensus_parser()
    args = parser.parse_args()

    # Convert string log level to logging constant
    log_level = getattr(logging, args.log_level.upper())

    # Configure logging formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Get the root logger to ensure all loggers in the hierarchy are configured
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Remove any existing handlers to avoid duplicates
    root_logger.handlers.clear()

    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # Add file handler if logfile is specified
    if args.logfile:
        file_handler = logging.FileHandler(str(args.logfile), mode="w")
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    log.info(f"Starting consensus generation with log level {args.log_level}")
    run_consensus_script(args)
    log.info("Consensus generation completed")
    return


if __name__ == "__main__":
    main()
# %%
