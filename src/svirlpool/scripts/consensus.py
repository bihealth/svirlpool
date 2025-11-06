"""
Consensus Assembly Module

This script creates consensus assemblies from sets of sequencing reads.
It provides functionality for alignment processing, read clustering,
consensus generation, and scoring.
"""

# =============================================================================
# IMPORTS AND CONFIGURATION
# =============================================================================

# Standard library imports
import argparse
import hashlib
import json
import shlex
import shutil
import sqlite3
import subprocess
import tempfile
from copy import deepcopy
from pathlib import Path

# Third-party imports
import attrs
import cattrs
import matplotlib
import numpy as np
import pysam
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree
from stopit import TimeoutException

matplotlib.use("Agg")  # Use non-interactive backend
# Configure logging
import logging

import matplotlib.pyplot as plt
import networkx as nx

# Local imports
from . import (
    alignments_to_rafs,
    consensus_class,
    copynumber_tracks,
    datatypes,
    util,
)

log = logging.getLogger(__name__)


# =============================================================================
# FILE I/O AND BAM/SAM PROCESSING
# =============================================================================


def subsample_alignments(
    input_bamfile: Path, output_samfile: Path, number: int = 20, verbose: bool = False
) -> list[str]:
    """Filter BAM file to include only the longest reads up to 'number' and write to SAM format."""
    # Read all alignments from bam file
    alignments = list(pysam.AlignmentFile(input_bamfile))

    alignments_filtered = [a for a in alignments if a.query_name != a.reference_name]
    alignments_by_length = sorted(
        alignments_filtered,
        key=lambda a: a.reference_end - a.reference_start,
        reverse=True,
    )
    alignments_sampled = alignments_by_length[: min(number, len(alignments_by_length))]

    if len(alignments_by_length) == 0:
        # No valid alignments, create empty output file
        with pysam.AlignmentFile(
            output_samfile, "w", template=pysam.AlignmentFile(input_bamfile, "rb")
        ) as f:
            pass
        return []

    # Get readnames of selected alignments
    readnames = [a.query_name for a in alignments_sampled]

    # Write selected alignments to samfile
    with pysam.AlignmentFile(
        output_samfile, "w", template=pysam.AlignmentFile(input_bamfile, "rb")
    ) as f:
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
    dict_read_sequences: dict[str, SeqRecord] = (
        dict()
    )  # {(readname:(aln_reverse:bool,query_sequence:str))}
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
        for chrom in set([
            chr for l in dict_supplementary_positions.values() for chr, pos in l
        ])
    }
    for readname, l in dict_supplementary_positions.items():
        for chr, pos in l:
            dict_positions[chr].append(pos)
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
            # check if all keys of dict_read_sequences have a SeqRecord
            # if so, break this loop
            if all([
                readname in dict_read_sequences
                for readname in dict_supplementary_positions.keys()
            ]):
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
    dict_alignments: dict[int, list[pysam.AlignedSegment]] = dict()
    dict_alignments_wt: dict[int, list[pysam.AlignedSegment]] = dict()
    readnames_in_signals: dict[str] = {
        signal.readname for cr in crs for signal in cr.sv_signals
    }
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
                if aln.query_name in readnames_in_signals:
                    if cr.crID not in dict_alignments:
                        dict_alignments[cr.crID] = []
                    dict_alignments[cr.crID].append(aln)
                else:  # add reads as wt if they don't appear in the signals
                    if cr.crID not in dict_alignments_wt:
                        dict_alignments_wt[cr.crID] = []
                    dict_alignments_wt[cr.crID].append(aln)
    return dict_alignments, dict_alignments_wt


def get_read_alignment_intervals_in_region(
    region_start: int,
    regions_end: int,
    alignments: list[pysam.AlignedSegment],
    buffer_clipped_length: int,
) -> dict[str, tuple[int, int, str, int, int]]:
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


def get_read_alignment_intervals_in_cr(
    crs: list[datatypes.CandidateRegion],
    buffer_clipped_length: int,
    dict_alignments: dict[int, list[pysam.AlignedSegment]],
) -> dict[str, tuple[int, int, str, int, int]]:
    """
    Extract read alignment intervals in candidate regions.

    dict_alignments is a dict of the form {crID:[pysam.AlignedSegment]}.
    Returns a dict of the form {readname:(read_start,read_end,ref_chr,ref_start,ref_end)}
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
    dict_all_intervals: dict[str, tuple[int, int, str, int, int]] = {}
    cr_extents = {cr.crID: (cr.referenceStart, cr.referenceEnd) for cr in crs}
    # get the maximum insertion size of all original alignments in the candidate regions
    max_insertion_size = max(
        [sv.size for cr in crs for sv in cr.sv_signals if sv.sv_type == 0], default=0
    )
    used_buffer_clipped_length = max(buffer_clipped_length, max_insertion_size)
    for crID in dict_alignments.keys():
        dict_all_intervals.update(
            get_read_alignment_intervals_in_region(
                alignments=dict_alignments[crID],
                buffer_clipped_length=used_buffer_clipped_length,
                region_start=cr_extents[crID][0],
                regions_end=cr_extents[crID][1],
            )
        )
    return dict_all_intervals


def get_max_extents_of_read_alignments_on_cr(
    dict_all_intervals: dict[str, tuple[int, int, str, int, int]],
) -> dict[str, tuple[int, int, str, int, str, int]]:
    """Returns a dict of the form {readname:(read_start,read_end,ref_start,ref_end)}"""
    dict_max_extents = {}
    for readname in dict_all_intervals.keys():
        intervals = dict_all_intervals[readname]
        min_start = min([
            (start, ref_start, ref_chr)
            for (start, end, ref_chr, ref_start, ref_end) in intervals
        ])
        max_end = max([
            (end, ref_end, ref_chr)
            for (start, end, ref_chr, ref_start, ref_end) in intervals
        ])
        dict_max_extents[readname] = (
            min_start[0],
            max_end[0],
            min_start[2],
            min_start[1],
            max_end[2],
            max_end[1],
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
        for key, (
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

    cut_reads: dict[str, SeqRecord] = dict()
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
# READ CLUSTERING AND REPRESENTATIVE SELECTION
# =============================================================================


def reads_that_cover_the_region(
    readnames_in_cluster: list[str],
    max_intervals: dict[str, tuple[int, int, int, int]],
    cr_reference_start: int,
    cr_reference_end: int,
    tolerance: int = 30,
) -> list[str]:
    """Find reads that cover a specific region within tolerance."""
    # limit max_intervals to those that are in readnames
    cr_interval = Interval(cr_reference_start, cr_reference_end)
    max_intervals_in_cluster = {
        readname: Interval(min(interval[2], interval[3]), max(interval[2], interval[3]))
        for readname, interval in max_intervals.items()
        if readname in readnames_in_cluster
    }
    # for each read, check if it covers the region
    result = [
        readname
        for readname, i in max_intervals_in_cluster.items()
        if i.overlap_size(cr_interval)
        >= cr_reference_end - cr_reference_start - tolerance
    ]
    return result


def merge_signals(
    signals: list[datatypes.ExtendedSVsignal],
) -> datatypes.ExtendedSVsignal:
    """Merge multiple SV signals into a single signal."""
    # change size, start and end of the signal. pick max of strength
    merged = deepcopy(signals[0])
    for signal in signals[1:]:
        merged.size += signal.size if signal.sv_type == 0 else -signal.size
        merged.start = min(merged.start, signal.start)
        merged.end = max(merged.end, signal.end)
        merged.strength = max(merged.strength, signal.strength)
    return merged


def find_representing_read_per_cr(
    candidate_region: datatypes.CandidateRegion,
    dict_max_extents: dict[str, tuple[int, int, int, int]],
    pool: set[str],
    cr_reference_start: int,
    cr_reference_end: int,
    verbose: bool = False,
) -> tuple[str, float]:
    """
    Find the best representative read for a candidate region.

    To make a consensus with racon, a representing read is chosen. This read should
    stretch the whole candidate region and have similar extents to others.

    Args:
        dict_max_extents: dictionary of readnames and their max extents (read_start,read_end,ref_start,ref_end)
        pool: the set of readnames to choose from

    Returns:
        tuple: (readname of the chosen read, score)
    """
    if len(pool) == 0:
        raise ValueError("pool must contain at least one read")
    if len(pool) == 1:
        return (pool.pop(), 0)
    # check if all readnames of pool are in dict_max_extents
    sorted_pool = sorted(pool)
    if not all([readname in dict_max_extents for readname in pool]):
        raise ValueError(
            f"not all readnames in pool are in dict_max_extents. pool={pool}, dict_max_extents={dict_max_extents}"
        )
    # compute a minimizing ranksum matrix for all reads in the pool
    # A scores matrix has entries:
    #   1) sum of indel signals of each read in the given candidate region
    indel_sums = np.zeros(len(sorted_pool))
    for i, readname in enumerate(sorted_pool):
        indel_sums[i] = sum([
            signal.size
            for signal in candidate_region.sv_signals
            if signal.readname == readname and signal.sv_type < 3
        ])
    # calculate the difference to the mean indel sum for each read
    mean_indel_sum = np.median(indel_sums)
    indel_differences = np.abs(indel_sums - mean_indel_sum)
    #   2) difference of the read coverage to the candidate region
    distances_to_start = np.array([
        abs(dict_max_extents[readname][2] - cr_reference_start)
        for readname in sorted_pool
    ])
    distances_to_end = np.array([
        abs(dict_max_extents[readname][3] - cr_reference_end)
        for readname in sorted_pool
    ])
    distances_to_start_end = distances_to_start + distances_to_end
    # concat the two distance arrays and the read indices
    # then rank the distances_to_start_end and then the indel_differences
    # calculate the sum of both ranks
    # choose the read with the smallest sum of ranks
    ranks_start_end = np.argsort(np.argsort(distances_to_start_end))
    ranks_indel = np.argsort(np.argsort(indel_differences))
    sum_ranks = ranks_start_end + ranks_indel
    chosen = np.argmin(sum_ranks)
    if verbose:  ## print the read ins rank sum order. for each read the name, then the distance to start+end, then the indel difference, then the sum of ranks
        log.info(
            f"rank sum order for candidate region {candidate_region.crID} ({candidate_region.chr}:{candidate_region.referenceStart}-{candidate_region.referenceEnd}):"
        )
        log.info("readname\tdistance_to_start_end\tindel_difference\tsum_of_ranks")
        for i in np.argsort(sum_ranks):
            log.info(
                f"{sorted_pool[i]}\t{distances_to_start_end[i]}\t{indel_differences[i]}\t{sum_ranks[i]}"
            )

    return sorted_pool[int(chosen)], sum_ranks[chosen]


def find_representing_read(
    crs: dict[int, datatypes.CandidateRegion],
    all_intervals: dict[str, list[tuple[int, int, str, int, int]]],
    pool: set[str],
    verbose: bool = False,
) -> str:
    """Find the best representative read across multiple candidate regions. Returns the name of the representative read."""
    # construct dict of intervaltrees for all intervals
    dict_interval_trees = {cr.chr: IntervalTree() for cr in crs.values()}
    for readname, l in all_intervals.items():
        for start, end, ref_chr, ref_start, ref_end in l:
            acutal_ref_start = min(ref_start, ref_end)
            actural_ref_end = max(ref_start, ref_end)
            actual_read_start = min(start, end)
            actural_read_end = max(start, end)
            dict_interval_trees[ref_chr][acutal_ref_start:actural_ref_end] = (
                readname,
                actual_read_start,
                actural_read_end,
            )
    # for each candidate region, find the representative read.
    # if there are more than 1 different representatives, choose the one that minimizes the ranks sum
    chosen_reads: list[tuple[str, float]] = []
    for cr in crs.values():
        dict_max_extents = {}
        for interval in dict_interval_trees[cr.chr][
            cr.referenceStart : cr.referenceEnd
        ]:
            if interval.data[0] not in pool:
                continue
            dict_max_extents[interval.data[0]] = (
                interval.data[1],
                interval.data[2],
                interval.begin,
                interval.end,
            )
        if dict_max_extents == {}:
            continue
        chosen_read = find_representing_read_per_cr(
            candidate_region=cr,
            dict_max_extents=dict_max_extents,
            pool=sorted(dict_max_extents.keys()),
            cr_reference_start=cr.referenceStart,
            cr_reference_end=cr.referenceEnd,
            verbose=verbose,
        )
        chosen_reads.append(chosen_read)
    if len(chosen_reads) == 0:
        raise ValueError("no representative read could be found")
    # if there are more than 1 different representatives, choose the one that is most often. If two are equally often, choose the one with the smallest ranks sum
    dict_chosen_reads = {readname: 0 for readname, score in chosen_reads}
    for readname, score in chosen_reads:
        dict_chosen_reads[readname] += 1
    max_count = max(dict_chosen_reads.values())
    chosen_read = sorted(
        [
            (readname, score)
            for readname, score in chosen_reads
            if dict_chosen_reads[readname] == max_count
        ],
        key=lambda x: x[1],
    )[0]
    return chosen_read[0]


# =============================================================================
# CONSENSUS GENERATION WITH RACON
# =============================================================================


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
    except:
        log.warning(
            f"racon failed to produce a consensus for {name} with command\n{cmd_racon}"
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

import os


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

    # Parse alignments
    cut_read_alns: list[pysam.AlignedFragment] = list(
        pysam.AlignmentFile(tmp_alignments.name, mode="r")
    )
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
        util.display_ascii_alignments(alignments=cut_read_alns, terminal_width=125)

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
            prefix="reads.", suffix=".fastq", dir=tmp_dir, delete=False
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
        prefix="reads.", suffix=".fastq", dir=tmp_dir, delete=False
    )
    with open(tmp_reads.name, "w") as f:
        SeqIO.write([cutreads[readname] for readname in readnames], f, "fastq")
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
        from .util import display_ascii_alignments

        display_ascii_alignments(
            alignments=list(pysam.AlignmentFile(tmp_alignments_sam.name, mode="r"))
        )

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
    """Filter a fastq file by read names."""
    with open(output_fastq, "w") as out_f:
        for record in SeqIO.parse(input_fastq, "fastq"):
            if record.id in readnames:
                SeqIO.write(record, out_f, "fastq")


# This function is run after cut reads of one cluster have been aligned to their representative. This means that there are
# cut reads of one cluster written to a file (use for lamassemble or racon), and alignments in sam format to be used with
# racon, and the representative read in fasta format - to be used with racon.
# at first, lamassemble assembly is tried, then racon assembly if lamassemble fails.
def assemble_consensus(
    lamassemble_mat: Path,
    name: str,
    reads_fasta: Path,
    consensus_fasta_path: Path,
    threads: int,
    timeout: int,
    verbose: bool = False,
) -> str | None:
    # create a temporary fasta file with the read sequences filtered by chosen reads. only the chosen reads can be kept.
    with tempfile.TemporaryDirectory() as tmp_dir:
        # ====== try lamassemble first ====== #
        consensus_sequence: str | None = make_consensus_with_lamassemble(
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


from scipy.sparse import csr_matrix
from sklearn.cluster import SpectralClustering


def partition_reads_spectral(
    penalty_matrix: dict[tuple[str, str], float], all_reads: list[str], n_clusters: int
) -> dict[str, int]:
    """
    Use spectral clustering for graph partitioning.
    Works well for finding balanced clusters with minimum cut.
    """
    # Create adjacency matrix (convert penalties to similarities)
    max_penalty = max(penalty_matrix.values()) if penalty_matrix else 1.0

    # Build sparse similarity matrix
    n_reads = len(all_reads)
    read_to_idx = {read: i for i, read in enumerate(all_reads)}

    # Convert penalties to similarities: similarity = max_penalty - penalty
    similarities = []
    row_indices = []
    col_indices = []

    for (read1, read2), penalty in penalty_matrix.items():
        if read1 in read_to_idx and read2 in read_to_idx:
            similarity = max_penalty - penalty
            i, j = read_to_idx[read1], read_to_idx[read2]

            similarities.extend([similarity, similarity])  # symmetric
            row_indices.extend([i, j])
            col_indices.extend([j, i])

    # Create sparse matrix
    similarity_matrix = csr_matrix(
        (similarities, (row_indices, col_indices)), shape=(n_reads, n_reads)
    )

    # Apply spectral clustering
    clustering = SpectralClustering(
        n_clusters=n_clusters,
        affinity="precomputed",
        assign_labels="kmeans",  # or 'discretize'
        random_state=42,
    )

    labels = clustering.fit_predict(similarity_matrix)
    return {read: int(label) for read, label in zip(all_reads, labels)}


def handle_isolated_reads(
    penalty_matrix: dict[tuple[str, str], float],
    all_reads: list[str],
    connectivity_threshold: float = 0.1,
) -> tuple[list[str], list[str]]:
    """
    Separate well-connected reads from isolated ones.
    """
    # Count connections per read
    read_connections = dict.fromkeys(all_reads, 0)

    for (read1, read2), penalty in penalty_matrix.items():
        # Consider connected if penalty is below threshold
        if penalty < connectivity_threshold:
            read_connections[read1] += 1
            read_connections[read2] += 1

    # Separate based on connectivity
    min_connections = 1  # Adjust based on your needs
    well_connected = [
        read for read, count in read_connections.items() if count >= min_connections
    ]
    isolated = [
        read for read, count in read_connections.items() if count < min_connections
    ]

    return well_connected, isolated


def consensus_while_clustering(
    samplename: str,
    lamassemble_mat: Path | str,
    pool: dict[str, SeqRecord],
    candidate_regions: dict[int, datatypes.CandidateRegion],
    partitions: int,
    timeout: int,
    threads: int = 1,
    tmp_dir_path: Path | None = None,
    verbose: bool = False,
) -> dict[str, consensus_class.Consensus] | None:
    # all-vs-all alignments with minimap2
    # calculate penalties (called score_ras_from_alignments, but its really not a score, but a penalty)
    # construct a graph of reads with edges between reads given by the alignment penalties. If no penalty is found, there is no edge.
    # find a partition of the graph given N clusters that minimizes the sum of edge weights between clusters and contains all reads (min coverage set)
    crIDs = [cr.crID for cr in candidate_regions.values()]
    result: dict[str, consensus_class.Consensus] | None = None
    try:
        with tempfile.TemporaryDirectory(
            dir=tmp_dir_path, delete=False if tmp_dir_path else True
        ) as tmp_dir:
            # 1) write all read sequences to a temporary fasta file
            tmp_all_reads = tempfile.NamedTemporaryFile(
                dir=tmp_dir,
                prefix="all_reads.",
                suffix=".fastq",
                delete=False if tmp_dir_path else True,
            )
            with open(tmp_all_reads.name, "w") as f:
                SeqIO.write(pool.values(), f, "fastq")
            # 2) align all reads to each other with minimap2
            tmp_all_vs_all_sam = tempfile.NamedTemporaryFile(
                dir=tmp_dir,
                prefix="all_vs_all.",
                suffix=".sam",
                delete=False if tmp_dir_path else True,
            )
            try:
                util.align_reads_with_minimap(
                    timeout=timeout,
                    bamout=tmp_all_vs_all_sam.name,
                    reads=tmp_all_reads.name,
                    reference=tmp_all_reads.name,
                    aln_args=" --sam-hit-only --secondary=yes -U 25,75 -H",
                    tech="ava-ont",
                    threads=threads,
                )
            except TimeoutException:
                log.error(f"Alignment with minimap2 timed out after {timeout} seconds.")
                return None

            # 3) parse alignments
            all_vs_all_alignments = list(
                pysam.AlignmentFile(tmp_all_vs_all_sam.name, mode="r")
            )
            if len(all_vs_all_alignments) == 0:
                log.warning(
                    "no alignments found in all-vs-all alignment of reads. Returning None."
                )
                return None
            # 4) score alignments
            penalty_matrix: dict[tuple[str, str], float] = {}
            #  - find all reference names in the alignments
            #  - for each reference name in the alignments, do scoring
            all_raf_scores: dict[str, dict[str, float]] = {}
            for reference_name in set([
                aln.reference_name
                for aln in all_vs_all_alignments
                if not aln.is_unmapped
            ]):
                raf_alignments = [
                    aln
                    for aln in all_vs_all_alignments
                    if aln.reference_name == reference_name and not aln.is_unmapped
                ]
                raf_scores = score_ras_from_alignments(
                    samplename=samplename,
                    alignments=raf_alignments,
                    read_sequences=pool,
                    reference_sequence=(
                        str(pool[reference_name].seq) if reference_name in pool else ""
                    ),
                    min_signal_size=12,
                    min_bnd_size=100,
                    reflen=len(pool[reference_name].seq),
                    verbose=False,
                )
                all_raf_scores[reference_name] = raf_scores

            percentile_values = [
                float(np.percentile(list(raf_scores.values()), 100 / partitions))
                for raf_scores in all_raf_scores.values()
                if len(raf_scores) > 0
            ]
            dynamic_cutoff = float(np.median(percentile_values))
            dynamic_cutoff = max(dynamic_cutoff, 3.0)  # at least 3.0
            if verbose:
                print(
                    f"    DYNAMIC CUTOFF across all reference reads: {dynamic_cutoff:.2f}"
                )

            for reference_name, raf_scores in all_raf_scores.items():
                if verbose:
                    print(f"+++ {reference_name} +++")
                    for read, value in raf_scores.items():
                        print(read, value, sep="\t")
                penalty_matrix.update({
                    (reference_name, read): score
                    for read, score in raf_scores.items()
                    if score <= dynamic_cutoff
                })

            if verbose:
                # print penalty matrix as a matrix with rows and columns sorted by read names
                sorted_reads = sorted(pool.keys())
                log.info("Penalty matrix:")
                log.info("\t" + "\t".join(sorted_reads))
                for read1 in sorted_reads:
                    row = [
                        (
                            f"{penalty_matrix.get((read1, read2), float('inf')):.2f}"
                            if (read1, read2) in penalty_matrix
                            else "inf"
                        )
                        for read2 in sorted_reads
                    ]
                    log.info(f"{read1}\t" + "\t".join(row))

            all_reads = list(pool.keys())

            # 1. Handle isolated reads first - they will end up in the unused reads
            well_connected, isolated = handle_isolated_reads(penalty_matrix, all_reads)
            log.debug(
                f"Well-connected reads: {len(well_connected)}, Isolated reads: {len(isolated)}"
            )

            if len(well_connected) < partitions:
                # Not enough well-connected reads for desired partitions
                # Fall back to simple clustering or adjust partition count
                effective_partitions = max(1, len(well_connected))
            else:
                effective_partitions = partitions

            # 2. Cluster well-connected reads
            if len(well_connected) > 1:
                clustering_result = partition_reads_spectral(
                    penalty_matrix=penalty_matrix,
                    all_reads=well_connected,
                    n_clusters=effective_partitions,
                )
            else:
                clustering_result = {well_connected[0]: 0} if well_connected else {}

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
                    suffix=".fastq",
                    delete=False if tmp_dir else True,
                )
                with open(reads_fasta.name, "w") as f:
                    SeqIO.write(
                        [pool[readname] for readname in chosen_reads], f, "fastq"
                    )
                consensus_name = f"{min(crIDs)}.{cluster_id}"

                consensus_sequence: str | None = assemble_consensus(
                    lamassemble_mat=lamassemble_mat,
                    name=consensus_name,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    threads=threads,
                    timeout=timeout,
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
                    suffix=".fastq",
                    delete=False if tmp_dir else True,
                )
                with open(reads_fasta.name, "w") as f:
                    SeqIO.write(
                        [pool[readname] for readname in isolated if readname in pool],
                        f,
                        "fastq",
                    )
                consensus_name = f"{min(crIDs)}.isolated"

                consensus_sequence: str | None = assemble_consensus(
                    lamassemble_mat=lamassemble_mat,
                    name=consensus_name,
                    reads_fasta=Path(reads_fasta.name),
                    consensus_fasta_path=Path(consensus_fasta.name),
                    threads=threads,
                    timeout=timeout,
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
        alns = list(pysam.AlignmentFile(tmp_alignments.name, mode="r"))
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
    N_samples = len(set([ras.read_name for ras in rass]))
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
            if 0.1 < SeqUtils.gc_fraction(Seq(dna_string)) <= 0.9:
                filtered_signals.append(sv_signal)
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
    rafs_scores_dict = dict()
    for ras, score in zip(rass, rass_scores):
        if ras.read_name in rafs_scores_dict:
            if rafs_scores_dict[ras.read_name] < score:
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
#  Visualization functions
# =======================================================================================================================================================


def visualize_raf_score_graph(
    all_raf_scores: dict[str, dict[str, float]],
    output_filepath: Path | str | None = None,
    figsize: tuple[float, float] = (5.12, 5.12),
    node_size: int = 300,
    font_size: int = 8,
    edge_width_scale: float = 2.0,
) -> None:
    """
    Visualize a weighted graph based on RAF (Read Alignment Fragment) scores.

    This function creates a graph where:
    - Nodes represent reference sequences (from all_raf_scores keys)
    - Edges represent alignment scores between references and reads
    - Edge colors: red if weight > weighted_median, black otherwise
    - Edge thickness is proportional to the weight

    Args:
        all_raf_scores: Dictionary mapping reference names to dictionaries of
                       {read_name: score} pairs. The structure is:
                       {reference_name: {read_name: alignment_score}}
        output_filepath: Path to save the figure as PNG. If None, no figure is saved.
        figsize: Figure size in inches (width, height). Default (5.12, 5.12) = 512x512 pixels at 100 DPI.
        node_size: Size of the nodes in the graph.
        font_size: Font size for node labels.
        edge_width_scale: Scaling factor for edge widths based on weights.

    Returns:
        None. Saves a PNG file if output_filepath is provided.

    Example:
        >>> all_raf_scores = {
        ...     'ref1': {'read1': 10.5, 'read2': 5.2},
        ...     'ref2': {'read1': 3.1, 'read3': 12.0}
        ... }
        >>> visualize_raf_score_graph(all_raf_scores, 'output.png')
    """

    # Input validation
    if not all_raf_scores:
        log.warning("all_raf_scores is empty. No graph to visualize.")
        return

    # Calculate weighted median across all scores
    weights = []
    values = []
    for ref_scores in all_raf_scores.values():
        if len(ref_scores) > 0:
            weights.append(len(ref_scores))
            values.append(float(np.median(list(ref_scores.values()))))

    if not weights:
        log.warning("No valid scores found in all_raf_scores. No graph to visualize.")
        return

    weighted_median = float(np.average(values, weights=weights))

    # Create directed graph
    G = nx.DiGraph()

    # Add nodes and edges
    for reference_name, read_scores in all_raf_scores.items():
        # Add reference node if not already present
        if not G.has_node(reference_name):
            G.add_node(reference_name, node_type="reference")

        # Add edges from reference to reads with scores as weights
        for read_name, score in read_scores.items():
            # Add read node if not already present
            if not G.has_node(read_name):
                G.add_node(read_name, node_type="read")

            # Add edge with weight
            G.add_edge(reference_name, read_name, weight=score)

    if G.number_of_nodes() == 0:
        log.warning("Graph has no nodes. No graph to visualize.")
        return

    # Prepare edge colors and widths based on weighted_median threshold
    edge_colors = []
    edge_widths = []

    for u, v, data in G.edges(data=True):
        weight = data["weight"]

        # Color: red if weight > weighted_median, black otherwise
        if weight > weighted_median:
            edge_colors.append("red")
        else:
            edge_colors.append("black")

        # Width proportional to weight (normalized)
        # Avoid zero or negative widths
        normalized_weight = max(
            0.5, weight / max(1.0, weighted_median) * edge_width_scale
        )
        edge_widths.append(normalized_weight)

    # Create figure with specified DPI to ensure 512x512 pixel resolution
    dpi = 100
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Use spring layout for better visualization
    try:
        pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    except:
        # Fallback to circular layout if spring layout fails
        pos = nx.circular_layout(G)

    # Separate nodes by type for different visual styles
    reference_nodes = [
        node
        for node, attr in G.nodes(data=True)
        if attr.get("node_type") == "reference"
    ]
    read_nodes = [
        node for node, attr in G.nodes(data=True) if attr.get("node_type") == "read"
    ]

    # Draw reference nodes (larger, different color)
    if reference_nodes:
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=reference_nodes,
            node_color="lightblue",
            node_size=node_size * 1.5,
            node_shape="s",  # square
            ax=ax,
            label="References",
        )

    # Draw read nodes (smaller)
    if read_nodes:
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=read_nodes,
            node_color="lightgreen",
            node_size=node_size,
            node_shape="o",  # circle
            ax=ax,
            label="Reads",
        )

    # Draw edges with colors and widths
    nx.draw_networkx_edges(
        G,
        pos,
        edge_color=edge_colors,
        width=edge_widths,
        alpha=0.6,
        arrows=True,
        arrowsize=10,
        ax=ax,
    )

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=font_size, font_weight="bold", ax=ax)

    # Add title and legend
    ax.set_title(
        f"RAF Score Graph\n(Weighted Median: {weighted_median:.2f})",
        fontsize=font_size + 2,
        fontweight="bold",
    )
    ax.legend(loc="upper right", fontsize=font_size)
    ax.axis("off")

    plt.tight_layout()

    # Save figure if filepath provided
    if output_filepath is not None:
        output_path = Path(output_filepath)
        # Create parent directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight", format="png")
        log.info(f"Graph visualization saved to {output_path}")

    plt.close(fig)


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
        padding_sizes_per_read=padding_sizes_per_read,
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


def _get_padding_sizes_per_read(
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


def _get_padding_read_names_of_consensus(
    padding_sizes_per_read: dict[str, tuple[int, int]],
) -> tuple[str, str]:
    """Returns the read names of the left and right padding reads."""
    return (
        max(padding_sizes_per_read.items(), key=lambda x: x[1][0])[0],
        max(padding_sizes_per_read.items(), key=lambda x: x[1][1])[0],
    )


def _get_read_padding_intervals(
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


def _create_padding_object(
    cons: consensus_class.Consensus,
    read_paddings_for_consensus: dict[str, tuple[int, int, int, int]],
    padding_sizes_per_read: dict[str, tuple[int, int]],
    padding_reads: tuple[str, str],
    padding_intervals: tuple[tuple[int, int], tuple[int, int]],
    read_records: dict[str, SeqRecord],
) -> consensus_class.ConsensusPadding:
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

    new_padding = consensus_class.ConsensusPadding(
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
        if not all([isinstance(crID, int) for crID in crIDs]):
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
    containers = dict()
    for crID, data in c.fetchall():
        crs_container = json.loads(data)
        crs_container = deserialize_crs_container(crs_container)
        containers[crID] = crs_container
    c.close()
    conn.close()
    return containers


def summed_indel_distribution(
    alns: dict[int, list[pysam.AlignedSegment]],
) -> dict[str, list[int]]:
    rafs: list[datatypes.ReadAlignmentFragment] = []
    for crID, alnlist in alns.items():
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
    for raf in rafs:
        insertions = [signal.size for signal in raf.SV_signals if signal.sv_type == 0]
        deletions = [signal.size for signal in raf.SV_signals if signal.sv_type == 1]
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


def can_stop_meta_iteration(
    consensus_objects: dict[str, consensus_class.Consensus] | None,
    expected_clusters: int,
    initial_pool_readnames: set[str],
) -> bool:
    if consensus_objects is None:
        return False
    n_results = (
        0
        if (consensus_objects is None or len(consensus_objects.keys()) == 0)
        else len(consensus_objects)
    )
    if n_results < expected_clusters:
        return True
    # if many reads are unused, then a re-run is necessary
    used_reads: set[str] = set()
    for co in consensus_objects.values():
        used_reads.update(co.get_used_readnames())
    if len(used_reads) < 0.9 * len(initial_pool_readnames):
        return False
    # the number of resulting objects can also be higher, but the surplus objects should
    # be outliers that have just very few reads
    # to check this, fill a list with the maximum coverage of each consenus
    # and sort it
    # count the number of consensus objects that have less than 25% of the expected reads per haplotype
    # those are outliers. There should only be as many outliers as the number of expected clusters
    outlier_count = sum(
        1
        for co in consensus_objects.values()
        if co.get_max_cutread_coverage() < 0.15 * len(used_reads)
    )
    real_count = len(consensus_objects) - outlier_count
    if real_count > expected_clusters:
        return False
    if outlier_count > expected_clusters * int(round(0.51)):
        return False
    return True


def process_consensus_container(
    samplename: str,
    crs_dict: dict[int, datatypes.CandidateRegion],
    path_alignments: Path,
    copy_number_tracks: Path,
    lamassemble_mat: Path | str,
    timeout: int,
    threads: int = 1,
    buffer_clipped_length: int = 10000,
    tmp_dir_path: Path | None = None,
    verbose: bool = False,
) -> tuple[
    dict[str, consensus_class.Consensus], dict[int, list[datatypes.SequenceObject]]
]:
    # alns: dict[crID:list[AlignedSegment]]
    # compute max_copy_number from the intervals of the candidate regions and by querying the bgzipped and tabix indexed copy number tracks
    regions_for_cn_query = [
        (cr.chr, cr.referenceStart, cr.referenceEnd) for cr in crs_dict.values()
    ]
    max_copy_number = copynumber_tracks.query_copynumber_from_regions(
        bgzip_bed=copy_number_tracks, regions=regions_for_cn_query
    )
    log.info(f"Maximum copy number for this container: {max_copy_number}")
    if max_copy_number > 4:
        # don't process this container. Too complex.
        log.warning(
            f"Maximum copy number {max_copy_number} exceeds threshold of 4. Skipping consensus building for this container."
        )
        # parse unused reads to datatypes.SequenceObject
        dict_unused_read_names: dict[int, set[str]] = {
            cr.crID: set([signal.readname for signal in cr.sv_signals])
            for cr in crs_dict.values()
        }
        dict_unused_reads: dict[int, list[datatypes.SequenceObject]] = {
            crID: [] for crID in dict_unused_read_names
        }
        for crID, readnames in dict_unused_read_names.items():
            for readname in readnames:
                # get the read sequence from the alignment file
                with pysam.AlignmentFile(path_alignments, "rb") as samfile:
                    read_seqRecord = samfile.fetch(readname=readname)
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
        return dict(), dict_unused_reads

    alns, alns_wt = get_read_alignments_for_crs(
        crs=list(crs_dict.values()), alignments=path_alignments
    )
    dict_all_intervals = get_read_alignment_intervals_in_cr(
        crs=list(crs_dict.values()),
        dict_alignments=alns,
        buffer_clipped_length=buffer_clipped_length,
    )
    max_intervals = get_max_extents_of_read_alignments_on_cr(
        dict_all_intervals=dict_all_intervals
    )
    read_records: dict[str, SeqRecord] = get_full_read_sequences_of_alignments(
        dict_alignments=alns, path_alignments=path_alignments
    )
    log.info("cutting reads from alignments")
    cutreads: dict[str, SeqRecord] = trim_reads(
        dict_alignments=alns, intervals=max_intervals, read_records=read_records
    )
    # if there is only one cutread, return an empty consensus dict and add it to the dict_unused_reads
    # print the distribution of ins-dels of all alignments
    # parse all alignments to sv signals

    # dict_summed_indels is of form readname:[sum_inss:int,sum_dels:int]
    if verbose:
        dict_summed_indels: dict[str, list[int]] = summed_indel_distribution(alns=alns)
        print_indel_distribution(dict_summed_indels)

    # =========================================== CONSENSUS BUILDING =========================================== #

    consensus_objects: dict[str, consensus_class.Consensus] = dict()

    res: dict[str, consensus_class.Consensus] | None = None
    res = consensus_while_clustering(
        samplename=samplename,
        lamassemble_mat=lamassemble_mat,
        pool=cutreads,
        candidate_regions=crs_dict,
        partitions=max_copy_number,
        timeout=timeout,
        threads=threads,
        tmp_dir_path=tmp_dir_path,
        verbose=verbose,
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
            cr.crID: set([signal.readname for signal in cr.sv_signals])
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
        return dict(), dict_unused_reads
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

    unused_seqobjects: dict[str, datatypes.SequenceObject] = dict()
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
    dict_unused_reads: dict[int, list[datatypes.SequenceObject]] = dict()
    for crID, cr in crs_dict.items():
        # save all reads to dict_unused_reads that are in the cr and in set_unused_readnames
        readnames_in_cr = set([signal.readname for signal in cr.sv_signals])
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
    lamassemble_mat: Path | str,
    path_alignments: Path,
    threads: int,
    buffer_clipped_sequence: int,
    timeout: int,
    crIDs: list[int] | None = None,
    tmp_dir_path: Path | str | None = None,
    verbose: bool = False,
) -> None:
    assert Path(lamassemble_mat).exists(), (
        f"lamassemble matrix file {lamassemble_mat} does not exist."
    )
    assert input.exists(), f"Input file {input} does not exist."
    assert input.is_file(), f"Input file {input} is not a file."
    assert output.parent.exists(), f"Output directory {output.parent} does not exist."
    assert threads > 0, "threads must be greater than 0."
    if crIDs:
        assert type(crIDs) == list, "crIDs must be a list of integers."
        assert len(crIDs) > 0, (
            "crIDs must be a list of integers with at least one element."
        )
        assert all([isinstance(crID, int) for crID in crIDs]), (
            "crIDs must be a list of integers."
        )
    log.info("load data")
    if crIDs:
        existing_crIDs = load_crIDs_from_containers_db(path_db=input)
        # log all crIDs that are missing in a warning
        if not all([crID in existing_crIDs for crID in crIDs]):
            log.warning(
                f"crIDs {set(crIDs) - set(existing_crIDs)} are not in the database. They are skipped."
            )
        # if no real crIDs are left, return with warning
        if not any([crID in existing_crIDs for crID in crIDs]):
            log.error("No crIDs to process. Exiting.")
            return
    log.info(f"Loading crs containers from database at {input}.")
    containers = load_crs_containers_from_db(path_db=input, crIDs=crIDs)
    log.info("loaded data")
    all_consensuses: dict[str, consensus_class.Consensus] = dict()
    all_unused_reads: dict[int, list[datatypes.SequenceObject]] = dict()
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
            verbose=verbose,
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

    log.info(f"Writing result to {output}.")
    with open(output, "w") as f:
        print(json.dumps(result.unstructure()), file=f)
    log.info("done")


def run_consensus_script(args, **kwargs):
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
        verbose=args.verbose,
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
        required=True,
        help="Path to the lamassemble mat file.",
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
        "--buffer-clipped-sequence",
        type=int,
        required=False,
        default=5000,
        help="If parts of reads are clipped, then a maximum of this number of bases are kept to cut the read.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=480,
        help="Timeout for consensus generation. Default is 480 (seconds)",
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
        "--verbose",
        action="store_true",
        default=False,
        help="prints the alignments in each clustering iteration to the terminal.",
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

    if args.logfile:
        # Configure logging to file
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.FileHandler(str(args.logfile)), logging.StreamHandler()],
        )
    else:
        # Configure logging to console only
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.StreamHandler()],
        )
    run_consensus_script(args)
    return


if __name__ == "__main__":
    main()
# %%
