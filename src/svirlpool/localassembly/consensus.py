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
import io
import json
import logging
import os
import shlex
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import cattrs
import matplotlib
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree
from sklearn.cluster import (
    AgglomerativeClustering,
    KMeans,
    SpectralClustering,
)
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

from ..signalprocessing import alignments_to_rafs, copynumber_tracks
from ..util import datatypes, util
from ..util.signal_loss_logger import get_signal_loss_logger
from . import consensus_class, consensus_lib
from . import read_cache as read_cache_mod

matplotlib.use("Agg")  # Use non-interactive backend
logging.getLogger("matplotlib").setLevel(logging.WARNING)
log = logging.getLogger(__name__)


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
    if dict_supplementary_positions:
        total = len(dict_supplementary_positions)
        missing = [
            rn for rn in dict_supplementary_positions if rn not in dict_read_sequences
        ]
        found = total - len(missing)

        if total <= 3 and found == 0:
            log.warning(
                f"No read sequences could be acquired for {total} hard-clipped read(s): "
                f"{', '.join(missing)}. Terminating assembly for this region."
            )
            return {}

        if found / total < 0.8:
            raise ValueError(
                f"Only {found}/{total} ({100 * found / total:.0f}%) read sequences could be acquired from "
                f"supplementary alignments. Missing reads: {', '.join(missing)}. "
                f"Make sure that at least one alignment of each read brings the full read sequence with it."
            )

        if missing:
            log.warning(
                f"Not all read sequences could be acquired: {found}/{total} found. "
                f"Missing read(s): {', '.join(missing)}."
            )
    # for all other alignments, just add the sequence and reverse flag.
    # Only accept alignments whose query_sequence length matches the inferred
    # read length — otherwise the stored sequence would be a hard-clipped
    # fragment rather than the full read DNA.
    expected_readnames: set[str] = {
        aln.query_name for crID in crIDs for aln in dict_alignments[crID]
    }
    for crID in crIDs:
        for aln in dict_alignments[crID]:
            if aln.query_name in dict_read_sequences:
                continue
            if not aln.query_sequence:
                continue
            if aln.infer_read_length() != len(aln.query_sequence):
                continue
            seq = (
                Seq(aln.query_sequence).reverse_complement()
                if aln.is_reverse
                else Seq(aln.query_sequence)
            )
            qualities = (
                {"phred_quality": aln.query_qualities} if aln.query_qualities else None
            )
            dict_read_sequences[aln.query_name] = SeqRecord(
                seq=seq,
                letter_annotations=qualities,
                id=aln.query_name,
                name=aln.query_name,
            )

    unresolved = [rn for rn in expected_readnames if rn not in dict_read_sequences]
    if unresolved:
        log.warning(
            f"Could not acquire full-length read sequence for {len(unresolved)}/"
            f"{len(expected_readnames)} read(s): {', '.join(unresolved)}. "
            "These reads will be excluded from the trimmed read set."
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
                if aln.query_name in cut_reads:
                    continue
                if aln.query_name not in read_records:
                    # full-length read sequence could not be acquired; skip
                    continue
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
# SIGNAL PARSING
# =============================================================================


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
    min_padding_size:int,
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
        min_padding_size=min_padding_size,
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
    min_padding_size: int,
) -> consensus_class.ConsensusPadding:
    """Creates a new ConsensusPadding object with the padding sequence and the read names."""

    used_padding_size :int = max(len(cons.consensus_sequence) * 2, min_padding_size)

    left_padding: Seq = read_records[padding_reads[0]].seq[
        padding_intervals[0][0] : padding_intervals[0][1]
    ]
    # needs to be reverse complimented if the cutread alignment is reverse
    if not read_paddings_for_consensus[padding_reads[0]][4]:
        left_padding = left_padding.reverse_complement()
    # cap padding length; keep the portion adjacent to the consensus
    left_padding = left_padding[-used_padding_size:]

    right_padding: Seq = read_records[padding_reads[1]].seq[
        padding_intervals[1][0] : padding_intervals[1][1]
    ]
    # needs to be reverse complimented if the cutread alignment is reverse
    if not read_paddings_for_consensus[padding_reads[1]][4]:
        right_padding = right_padding.reverse_complement()
    # cap padding length; keep the portion adjacent to the consensus
    right_padding = right_padding[:used_padding_size]
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


# =============================================================================
# K-MER / PCA CLUSTERING AND CONSENSUS
# =============================================================================
#
# Reads are clustered by turning each trimmed read and each racon-polished read
# into a canonical k-mer count vector, reducing the combined matrix with PCA and
# clustering the points into a fixed number of groups (one per expected allele).
# For each cluster the polished sequence closest to the cluster centroid is used
# as the consensus, and the trimmed reads whose cluster label matches are aligned
# back to that consensus by ``final_consensus``.
#
# minimap2 and racon are external binaries that require seekable files; those are
# staged on a RAM-backed tmpfs (``/dev/shm``) when available. racon writes its
# polished FASTA to stdout, which is captured and parsed in memory.

CLUSTERING_ALGORITHMS = ("kmeans", "agglomerative", "spectral", "gmm")


def _select_ramdisk(preferred: str = "/dev/shm") -> str | None:
    """Return a RAM-backed tmpfs directory if usable, else ``None``."""
    p = Path(preferred)
    if p.is_dir() and os.access(p, os.W_OK):
        return str(p)
    return None


def _dna_to_kmer_counts(
    dna_sequence: str, k: int, letter_dict: dict
) -> dict[int, int]:
    """Return canonical k-mer counts of ``dna_sequence`` as ``{kmer_hash: count}``."""
    return dict(
        util.kmer_counter_from_string(
            dna_sequence.upper(), letter_dict=letter_dict, k=k
        )
    )


def kmer_count_matrix(sequences: list[str], k: int) -> tuple[np.ndarray, list[int]]:
    """Build a dense canonical k-mer count matrix for a list of sequences.

    Returns ``(matrix, feature_hashes)`` where ``matrix`` has shape
    ``(n_sequences, n_features)``.
    """
    letter_dict = util.compute_letter_dict("ACGT")
    per_seq_counts = [
        _dna_to_kmer_counts(seq, k=k, letter_dict=letter_dict) for seq in sequences
    ]
    all_hashes: set[int] = set()
    for counts in per_seq_counts:
        all_hashes.update(counts.keys())
    feature_hashes = sorted(all_hashes)
    hash_to_col = {h: i for i, h in enumerate(feature_hashes)}
    matrix = np.zeros((len(sequences), len(feature_hashes)), dtype=np.float32)
    for row, counts in enumerate(per_seq_counts):
        for h, c in counts.items():
            matrix[row, hash_to_col[h]] = c
    return matrix, feature_hashes


def run_minimap_ava(reads_fasta: Path, paf_out: Path, threads: int = 1) -> Path:
    """Run ``minimap2 -H -U 25,35 -x ava-ont`` (all-vs-all) and write PAF overlaps."""
    cmd = shlex.split(
        f"minimap2 -H -U 25,35 -x ava-ont -t {threads} {reads_fasta} {reads_fasta}"
    )
    log.info("running minimap2 (ava): %s", " ".join(cmd))
    with open(paf_out, "w") as f:
        subprocess.check_call(cmd, stdout=f)
    return paf_out


def run_racon(
    reads_fasta: Path, paf: Path, threads: int = 1
) -> tuple[list[str], list[str]]:
    """Polish reads with racon, returning ``(names, sequences)`` parsed from stdout."""
    cmd = shlex.split(f"racon -t {threads} {reads_fasta} {paf} {reads_fasta}")
    log.info("running racon: %s", " ".join(cmd))
    completed = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    names: list[str] = []
    sequences: list[str] = []
    for record in SeqIO.parse(io.StringIO(completed.stdout.decode()), "fasta"):
        names.append(record.id)
        sequences.append(str(record.seq))
    return names, sequences


def reduce_and_cluster(
    matrix: np.ndarray,
    n_clusters: int,
    n_components: int,
    algorithm: str = "kmeans",
    random_state: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Standardize, PCA-reduce and cluster the k-mer count matrix.

    Returns ``(coords, labels)`` where ``coords`` are the PCA coordinates and
    ``labels`` is the cluster label per sample.
    """
    if algorithm not in CLUSTERING_ALGORITHMS:
        raise ValueError(
            f"unknown algorithm '{algorithm}', choose from {CLUSTERING_ALGORITHMS}"
        )
    n_samples = matrix.shape[0]
    if n_samples < n_clusters:
        raise ValueError(
            f"cannot form {n_clusters} clusters from {n_samples} sequences"
        )
    scaled = StandardScaler().fit_transform(matrix)
    n_used = max(1, min(n_components, n_samples, scaled.shape[1]))
    coords = PCA(n_components=n_used, random_state=random_state).fit_transform(scaled)

    if algorithm == "kmeans":
        labels = KMeans(
            n_clusters=n_clusters, random_state=random_state, n_init=10
        ).fit_predict(coords)
    elif algorithm == "agglomerative":
        labels = AgglomerativeClustering(
            n_clusters=n_clusters, linkage="ward"
        ).fit_predict(coords)
    elif algorithm == "spectral":
        # n_neighbors must stay below the number of samples
        n_neighbors = max(2, min(10, n_samples - 1))
        labels = SpectralClustering(
            n_clusters=n_clusters,
            random_state=random_state,
            affinity="nearest_neighbors",
            n_neighbors=n_neighbors,
            assign_labels="kmeans",
        ).fit_predict(coords)
    else:  # gmm
        labels = GaussianMixture(
            n_components=n_clusters, random_state=random_state
        ).fit_predict(coords)
    return coords, labels


def cluster_reads_kmer(
    reads: dict[str, SeqRecord],
    n_clusters: int,
    k: int = 8,
    n_components: int = 2,
    algorithm: str = "kmeans",
    threads: int = 1,
    tmp_dir_path: Path | str | None = None,
) -> dict[int, tuple[str, list[str]]]:
    """Cluster trimmed reads via k-mer/PCA and pick one consensus per cluster.

    Pipeline: minimap2 all-vs-all overlaps -> racon polish (in memory) ->
    canonical k-mer count matrix of trimmed reads + polished reads -> PCA ->
    clustering. For each cluster the polished sequence closest to the cluster
    centroid is the consensus; the trimmed reads whose label matches the cluster
    are its members.

    Returns ``{cluster_label: (consensus_sequence, [member_read_names])}``.
    """
    read_names = list(reads.keys())
    read_seqs = [str(reads[name].seq) for name in read_names]
    n_reads = len(read_names)
    if n_reads == 0:
        return {}
    # Cap the number of clusters at the number of reads.
    effective_clusters = max(1, min(n_clusters, n_reads))

    # A single read trivially forms one consensus.
    if n_reads == 1:
        return {0: (read_seqs[0], [read_names[0]])}

    ramdisk = _select_ramdisk() if tmp_dir_path is None else str(tmp_dir_path)
    with tempfile.TemporaryDirectory(prefix="kmer_clustering.", dir=ramdisk) as tmp:
        tmp_dir = Path(tmp)
        reads_ram = tmp_dir / "reads.fasta"
        with open(reads_ram, "w") as f:
            for name, seq in zip(read_names, read_seqs, strict=True):
                f.write(f">{name}\n{seq}\n")
        paf = run_minimap_ava(reads_ram, tmp_dir / "overlaps.paf", threads=threads)
        try:
            polished_names, polished_seqs = run_racon(
                reads_ram, paf, threads=threads
            )
        except subprocess.CalledProcessError as e:
            log.warning(
                "racon failed (%s); falling back to trimmed reads as consensus "
                "candidates.",
                e,
            )
            polished_names, polished_seqs = [], []

    # Combine trimmed reads + polished reads into one feature space.
    names = read_names + [f"polished:{n}" for n in polished_names]
    sequences = read_seqs + polished_seqs
    is_polished = np.array(
        [False] * len(read_seqs) + [True] * len(polished_seqs), dtype=bool
    )

    matrix, _ = kmer_count_matrix(sequences, k=k)
    coords, labels = reduce_and_cluster(
        matrix,
        n_clusters=effective_clusters,
        n_components=n_components,
        algorithm=algorithm,
    )

    read_labels = labels[:n_reads]
    result: dict[int, tuple[str, list[str]]] = {}
    for label in sorted(set(labels.tolist())):
        member_mask = labels == label
        member_reads = [
            read_names[i] for i in range(n_reads) if read_labels[i] == label
        ]
        if len(member_reads) == 0:
            continue
        # consensus = polished seq closest to cluster centroid; if the cluster
        # has no polished seq, fall back to the longest member read.
        centroid = coords[member_mask].mean(axis=0)
        polished_idx = np.where(member_mask & is_polished)[0]
        if polished_idx.size > 0:
            distances = np.linalg.norm(coords[polished_idx] - centroid, axis=1)
            best = int(polished_idx[int(np.argmin(distances))])
            consensus_sequence = sequences[best]
        else:
            consensus_sequence = max(
                (str(reads[rn].seq) for rn in member_reads), key=len
            )
        result[int(label)] = (consensus_sequence, member_reads)
    return result


def build_consensuses_kmer(
    samplename: str,
    cutreads: dict[str, SeqRecord],
    candidate_regions: dict[int, datatypes.CandidateRegion],
    n_clusters: int,
    k: int,
    n_components: int,
    algorithm: str,
    threads: int = 1,
    tmp_dir_path: Path | str | None = None,
    verbose: bool = False,
) -> dict[str, consensus_class.Consensus] | None:
    """Cluster ``cutreads`` with the k-mer/PCA method and build one Consensus per cluster."""
    if len(cutreads) == 0:
        return None
    crIDs = [cr.crID for cr in candidate_regions.values()]
    original_regions = [
        (cr.chr, cr.referenceStart, cr.referenceEnd)
        for cr in candidate_regions.values()
    ]

    try:
        clusters = cluster_reads_kmer(
            reads=cutreads,
            n_clusters=n_clusters,
            k=k,
            n_components=n_components,
            algorithm=algorithm,
            threads=threads,
            tmp_dir_path=tmp_dir_path,
        )
    except Exception as e:
        log.warning(f"k-mer clustering failed with exception: {e}. Returning None.")
        return None
    if not clusters:
        return None

    result: dict[str, consensus_class.Consensus] = {}
    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tmp_dir:
        for label, (consensus_sequence, member_reads) in sorted(clusters.items()):
            member_reads = [rn for rn in member_reads if rn in cutreads]
            if len(member_reads) == 0 or not consensus_sequence:
                continue
            consensus_name = f"{min(crIDs)}.{label}"
            consensus_fasta = tempfile.NamedTemporaryFile(
                dir=tmp_dir,
                prefix=f"consensus.kmer.{label}.",
                suffix=".fasta",
                delete=False if tmp_dir_path else True,
            )
            reads_fasta = tempfile.NamedTemporaryFile(
                dir=tmp_dir,
                prefix=f"reads.kmer.{label}.",
                suffix=".fasta",
                delete=False if tmp_dir_path else True,
            )
            with open(consensus_fasta.name, "w") as f:
                f.write(f">{consensus_name}\n{consensus_sequence}\n")
            with open(reads_fasta.name, "w") as f:
                SeqIO.write([cutreads[rn] for rn in member_reads], f, "fasta")

            consensus = final_consensus(
                samplename=samplename,
                min_indel_size=8,
                min_bnd_size=100,
                reads_fasta=Path(reads_fasta.name),
                consensus_fasta_path=Path(consensus_fasta.name),
                consensus_sequence=consensus_sequence,
                ID=consensus_name,
                crIDs=crIDs,
                original_regions=original_regions,
                threads=threads,
                verbose=verbose,
                tmp_dir_path=Path(tmp_dir),
            )
            if consensus is not None:
                result[consensus_name] = consensus
            else:
                log.warning(
                    f"final consensus generation failed for k-mer cluster {label} "
                    f"with {len(member_reads)} reads. Skipping this cluster."
                )
    return result if result else None


def process_consensus_container(
    samplename: str,
    crs_dict: dict[int, datatypes.CandidateRegion],
    read_cache: "read_cache_mod.ReadSequenceCache",
    copy_number_tracks: Path,
    timeout: int,
    buffer_clipped_length: int,
    min_padding_size: int,
    threads: int = 1,
    tmp_dir_path: Path | str | None = None,
    verbose: bool = False,
    cn_override: int | None = None,
    extra_clusters: int = 2,
    kmer_size: int = 8,
    pca_components: int = 20,
    clustering_algorithm: str = "kmeans",
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
        # Unused-reads aggregation is currently disabled downstream
        # (`crs_containers_to_consensus` does not propagate them), so we
        # return empty dicts and avoid the previously broken per-readname
        # BAM fetch that triggered for these high-CN regions.
        return {}, {}

    # Fetch alignments and full read sequences for every CR via the
    # shared sliding cache. The cache deduplicates across CRs in the
    # same batch and avoids reopening the BAM file.
    alns: dict[int, list[pysam.AlignedSegment]] = {}
    read_records: dict[str, SeqRecord] = {}
    for cr in crs_dict.values():
        cr_alns, cr_seqs = read_cache.fetch_for_cr(cr)
        alns[cr.crID] = cr_alns
        read_records.update(cr_seqs)
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
    log.info("cutting reads from alignments")
    cutreads: dict[str, SeqRecord] = trim_reads(
        dict_alignments=alns, intervals=max_intervals, read_records=read_records
    )
    log.info(f"number of reads: {len(cutreads)}")
    log.info(
        f"summed trimmed reads bp: {sum(len(read.seq) for read in cutreads.values())}"
    )

    # =========================================== CONSENSUS BUILDING =========================================== #
    # Cluster the trimmed reads with the k-mer / PCA approach and build one
    # consensus per cluster. The number of clusters defaults to the local copy
    # number plus a margin (extra_clusters).
    n_clusters = max_copy_number + extra_clusters
    log.info(
        f"k-mer clustering into up to {n_clusters} clusters "
        f"(copy number {max_copy_number} + {extra_clusters})."
    )

    consensus_objects: dict[str, consensus_class.Consensus] = {}
    res = build_consensuses_kmer(
        samplename=samplename,
        cutreads=cutreads,
        candidate_regions=crs_dict,
        n_clusters=n_clusters,
        k=kmer_size,
        n_components=pca_components,
        algorithm=clustering_algorithm,
        threads=threads,
        tmp_dir_path=tmp_dir_path,
        verbose=verbose,
    )
    if res is not None:
        consensus_objects.update(res)
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
            consensus_object=consensus,
            cutreads=cutreads,
            read_records=read_records,
            min_padding_size=min_padding_size,
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


def _load_crIDs_from_batch_tsv(batch_tsv: Path, batch_id: int) -> list[int]:
    """Return the comma-separated crID list of the requested batch row.

    The TSV is written by :mod:`svirlpool.candidateregions.crs_to_batches`
    with columns ``batch_id\tchr\tstart\tend\tcrIDs`` (header row).
    """
    if not batch_tsv.exists():
        raise FileNotFoundError(f"Batch TSV {batch_tsv} does not exist.")
    with open(batch_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            batch_id_col = header.index("batch_id")
            crIDs_col = header.index("crIDs")
        except ValueError as e:
            raise ValueError(
                f"Batch TSV {batch_tsv} is missing expected columns 'batch_id' / 'crIDs'."
            ) from e
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if int(fields[batch_id_col]) == batch_id:
                csv = fields[crIDs_col].strip()
                if not csv:
                    return []
                return [int(x) for x in csv.split(",")]
    raise ValueError(f"Batch id {batch_id} not found in {batch_tsv}.")


def crs_containers_to_consensus(
    samplename: str,
    input: Path,
    copy_number_tracks: Path,
    output: Path,
    path_alignments: Path,
    threads: int,
    buffer_clipped_sequence: int,
    timeout: int,
    reference: Path,
    crIDs: list[int] | None = None,
    tmp_dir_path: Path | str | None = None,
    verbose: bool = False,
    fasta_debug_path: Path | None = None,
    cn_override: int | None = None,
    min_padding_size: int = 30000,
    extra_clusters: int = 2,
    kmer_size: int = 8,
    pca_components: int = 20,
    clustering_algorithm: str = "kmeans",
) -> None:
    """Batch driver: process a list of containers and stream JSONL results.

    Containers are sorted by genomic coordinate within their chromosome
    and processed in order so that a single
    :class:`read_cache_mod.ReadSequenceCache` (one open BAM handle) can
    be shared across the whole batch. Each container produces one line
    of JSON (a serialised :class:`consensus_class.CrsContainerResult`)
    that is written incrementally to ``output``; nothing accumulates in
    memory.
    """
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
            # still create an empty output file so snakemake sees the target
            with open(output, "w"):
                pass
            return
    log.info(f"Loading crs containers from database at {input}.")
    containers = load_crs_containers_from_db(path_db=input, crIDs=crIDs)
    log.info("loaded data")

    # Sort containers by genomic coordinate so the sliding cache stays
    # efficient. Key by (chromosome, leftmost CR start, representative
    # crID) for deterministic order.
    def _container_sort_key(item):
        rep_crID, container = item
        crs = container["crs"]
        chrs = [cr.chr for cr in crs]
        min_start = min(cr.referenceStart for cr in crs)
        primary_chr = min(chrs)
        return (primary_chr, min_start, rep_crID)

    sorted_containers = sorted(containers.items(), key=_container_sort_key)
    n_containers = len(sorted_containers)

    _ref_fasta: pysam.FastaFile | None = (
        pysam.FastaFile(str(reference)) if reference is not None else None
    )
    cache = read_cache_mod.ReadSequenceCache(path_alignments=path_alignments)
    n_written = 0
    try:
        with open(output, "w") as out_f:
            for idx, (rep_crID, container) in enumerate(sorted_containers):
                log.info(
                    f"PROGRESS [{idx + 1}/{n_containers}] Processing container (representative crID {rep_crID})."
                )
                crs_dict = {cr.crID: cr for cr in container["crs"]}
                consensuses, _unused = process_consensus_container(
                    samplename=samplename,
                    crs_dict=crs_dict,
                    read_cache=cache,
                    tmp_dir_path=tmp_dir_path,
                    copy_number_tracks=copy_number_tracks,
                    threads=threads,
                    buffer_clipped_length=buffer_clipped_sequence,
                    timeout=timeout,
                    verbose=verbose,
                    cn_override=cn_override,
                    extra_clusters=extra_clusters,
                    kmer_size=kmer_size,
                    pca_components=pca_components,
                    clustering_algorithm=clustering_algorithm,
                    min_padding_size=min_padding_size,
                )

                # Validate consensuses immediately so the offending container
                # is identified in the log if validation fails.
                for consensusID, consensus in consensuses.items():
                    if len(consensus.original_regions) == 0:
                        raise ValueError(
                            f"Consensus {consensusID} has no original regions."
                        )

                # Stream this container's result as one JSONL record. Even
                # when no consensus could be built we emit an empty result
                # so the line count matches the container count.
                container_result = consensus_class.CrsContainerResult(
                    consensus_dicts=consensuses,
                    unused_reads={},
                )
                out_f.write(json.dumps(container_result.unstructure()) + "\n")
                out_f.flush()
                n_written += 1

                if verbose and consensuses:
                    log.info(
                        f"Container {rep_crID}: {len(consensuses)} consensuses; cache holds {cache.current_size_reads()} reads ({cache.current_size_bp()} bp)."
                    )

                if fasta_debug_path is not None:
                    _mode = "w" if idx == 0 else "a"
                    with open(fasta_debug_path, _mode) as f_debug:
                        for consensusID, consensus in consensuses.items():
                            if consensus.consensus_padding is not None:
                                f_debug.write(
                                    f">{consensusID}\n{consensus.consensus_padding.sequence}\n"
                                )
                if verbose and tmp_dir_path is not None:
                    fasta_path = (
                        Path(tmp_dir_path)
                        / f"{samplename}_consensus_padded_sequences.fasta"
                    )
                    _mode = "w" if idx == 0 else "a"
                    with open(fasta_path, _mode) as f_tmp:
                        for consensusID, consensus in consensuses.items():
                            if consensus.consensus_padding is not None:
                                f_tmp.write(
                                    f">{consensusID}\n{consensus.consensus_padding.sequence}\n"
                                )

                # Advance the sliding window so the cache can evict any
                # read whose alignment ends before the next container
                # starts on the same chromosome.
                next_window_start: int | float
                if idx + 1 < n_containers:
                    next_rep, next_container = sorted_containers[idx + 1]
                    next_chrs = {cr.chr for cr in next_container["crs"]}
                    cur_chrs = {cr.chr for cr in container["crs"]}
                    if next_chrs == cur_chrs and len(cur_chrs) == 1:
                        next_window_start = min(
                            cr.referenceStart for cr in next_container["crs"]
                        )
                    else:
                        # chromosome change is handled by the cache itself
                        # on the next fetch_for_cr call; evict everything
                        # now to bound memory.
                        next_window_start = float("inf")
                else:
                    next_window_start = float("inf")
                cache.advance(window_start=next_window_start)
    finally:
        if _ref_fasta is not None:
            _ref_fasta.close()
        cache.close()

    log.info(
        f"Wrote {n_written} container results to {output}; "
        f"cache stats: {cache.fetch_calls} fetch calls, "
        f"{cache.cache_misses_reads} new reads loaded, "
        f"{cache.cache_hits_reads} reused from cache."
    )
    log.info("done")


def run_consensus_script(args, **kwargs):
    crIDs = args.crIDs
    if getattr(args, "batch_tsv", None) is not None:
        if getattr(args, "batch_id", None) is None:
            raise ValueError("--batch-id is required when --batch-tsv is given.")
        crIDs = _load_crIDs_from_batch_tsv(Path(args.batch_tsv), int(args.batch_id))
        log.info(
            f"Loaded {len(crIDs)} container crIDs for batch {args.batch_id} from {args.batch_tsv}."
        )
    crs_containers_to_consensus(
        samplename=args.samplename,
        input=args.input,
        path_alignments=args.alignments,
        copy_number_tracks=args.copy_number_tracks,
        output=args.output,
        reference=args.reference,
        threads=args.threads,
        crIDs=crIDs,
        buffer_clipped_sequence=args.buffer_clipped_sequence,
        timeout=args.timeout,
        tmp_dir_path=args.tmp_dir_path,
        verbose=args.verbose,
        fasta_debug_path=args.fasta_debug_path,
        cn_override=args.cn_override,
        min_padding_size=args.min_padding_size,
        extra_clusters=args.extra_clusters,
        kmer_size=args.kmer_size,
        pca_components=args.pca_components,
        clustering_algorithm=args.clustering_algorithm,
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
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference genome in FASTA(.gz) format.",
    )
    parser.add_argument(
        "--extra-clusters",
        type=int,
        required=False,
        default=2,
        help="Number of clusters in addition to the local copy number "
        "(n_clusters = copy_number + extra_clusters). Default is 2.",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        required=False,
        default=8,
        help="k-mer size used to build the read/consensus feature vectors. Default is 8.",
    )
    parser.add_argument(
        "--pca-components",
        type=int,
        required=False,
        default=20,
        help="Number of leading PCA components used for clustering. Default is 20.",
    )
    parser.add_argument(
        "--clustering-algorithm",
        type=str,
        choices=list(CLUSTERING_ALGORITHMS),
        default="kmeans",
        help="Clustering algorithm used on the PCA-reduced k-mer matrix: "
        "'kmeans' (default), 'agglomerative', 'spectral' or 'gmm'.",
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
        help="List of crIDs to process. If not given, all crIDs are processed. Ignored when --batch-tsv is set.",
    )
    parser.add_argument(
        "--batch-tsv",
        type=Path,
        required=False,
        default=None,
        help="Path to a batches TSV produced by svirlpool.candidateregions.crs_to_batches. "
        "When given together with --batch-id, the crIDs to process are loaded from the matching row.",
    )
    parser.add_argument(
        "--batch-id",
        type=int,
        required=False,
        default=None,
        help="Integer id of the batch in --batch-tsv to process.",
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
        "--diag-logfile",
        default=None,
        required=False,
        help="Path to a separate diagnostic/technical logfile. Captures runtime "
        "environment (host, SLURM vars, cgroup memory limit), periodic RSS "
        "samples, received signals (SIGTERM/SIGUSR1/SIGXCPU/...), uncaught "
        "exceptions and peak memory at exit. Intended for diagnosing OOM "
        "kills, time-limit terminations and similar HPC-side failures "
        "without polluting --logfile.",
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
    parser.add_argument(
        "--fasta-debug-path",
        type=Path,
        required=False,
        default=None,
        help="If provided, write the final padded consensus sequences to a FASTA file at this path.",
    )
    parser.add_argument(
        "--min-padding-size",
        type=int,
        required=False,
        default=10000,
        help="Minimum number of bases to use for padding flanks (default: 10000).",
    )
    return parser


# =============================================================================
# DIAGNOSTIC LOGGING HELPERS
# =============================================================================
#
# Two log streams are maintained:
#   * ``log`` (this module's standard logger, propagating to root) — the
#     algorithm-oriented log that goes to --logfile.
#   * ``diag_log`` (dedicated, non-propagating) — a technical/diagnostic log
#     for environment info, signal receipt, periodic RSS samples and uncaught
#     exceptions. It writes to --diag-logfile (when given) and to stderr.
#
# Two streams keep the algorithm log readable while still preserving the
# information needed to identify OOM kills, time-limit terminations,
# illegal-instruction crashes, etc. on HPC nodes. SIGKILL leaves no Python
# traceback, so the diagnostic handlers flush after every record.


diag_log = logging.getLogger("svirlpool.consensus.diagnostics")
diag_log.propagate = False


def _read_cgroup_mem_limit_bytes():
    for path in (
        "/sys/fs/cgroup/memory.max",  # cgroup v2
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",  # cgroup v1
    ):
        try:
            with open(path) as fh:
                value = fh.read().strip()
            if value in ("", "max"):
                return None
            return int(value)
        except (OSError, ValueError):
            continue
    return None


def _current_rss_bytes():
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return None


class _FlushingFileHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        try:
            self.flush()
        except Exception:
            pass


class _FlushingStreamHandler(logging.StreamHandler):
    def emit(self, record):
        super().emit(record)
        try:
            self.flush()
        except Exception:
            pass


def _flush_all_log_handlers():
    for lg in (logging.getLogger(), diag_log):
        for h in lg.handlers:
            try:
                h.flush()
            except Exception:
                pass


def _log_environment_diagnostics():
    import socket

    diag_log.info("host=%s pid=%d", socket.gethostname(), os.getpid())
    slurm_keys = (
        "SLURM_JOB_ID",
        "SLURM_RESTART_COUNT",
        "SLURM_NODELIST",
        "SLURM_MEM_PER_NODE",
        "SLURM_MEM_PER_CPU",
        "SLURM_CPUS_PER_TASK",
        "SLURM_NTASKS",
        "SLURM_JOB_NAME",
    )
    for k in slurm_keys:
        if k in os.environ:
            diag_log.info("%s=%s", k, os.environ[k])
    cg = _read_cgroup_mem_limit_bytes()
    if cg is not None:
        diag_log.info(
            "cgroup memory limit: %d bytes (%.1f MiB)", cg, cg / (1024 * 1024)
        )
    else:
        diag_log.info("cgroup memory limit: unavailable")
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                key, _, rest = line.partition(":")
                if key.strip() in ("MemTotal", "MemAvailable"):
                    diag_log.info("/proc/meminfo %s: %s", key.strip(), rest.strip())
    except OSError:
        pass


def _install_excepthook():
    def _hook(exc_type, exc_value, exc_tb):
        diag_log.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_tb))
        _flush_all_log_handlers()

    sys.excepthook = _hook


def _install_signal_handlers():
    import signal

    # SLURM sends SIGUSR1 (when --signal=USR1@... is configured) before a job
    # is killed for time-limit; SIGTERM precedes the SIGKILL on cancellation;
    # SIGXCPU is raised when CPU time soft-limit is exceeded.
    signals = [signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGXCPU]
    for name in ("SIGUSR1", "SIGUSR2"):
        if hasattr(signal, name):
            signals.append(getattr(signal, name))

    def _handler(signum, frame):
        try:
            sig_name = signal.Signals(signum).name
        except (ValueError, AttributeError):
            sig_name = str(signum)
        rss = _current_rss_bytes()
        rss_str = f"{rss / (1024 * 1024):.1f} MiB" if rss else "unknown"
        diag_log.critical(
            "Received signal %s (%d). Current RSS: %s", sig_name, signum, rss_str
        )
        _flush_all_log_handlers()
        # Restore default disposition and re-raise so the exit status reflects
        # the signal (Snakemake/SLURM accounting then shows the real cause).
        signal.signal(signum, signal.SIG_DFL)
        os.kill(os.getpid(), signum)

    for sig in signals:
        try:
            signal.signal(sig, _handler)
        except (OSError, ValueError):
            pass


def _start_memory_monitor(interval_s: float = 10.0):
    """Periodically log RSS so we can correlate the last log line with a sudden
    SIGKILL (the OOM killer does not give the process a chance to log)."""
    import threading

    stop_event = threading.Event()

    def _run():
        while not stop_event.wait(interval_s):
            rss = _current_rss_bytes()
            if rss is not None:
                diag_log.info("memory-monitor RSS=%.1f MiB", rss / (1024 * 1024))
                _flush_all_log_handlers()

    threading.Thread(target=_run, name="memory-monitor", daemon=True).start()
    return stop_event


def _register_peak_memory_atexit():
    import atexit
    import resource as _resource

    def _log_peak():
        try:
            ru = _resource.getrusage(_resource.RUSAGE_SELF)
            # ru_maxrss is in KiB on Linux
            diag_log.info(
                "exit: peak_rss=%.1f MiB user_cpu=%.1fs sys_cpu=%.1fs",
                ru.ru_maxrss / 1024.0,
                ru.ru_utime,
                ru.ru_stime,
            )
        except Exception:
            pass
        _flush_all_log_handlers()

    atexit.register(_log_peak)


def _configure_diagnostic_logger(log_level: int, diag_logfile: str | None) -> None:
    formatter = logging.Formatter("%(asctime)s - DIAG - %(levelname)s - %(message)s")
    diag_log.setLevel(log_level)
    diag_log.handlers.clear()

    stderr_handler = _FlushingStreamHandler()
    stderr_handler.setLevel(log_level)
    stderr_handler.setFormatter(formatter)
    diag_log.addHandler(stderr_handler)

    if diag_logfile:
        diag_file_handler = _FlushingFileHandler(diag_logfile, mode="w")
        diag_file_handler.setLevel(log_level)
        diag_file_handler.setFormatter(formatter)
        diag_log.addHandler(diag_file_handler)


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    """Main entry point for the consensus script."""
    parser = get_consensus_parser()
    args = parser.parse_args()

    log_level = getattr(logging, args.log_level.upper())
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    root_logger.handlers.clear()

    console_handler = _FlushingStreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if args.logfile:
        file_handler = _FlushingFileHandler(str(args.logfile), mode="w")
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    _configure_diagnostic_logger(
        log_level=log_level,
        diag_logfile=str(args.diag_logfile) if args.diag_logfile else None,
    )

    _install_excepthook()
    _install_signal_handlers()
    _register_peak_memory_atexit()
    _log_environment_diagnostics()
    mem_monitor_stop = _start_memory_monitor(interval_s=10.0)

    log.info("Starting consensus generation with log level %s", args.log_level)
    diag_log.info("Starting consensus generation with log level %s", args.log_level)
    try:
        run_consensus_script(args)
    finally:
        mem_monitor_stop.set()
    log.info("Consensus generation completed")
    diag_log.info("Consensus generation completed")
    return


if __name__ == "__main__":
    main()
# %%
