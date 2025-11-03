import argparse
import json
import logging
import pickle
import shlex
import subprocess
import tempfile
import time
from copy import deepcopy
from pathlib import Path, PosixPath
from typing import Generator

import attrs
import cattrs
from Bio import SeqIO
from Bio.SeqRecord import Seq, SeqRecord
from intervaltree import IntervalTree
from pysam import AlignedSegment, AlignmentFile
from tqdm import tqdm

log = logging.getLogger(__name__)

from . import (
    SVpatterns,
    SVprimitives,
    consensus_align_lib,
    consensus_class,
    datatypes,
    util,
)


def parse_crs_container_results(
    path: PosixPath,
) -> Generator[consensus_class.CrsContainerResult, None, None]:
    with open(path, "r") as f:
        for line in f:
            crs_container_result = json.loads(line)
            yield cattrs.structure(
                crs_container_result, consensus_class.CrsContainerResult
            )


def load_alignments(
    path_alignments: Path, samplename: str | None = None, parse_DNA: bool = True
) -> dict[str, list[datatypes.Alignment]]:
    """produces a dict consensusID:list[Alignment]. If a consensus has no alignments, the list is empty"""
    # load all alignments and parse from pysam.alignedSegment to datatypes.Alignment
    dict_consensus_to_alignments: dict[str, list[datatypes.Alignment]] = {}
    with open(path_alignments, "rb") as f:
        for aln in AlignmentFile(f):
            if aln.is_mapped:  # no mapped consensuses -> empty list
                if aln.query_name not in dict_consensus_to_alignments:
                    dict_consensus_to_alignments[aln.query_name] = []
                if not parse_DNA:
                    aln.query_sequence = None
                dict_consensus_to_alignments[aln.query_name].append(
                    util.parse_alignment(
                        aln=aln,
                        samplename=samplename,
                        parse_qualities=False,
                        parse_sequence=False,
                    )
                )
    return dict_consensus_to_alignments


def trf_to_interval_tree(input_trf: Path) -> dict[str, IntervalTree]:
    """produces a dict chrom:IntervalTree"""
    # load the trf intervals to an intervaltree and set their index as their names
    trf_intervals = {}
    with open(input_trf, "r") as f:
        for i, line in enumerate(f):
            chrom, start, end = line.strip().split()[0:3]
            if chrom not in trf_intervals:
                trf_intervals[chrom] = IntervalTree()
            trf_intervals[chrom][int(start) : int(end)] = i
    return trf_intervals


def alignments_to_consensusAlignments(
    alignments: dict[str, list[datatypes.Alignment]],
) -> dict[str, list[consensus_class.ConsensusAlignment]]:
    """Converts a dict of alignments to a dict of ConsensusAlignment objects.
    Args:
        alignments (dict[str,list[datatypes.Alignment]]): dict of alignments
        reference_name (str): name of the reference, e.g. GRCh38 or HS1
    Returns:
        dict[str,list[consensus_class.ConsensusAlignment]]: dict of ConsensusAlignment objects
    """
    dict_results: dict[str, list[consensus_class.ConsensusAlignment]] = (
        {}
    )  # dict consensusID:list[ConsensusAlignment]
    current_uid: int = 0
    for consensusID, alns in alignments.items():
        if consensusID not in dict_results:
            dict_results[consensusID] = []
        for aln in alns:
            consensusAlignment: consensus_class.ConsensusAlignment = (
                consensus_class.ConsensusAlignment(
                    alignment=aln, reference_name=aln.reference_name, uid=current_uid
                )
            )
            current_uid += 1
            dict_results[consensusID].append(consensusAlignment)
    return dict_results


def add_trf_annotations_to_consensusAlignments_inplace(
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    trf_intervals: dict[str, IntervalTree],
) -> None:
    """Finds the trf annotations for each alignment and produces a dict consensusID:list[tuple[Alignment,(chrom,start,end)]"""
    for consensusAlns in dict_alignments.values():
        for consensusAlignment in consensusAlns:
            chrom = str(consensusAlignment.reference_name)  # type: ignore
            if chrom in trf_intervals:
                overlapping_intervals = trf_intervals[chrom][
                    consensusAlignment.core_reference_interval[
                        0
                    ] : consensusAlignment.core_reference_interval[1]
                ]
                # for interval in overlapping_intervals:
                # assert type(interval.data) == int, f"interval.data is not an int (repeatID), but {type(interval.data)}, data={interval.data}"
                consensusAlignment.trf_intervals = [
                    (it.begin, it.end, it.data) for it in overlapping_intervals
                ]  # .append((interval.begin,interval.end,interval.data))


def align_padded_consensus_sequences(
    input_consensus_container_results: Path | str,
    input_reference: Path | str,
    path_bamout: Path | str,
    path_fastaout: Path | str,
    threads: int,
    overwrite_fasta: bool = False,
    tmp_dir_path: Path | None = None,
) -> tuple[dict[str, bytes], dict[str, tuple[int, int]]]:
    # iterate all CrsContainerResults and write the padded DNA sequences to a fasta file
    core_sequences: dict[str, bytes] = (
        {}
    )  # save key = consensusID, value = sequence (compressed with pickle and simplified padding parts)
    core_intervals: dict[str, tuple[int, int]] = (
        {}
    )  # save key = consensusID, value = (start,end) of the core sequence

    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tmp_dir:

        log.debug(
            f"Writing core sequences and intervals to fasta file: {str(path_fastaout)} "
        )
        with open(path_fastaout, "w") as fasta_writer:
            for crs_container_result in parse_crs_container_results(
                input_consensus_container_results
            ):
                log.debug(
                    f"Processing CrsContainerResult with {len(crs_container_result.consensus_dicts)} consensuses..."
                )
                for consensus in crs_container_result.consensus_dicts.values():
                    core_sequences[consensus.ID] = pickle.dumps(
                        consensus.consensus_sequence
                    )
                    core_intervals[consensus.ID] = (
                        consensus.consensus_padding.consensus_interval_on_sequence_with_padding
                    )
                    seqRec = SeqRecord(
                        seq=Seq(consensus.consensus_padding.sequence),
                        id=consensus.ID,
                        name=consensus.ID,
                        description="",
                    )
                    log.debug(
                        f"Writing consensus {consensus.ID} with length {len(consensus.consensus_padding.sequence)} to fasta file"
                    )
                    SeqIO.write(seqRec, fasta_writer, "fasta")

        if Path(path_bamout).exists() and not overwrite_fasta:
            log.info(
                f"BAM file {path_bamout} already exists. Skipping alignment to reference. set overwrite_fasta to true if this is not desired."
            )
            return core_sequences, core_intervals

        aln_args = " --secondary=no "
        tmp_bam_out_unfiltered = tempfile.NamedTemporaryFile(
            dir=tmp_dir,
            delete=False if tmp_dir_path else True,
            suffix=".bam",
            prefix="consensus_to_reference_unfiltered.",
        )

        # Inline align_consensus_sequences_to_reference function
        log.debug("aligning consensus sequences to reference with minimap2")
        # check if a minimap2 index exists
        reference_path = Path(input_reference)
        if Path(f"{input_reference}.mmi").exists():
            log.debug("using existing minimap2 index")
            reference_path = Path(f"{input_reference}.mmi")
        util.align_reads_with_minimap(
            reference=reference_path,
            bamout=tmp_bam_out_unfiltered.name,
            reads=path_fastaout,
            tech="map-ont",
            threads=threads,
            aln_args=aln_args,
        )
        # filter the consensus alignments for local alignments that are not covering the core intervals
        log.debug(
            "filtering the consensus alignments for local alignments that are not covering the core intervals..."
        )
        with AlignmentFile(tmp_bam_out_unfiltered, "rb") as f:
            header = f.header
            log.debug(
                "filtering the consensus alignments for local alignments that are not covering the core intervals..."
            )
            with AlignmentFile(path_bamout, "wb", header=header) as g:
                for aln in f:
                    if aln.is_mapped:
                        # check if the alignment covers the core interval
                        traced_back_ref_start, traced_back_ref_end = (
                            util.get_interval_on_ref_in_region(
                                a=aln,
                                start=core_intervals[aln.query_name][0],
                                end=core_intervals[aln.query_name][1],
                            )
                        )
                        if traced_back_ref_start != traced_back_ref_end:
                            g.write(aln)
        cmd_index = ["samtools", "index", path_bamout]
        subprocess.check_call(cmd_index)
    return core_sequences, core_intervals


def svPrimitives_to_svPatterns(
    SVprimitives: list[SVprimitives.SVprimitive], max_del_size: int
) -> list[SVpatterns.SVpatternType]:
    """
    Converts a list of SVprimitives to a list of SVpatterns.
    This is done by parsing the SVprimitives and creating SVpatterns from them.
    """
    if len(SVprimitives) == 0:
        log.warning("No SVprimitives provided. Returning empty list of SVpatterns.")
        return []
    current_consensusID = SVprimitives[0].consensusID

    # Group SVprimitives by consensusID
    grouped_svprimitives = {}
    for svp in SVprimitives:
        if svp.consensusID not in grouped_svprimitives:
            grouped_svprimitives[svp.consensusID] = []
        grouped_svprimitives[svp.consensusID].append(svp)

    # Process each group of SVprimitives
    svpatterns = []
    for consensusID, group in grouped_svprimitives.items():
        svpatterns.extend(
            SVpatterns.parse_SVprimitives_to_SVpatterns(
                SVprimitives=group, max_del_size=max_del_size
            )
        )

    return svpatterns


# def add_depth_to_svPrimitives_genotypes_inplace(
#         svPrimitives: list[SVprimitives.SVprimitive],
#         _covtree:dict[str,IntervalTree]) -> None:
#     """Add depth genotypes to SVprimitives based on coverage trees."""
#     for svp in svPrimitives:
#         if not svp.chr in _covtree:
#             log.debug(f"Chromosome {svp.chr} not found in coverage tree. Skipping depth assignment for SVprimitive {svp.consensusID}.")
#             continue
#         intervals = _covtree[svp.chr][svp.ref_start]
#         if len(intervals) == 0:
#             log.debug(f"No coverage data found on position {svp.chr}:{svp.ref_start}. consensusID: {svp.consensusID}.")
#             continue
#         start_depth = next(iter(intervals)).data
#         if svp.sv_type != 1:  # not a deletion
#             svp.genotypeMeasurement.add_depths(start_depth)
#         else: # deletion
#             intervals = _covtree[svp.chr][svp.ref_end]
#             if len(intervals) == 0:
#                 log.debug(f"No coverage data found on position {svp.chr}:{svp.ref_end}. consensusID: {svp.consensusID}.")
#                 continue
#             end_depth = next(iter(intervals)).data
#             svp.genotypeMeasurement.add_depths(start_depth,end_depth)


def add_consensus_sequence_and_size_distortions_to_svPatterns(
    consensus_results_path: Path | str,
    svPatterns: list[SVpatterns.SVpatternType],
    distance_scale: float,
    falloff: float,
) -> list[SVpatterns.SVpatternType]:
    """Add size distortions and inserted sequences to SVpatterns based on consensus results.

    This function processes SVpatterns one by one to avoid heavy memory re-allocations.
    Each processed element is copied to a new output list and deleted from input to free memory.

    Returns:
        list[SVpatterns.SVpatternType]: New list with processed SVpatterns
    """
    # Group SVpatterns by consensusID for efficient lookup
    svPatterns_dict = {}
    for i, svp in enumerate(svPatterns):
        if svp.consensusID not in svPatterns_dict:
            svPatterns_dict[svp.consensusID] = []
        svPatterns_dict[svp.consensusID].append((i, svp))

    # Track processed indices to avoid processing the same SVpattern twice
    processed_indices = set()
    output_svPatterns = []

    for crs_container_result in parse_crs_container_results(consensus_results_path):
        for consensus in crs_container_result.consensus_dicts.values():
            consensusID = consensus.ID
            if consensusID not in svPatterns_dict:
                continue

            for idx, svp in svPatterns_dict[consensusID]:
                if idx in processed_indices:
                    continue

                # Process the SVpattern - create a copy to avoid modifying original
                processed_svp = (
                    svp  # SVpatterns are typically immutable enough for this
                )

                # Get the size distortions from the consensus
                size_distortions = SVpatterns.distortions_by_svPattern(
                    svPattern=processed_svp,
                    consensus=consensus,
                    distance_scale=distance_scale,
                    falloff=falloff,
                )
                log.debug(
                    f"Size distortions for SVpattern {processed_svp.consensusID} in consensus {consensusID}: {size_distortions}"
                )
                if size_distortions != None and len(size_distortions) > 0:
                    processed_svp.size_distortions = size_distortions
                else:
                    raise ValueError(
                        f"No size distortions found for SVpattern {processed_svp.consensusID} in consensus {consensusID}. This should not happen. Please check the input data. svPatterns: {processed_svp}\n{consensus}"
                    )

                # Set sequences based on SVpattern type
                if isinstance(processed_svp, SVpatterns.SVpatternInsertion):
                    processed_svp.set_sequence(
                        processed_svp.get_sequence_from_consensus(consensus=consensus)
                    )
                elif isinstance(processed_svp, SVpatterns.SVpatternInversion):
                    processed_svp.set_inserted_sequence(
                        processed_svp.get_sequence_from_consensus(consensus=consensus)
                    )
                elif isinstance(processed_svp, SVpatterns.SVpatternDeletion):
                    continue  # silent skip, as deletions dont have an alt sequence
                else:
                    log.warning(
                        f"SVpattern type {type(processed_svp)} is not supported. Skipping alt sequence assignment."
                    )
                # TODO: add sequence of other (complex) SVpatterns

                # Add to output and mark as processed
                output_svPatterns.append(processed_svp)
                processed_indices.add(idx)

                # Clear reference to help with garbage collection
                svPatterns[idx] = None

    # Add any unprocessed SVpatterns (those without matching consensus)
    for i, svp in enumerate(svPatterns):
        if i not in processed_indices and svp is not None:
            output_svPatterns.append(svp)

    # Clear the input list to free memory
    svPatterns.clear()

    return output_svPatterns


def add_reference_sequence_to_svPatterns(
    svPatterns: list[SVpatterns.SVpatternType], reference_sequence: str
) -> list[SVpatterns.SVpatternType]:
    """Add reference sequence to SVpatterns."""
    # strategy: 1) save all intervals of all deletions to a regions file
    # put all intervals as keys to a dict. the values are the svPatterns.
    # 2) use samtools faidx to extract the reference sequences
    # 3) iterate over the output of getfasta and assign the sequences to the svPatterns via set_sequence
    with tempfile.TemporaryDirectory(dir=None, delete=True) as tmp_dir:
        regions: dict[
            str, list[SVpatterns.SVpatternDeletion | SVpatterns.SVpatternInversion]
        ] = {}
        output: list[SVpatterns.SVpatternType] = []

        # Separate SVpatterns that need reference sequences from those that don't
        svpatterns_needing_ref_seq = []
        svpatterns_not_needing_ref_seq = []

        for svp in svPatterns:
            if isinstance(svp, SVpatterns.SVpatternDeletion) or isinstance(
                svp, SVpatterns.SVpatternInversion
            ):
                svpatterns_needing_ref_seq.append(svp)
                region_deletion = svp.get_reference_region()
                region_key = (
                    f"{region_deletion[0]}:{region_deletion[1]}-{region_deletion[2]}"
                )
                if region_key not in regions:
                    regions[region_key] = []
                regions[region_key].append(svp)
            else:
                # Copy other SVpatterns directly to output
                svpatterns_not_needing_ref_seq.append(deepcopy(svp))

        # Add all non-deletion/inversion SVpatterns to output
        output.extend(svpatterns_not_needing_ref_seq)

        # Process deletions and inversions if any exist
        if regions:
            regions_file = Path(tmp_dir) / "regions.bed"
            with open(regions_file, "w") as f:
                for region_key, svps in regions.items():
                    for svp in svps:
                        print(region_key, file=f)
            cmd_faidx = (
                f"samtools faidx --region-file {regions_file} {reference_sequence}"
            )
            tmp_output = Path(tmp_dir) / "reference_sequences.fasta"
            subprocess.check_call(shlex.split(cmd_faidx), stdout=open(tmp_output, "w"))

            # now read the fasta file and assign the sequences to the SVpatterns
            for record in SeqIO.parse(tmp_output, "fasta"):
                region_key = f"{record.id}"
                if region_key in regions:
                    for svp in regions[region_key]:
                        x = deepcopy(svp)
                        if x.get_sv_type() == "DEL":
                            x.set_sequence(str(record.seq))
                        elif x.get_sv_type() == "INV":
                            x.set_deleted_sequence(str(record.seq))
                        output.append(x)  # append the modified copy
                        svp = None  # remove reference to original to free memory
                    regions[region_key].clear()
                else:
                    log.warning(
                        f"Region {region_key} not found in SVpatterns. This should not happen. Please check the input data."
                    )

        return output


# def add_core_reference_intervals_inplace(
#         dict_alignments:dict[str,list[consensus_class.ConsensusAlignment]],
#         core_intervals:dict[str,tuple[int,int]],
#         pysam_cache:dict[int,AlignedSegment]) -> None:
#     # each consensusAlignment has a core_reference_interval:tuple[int,int]=(0,0)
#     # it is determined by tracing back the core interval on the reference via the alignment
#     for consensus_id, alignments in tqdm(dict_alignments.items()):
#         if consensus_id not in core_intervals:
#             log.warning(f"consensusID {consensus_id} not found in core_intervals. This should not happen. Please check the input data.")
#             continue
#         core_start, core_end = core_intervals[consensus_id]
#         for consensusAlignment in alignments:
#             traced_back_ref_start, traced_back_ref_end = util.get_interval_on_ref_in_region(
#                 a=pysam_cache[consensusAlignment.uid],
#                 start=core_start,
#                 end=core_end)
#             consensusAlignment.core_reference_interval = (min(traced_back_ref_start, traced_back_ref_end), max(traced_back_ref_start, traced_back_ref_end))


import multiprocessing as mp


def _process_alignment_batch_serialized(
    serialized_batch_data: dict,
) -> list[tuple[str, int, tuple[int, int]]]:
    """Process a batch of serialized alignments to compute core reference intervals.

    Args:
        serialized_batch_data: Dictionary containing serialized batch data

    Returns:
        List of tuples (consensus_id, alignment_index, core_reference_interval)
    """
    # Deserialize the batch data
    batch_dict_alignments = {}
    for consensus_id, serialized_alignments in serialized_batch_data[
        "dict_alignments"
    ].items():
        batch_dict_alignments[consensus_id] = []
        for serialized_alignment in serialized_alignments:
            consensus_alignment = cattrs.structure(
                serialized_alignment, consensus_class.ConsensusAlignment
            )
            batch_dict_alignments[consensus_id].append(consensus_alignment)

    batch_core_intervals = serialized_batch_data["core_intervals"]

    # Deserialize pysam alignments using the deserialize function
    batch_pysam_cache = {}
    for uid_str, serialized_alignment_data in serialized_batch_data[
        "pysam_cache"
    ].items():
        uid = int(uid_str)
        batch_pysam_cache[uid] = deserialize_pysam_AlignedSegment(
            serialized_alignment_data
        )

    results = []

    for consensus_id, alignments in batch_dict_alignments.items():
        if consensus_id not in batch_core_intervals:
            log.warning(
                f"consensusID {consensus_id} not found in core_intervals. Skipping."
            )
            continue

        core_start, core_end = batch_core_intervals[consensus_id]

        for alignment_idx, consensus_alignment in enumerate(alignments):
            if consensus_alignment.uid not in batch_pysam_cache:
                log.warning(
                    f"Alignment uid {consensus_alignment.uid} not found in pysam_cache. Skipping."
                )
                continue

            traced_back_ref_start, traced_back_ref_end = (
                util.get_interval_on_ref_in_region(
                    a=batch_pysam_cache[consensus_alignment.uid],
                    start=core_start,
                    end=core_end,
                )
            )

            core_ref_interval = (
                min(traced_back_ref_start, traced_back_ref_end),
                max(traced_back_ref_start, traced_back_ref_end),
            )
            results.append((consensus_id, alignment_idx, core_ref_interval))

    return results


def _create_serialized_batch(
    batch_consensus_ids: list[str],
    dict_alignments: dict,
    core_intervals: dict,
    pysam_cache: dict,
) -> dict:
    """Create a serialized batch for the given consensus IDs.

    Args:
        batch_consensus_ids: List of consensus IDs for this batch
        dict_alignments: Full dictionary of alignments
        core_intervals: Full dictionary of core intervals
        pysam_cache: Full pysam cache

    Returns:
        Serialized batch data ready for processing
    """
    # Create sub-containers for this batch
    batch_dict_alignments = {}
    batch_core_intervals = {}
    batch_pysam_cache = {}

    # Collect all UIDs needed for this batch
    required_uids = set()

    for consensus_id in batch_consensus_ids:
        if consensus_id in dict_alignments:
            batch_dict_alignments[consensus_id] = dict_alignments[consensus_id]
            for consensus_alignment in dict_alignments[consensus_id]:
                required_uids.add(consensus_alignment.uid)

        if consensus_id in core_intervals:
            batch_core_intervals[consensus_id] = core_intervals[consensus_id]

    # Create pysam cache subset for this batch using the serialize function
    for uid in required_uids:
        if uid in pysam_cache:
            batch_pysam_cache[uid] = serialize_pysam_AlignedSegment(pysam_cache[uid])

    # Serialize the batch data
    serialized_batch = {
        "dict_alignments": {
            consensus_id: [
                consensus_alignment.unstructure() for consensus_alignment in alignments
            ]
            for consensus_id, alignments in batch_dict_alignments.items()
        },
        "core_intervals": batch_core_intervals,
        "pysam_cache": {
            str(uid): serialized_alignment
            for uid, serialized_alignment in batch_pysam_cache.items()
        },
    }

    return serialized_batch


def add_core_reference_intervals_inplace_parallel(
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    core_intervals: dict[str, tuple[int, int]],
    pysam_cache: dict[int, AlignedSegment],
    num_workers: int,
    batch_size: int = 100,
) -> None:
    """Parallelized version of add_core_reference_intervals_inplace using multiprocessing.

    Args:
        dict_alignments: Dictionary mapping consensus IDs to lists of ConsensusAlignment objects
        core_intervals: Dictionary mapping consensus IDs to core interval tuples
        pysam_cache: Dictionary mapping alignment UIDs to AlignedSegment objects
        num_workers: Number of worker processes to use
        batch_size: Number of alignments to process per batch
    """
    # Get all consensus IDs and divide them into batches
    all_consensus_ids = list(dict_alignments.keys())

    if not all_consensus_ids:
        log.info("No consensus IDs to process.")
        return

    # Create batches of consensus IDs
    consensus_id_batches = []
    for i in range(0, len(all_consensus_ids), batch_size):
        batch_consensus_ids = all_consensus_ids[i : i + batch_size]
        consensus_id_batches.append(batch_consensus_ids)

    log.info(
        f"Processing {len(all_consensus_ids)} consensus IDs in {len(consensus_id_batches)} batches using {num_workers} workers..."
    )

    # Create a generator function that yields serialized batches on-demand
    def batch_generator():
        for batch_consensus_ids in consensus_id_batches:
            yield _create_serialized_batch(
                batch_consensus_ids, dict_alignments, core_intervals, pysam_cache
            )

    # Process batches in parallel using the generator
    with mp.Pool(processes=num_workers) as pool:
        batch_results = list(
            tqdm(
                pool.imap(_process_alignment_batch_serialized, batch_generator()),
                total=len(consensus_id_batches),
                desc="Processing alignment batches",
            )
        )

    # Apply results back to the original data structure
    results_applied = 0
    for batch_result in batch_results:
        for consensus_id, alignment_idx, core_ref_interval in batch_result:
            if consensus_id in dict_alignments and alignment_idx < len(
                dict_alignments[consensus_id]
            ):
                dict_alignments[consensus_id][
                    alignment_idx
                ].core_reference_interval = core_ref_interval
                results_applied += 1

    log.info(
        f"Applied {results_applied} core reference intervals to consensus alignments."
    )


def serialize_pysam_AlignedSegment(
    pysam_Aligned_segment: AlignedSegment,
) -> dict[str, dict]:
    """Serialize a pysam AlignedSegment to a dictionary."""
    if not hasattr(pysam_Aligned_segment, "header"):
        raise ValueError("pysam_Aligned_segment must have a header attribute.")
    samdict = pysam_Aligned_segment.to_dict()
    headerdict = pysam_Aligned_segment.header.to_dict()
    return {"sam": samdict, "header": headerdict}


from pysam import AlignmentHeader


def deserialize_pysam_AlignedSegment(
    serialized_data: dict[str, dict],
) -> AlignedSegment:
    """Deserialize a dictionary back to a pysam AlignedSegment."""
    header = AlignmentHeader.from_dict(serialized_data["header"])
    return AlignedSegment.from_dict(serialized_data["sam"], header)


def add_core_reference_intervals_inplace(
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    core_intervals: dict[str, tuple[int, int]],
    pysam_cache: dict[int, AlignedSegment],
    num_workers: int,
    batch_size: int = 100,
) -> None:
    """Add core reference intervals to consensus alignments in place.

    Args:
        dict_alignments: Dictionary mapping consensus IDs to lists of ConsensusAlignment objects
        core_intervals: Dictionary mapping consensus IDs to core interval tuples
        pysam_cache: Dictionary mapping alignment UIDs to AlignedSegment objects
        use_parallel: Whether to use parallel processing (default: True)
        num_workers: Number of worker processes to use for parallel processing
        batch_size: Number of alignments to process per batch
    """
    if num_workers > 1 and len(dict_alignments) > 1:
        add_core_reference_intervals_inplace_parallel(
            dict_alignments=dict_alignments,
            core_intervals=core_intervals,
            pysam_cache=pysam_cache,
            num_workers=num_workers,
            batch_size=batch_size,
        )
    else:
        # Fallback to sequential processing for small datasets or when parallel is disabled
        for consensus_id, alignments in tqdm(
            dict_alignments.items(), desc="Processing alignments sequentially"
        ):
            if consensus_id not in core_intervals:
                log.warning(
                    f"consensusID {consensus_id} not found in core_intervals. This should not happen. Please check the input data."
                )
                continue
            core_start, core_end = core_intervals[consensus_id]
            for consensusAlignment in alignments:
                traced_back_ref_start, traced_back_ref_end = (
                    util.get_interval_on_ref_in_region(
                        a=pysam_cache[consensusAlignment.uid],
                        start=core_start,
                        end=core_end,
                    )
                )
                consensusAlignment.core_reference_interval = (
                    min(traced_back_ref_start, traced_back_ref_end),
                    max(traced_back_ref_start, traced_back_ref_end),
                )


def svPatterns_from_consensus_sequences(
    samplename: str,
    input_consensus_container_results: Path | str,
    input_reference: Path | str,
    input_trf: Path | str,
    # covtrees_db:Path|str,
    output_consensus_fasta: Path | str,
    output_consensus_to_reference_alignments: Path | str,
    threads: int,
    min_signal_size: int,
    min_bnd_size: int,
    max_del_size: int,
    distance_scale: float,
    falloff: float,
    tmp_dir_path: Path | str | None = None,
    dont_merge_horizontally: bool = False,
    enable_profiling: bool = False,
    profiling_output_dir: Path | str | None = None,
) -> list[SVpatterns.SVpatternType]:
    # write all consensus sequences to a fasta file and align to the target reference
    # write all alignments to a file
    if dont_merge_horizontally:
        log.warning(
            "dont_merge_horizontally is set to True. This will not merge SVs within the same consensus alignment. This is not recommended for general use cases."
        )
    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tmp_dir:
        if output_consensus_to_reference_alignments is None:
            output_consensus_to_reference_alignments = (
                Path(tmp_dir) / "consensus_to_reference.bam"
            )

        sequences_core, intervals_core = align_padded_consensus_sequences(
            input_consensus_container_results=input_consensus_container_results,
            path_fastaout=output_consensus_fasta,
            path_bamout=output_consensus_to_reference_alignments,
            input_reference=input_reference,
            threads=threads,
            tmp_dir_path=tmp_dir,
        )

        log.debug(
            f"loading alignments from {output_consensus_to_reference_alignments}..."
        )
        # consensusID (read name of aligned consensus) : list[Alignment]
        consensus_alignments: dict[str, list[datatypes.Alignment]] = load_alignments(
            path_alignments=output_consensus_to_reference_alignments, parse_DNA=False
        )
        n_alignments = sum(
            [len(alignments) for alignments in consensus_alignments.values()]
        )
        if n_alignments == 0:
            raise ValueError(
                f"No consensus alignments found in {output_consensus_to_reference_alignments}!"
            )
        log.info(f"{n_alignments} consensus alignments generated.")

        log.info("Converting consensus alignments to ConsensusAlignment objects...")
        dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]] = (
            alignments_to_consensusAlignments(alignments=consensus_alignments)
        )
        log.info(
            f"{sum([len(alignments) for alignments in dict_alignments.values()])} consensus alignments converted to ConsensusAlignment objects."
        )

        log.info("Caching pysam alignments for tracing back core intervals...")
        pysam_cache: dict[int, AlignedSegment] = {}
        # fill the pysam cache
        for consensusID, consensusAlignments in tqdm(dict_alignments.items()):
            for consensusAlignment in consensusAlignments:
                # Cache the pysam alignment using the alignment hash
                if consensusAlignment.uid not in pysam_cache:
                    pysam_cache[consensusAlignment.uid] = (
                        consensusAlignment.alignment.to_pysam()
                    )
        log.info(
            f"{len(pysam_cache)} pysam alignments cached for tracing back core intervals."
        )

        log.info("Adding core reference intervals to consensus alignments...")
        # TODO: slow, needs optimization
        add_core_reference_intervals_inplace(
            dict_alignments=dict_alignments,
            core_intervals=intervals_core,
            pysam_cache=pysam_cache,
            num_workers=threads,
            batch_size=500,
        )

        log.debug(
            "adding tandem repeat annotations to the consensus-to-reference alignments..."
        )
        add_trf_annotations_to_consensusAlignments_inplace(
            dict_alignments=dict_alignments,
            trf_intervals=trf_to_interval_tree(input_trf),
        )
        consensus_alignments = dict()  # has no real use anymore, free memory

        unaligned_consensusIDs = set(sequences_core.keys()) - set(
            dict_alignments.keys()
        )
        if len(unaligned_consensusIDs) > 0:
            if len(unaligned_consensusIDs) == len(dict_alignments):
                raise ValueError(
                    "all consensus sequences are unaligned to the reference. Either the path to the alignments or the fasta file is corrupted or you need to allocate more memory."
                )
            log.info(
                f"{len(unaligned_consensusIDs)} consensus sequences were not aligned to the reference."
            )
            log.debug(f"they are: {unaligned_consensusIDs}")

        log.info("Parsing SV signals from consensus alignments...")
        parse_sv_signals_total_time = 0.0
        total_calls = 0

        # Setup profiling if enabled
        profile_counter = 0

        log.info(
            "Parsing SV signals from consensus alignments and adding ALT sequences to them..."
        )

        @attrs.define
        class LastConsensus:
            consensusID: str
            sequence: str
            core_interval: tuple[int, int]

        # inplace modifying consensus_aligments is a bad problem, as it forces memory re-allocations of the whole array each time.
        # a better solution is to work in batches and to copy the modified objects to a new container

        last_consensus: LastConsensus | None = None
        dict_alignments_processed: dict[
            str, list[consensus_class.ConsensusAlignment]
        ] = {}
        consensusIDs = list(dict_alignments.keys())
        for consensusID, consensusAlignments in tqdm(dict_alignments.items()):
            if last_consensus is None or last_consensus.consensusID != consensusID:
                last_consensus = LastConsensus(
                    consensusID=consensusID,
                    sequence=pickle.loads(sequences_core[consensusID]),
                    core_interval=intervals_core[consensusID],
                )

            # Create a deep copy to avoid memory re-allocations in the original dict
            consensusAlignments_copy = deepcopy(consensusAlignments)

            for consensusAlignment in consensusAlignments_copy:
                # add proto_svs (MergedSVSignal)s to all ConsensusAlignments in dict_alignments
                # contains the alt sequences
                # parse signals (contain ALT sequences already!)
                parse_start_time = time.time()

                # Use profiling for first few calls if enabled
                if enable_profiling and profile_counter < 30:
                    profile_output_path = None
                    if profiling_output_dir:
                        profile_output_path = (
                            Path(profiling_output_dir)
                            / f"profile_consensus_{consensusID}_{profile_counter}.txt"
                        )

                    merged_svs, profiling_output = (
                        consensus_align_lib.profile_parse_sv_signals_from_consensus(
                            samplename=samplename,
                            consensusAlignment=consensusAlignment,
                            consensus_sequence=last_consensus.sequence,
                            interval_core=last_consensus.core_interval,
                            profile_output_path=profile_output_path,
                        )
                    )

                    log.info(
                        f"Profiling results for consensus {consensusID} (call {profile_counter}):"
                    )
                    # Print just the top 5 functions
                    lines = profiling_output.split("\n")
                    header_found = False
                    line_count = 0
                    for line in lines:
                        if "ncalls" in line and "tottime" in line:
                            header_found = True
                            log.info(line)
                            continue
                        if header_found and line.strip() and line_count < 5:
                            log.info(line)
                            line_count += 1
                        elif header_found and line_count >= 5:
                            break

                    profile_counter += 1
                else:
                    merged_svs: list[datatypes.MergedSVSignal] = (
                        consensus_align_lib.parse_sv_signals_from_consensus(
                            samplename=samplename,
                            consensusAlignment=consensusAlignment,
                            consensus_sequence=last_consensus.sequence,
                            interval_core=last_consensus.core_interval,
                            _pysam_cache=pysam_cache,
                            min_signal_size=min_signal_size,
                            min_bnd_size=min_bnd_size,
                        )
                    )

                parse_end_time = time.time()
                parse_sv_signals_total_time += parse_end_time - parse_start_time
                total_calls += 1
                consensusAlignment.proto_svs = merged_svs
                # TODO: I see no reason how this could happen. This check seems obsolete.
                # check all merged_svs (proto_svs) if they have different chromosomes and if so, if they share the same repeats
                for i in range(len(merged_svs)):
                    for j in range(i + 1, len(merged_svs)):
                        if merged_svs[i].chr != merged_svs[j].chr:
                            raise ValueError(
                                f"SVs {i} and {j} have different chromosomes: {merged_svs[i].chr} and {merged_svs[j].chr}. This is not allowed.\nsv a: {merged_svs[i]}\nsv b: {merged_svs[j]}"
                            )

            # Store the processed copy and clear the original to free memory
            dict_alignments_processed[consensusID] = consensusAlignments_copy
            consensusAlignments.clear()

        # Clear the original dict to free memory
        dict_alignments.clear()

        unprocessed_alignments = set(dict_alignments_processed.keys()) - set(
            sequences_core.keys()
        )
        if len(unprocessed_alignments) > 0:
            log.warning(
                f"{len(unprocessed_alignments)} consensus sequences are aligned but not present in the fasta output. This should not happen. Please check the input data and the processing steps."
            )
            for consensusID in unprocessed_alignments:
                log.warning(f" - {consensusID}")

        # =====================================================================
        #  generate all SV primitives
        # iterate all consensusResults again
        log.info(
            "Generating SV primitives from consensus objects and their alignments..."
        )

        # Setup performance logging
        performance_log_path = Path(
            "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/bndtest/performance.log"
        )
        performance_log_path.parent.mkdir(parents=True, exist_ok=True)

        # Configure performance logger
        perf_logger = logging.getLogger("performance")
        perf_logger.setLevel(logging.INFO)

        # Remove any existing handlers to avoid duplicates
        for handler in perf_logger.handlers[:]:
            perf_logger.removeHandler(handler)

        # Create file handler for performance logging
        perf_handler = logging.FileHandler(performance_log_path, mode="w")
        perf_formatter = logging.Formatter(
            "%(asctime)s,%(message)s", datefmt="%Y-%m-%d %H:%M:%S.%f"
        )
        perf_handler.setFormatter(perf_formatter)
        perf_logger.addHandler(perf_handler)
        perf_logger.propagate = False

        # Write CSV header
        perf_logger.info(
            "timestamp,consensus_id,sv_primitive_idx,consensus_length,num_sv_primitives,sum_sv_primitive_sizes,sv_type,sv_size,execution_time_ms"
        )

        # TODO: continue here: fix the multiprocessed verson of svPrimitives generation
        svPrimitives: list[SVprimitives.SVprimitive] = (
            SVprimitives.generate_SVprimitives_parallel(
                crs_container_results_iter=parse_crs_container_results(
                    input_consensus_container_results
                ),
                dict_alignments=dict_alignments_processed,
                pysam_cache=pysam_cache,
                intervals_core=intervals_core,
                samplename=samplename,
                n_workers=1,
            )
        )  # TODO: fix parallelization and change this parameter to threads
        # debug: print all types of svPrimitives that occurr in svPrimitives
        svp_types = set([svp.sv_type for svp in svPrimitives])
        log.warning(f"{len(svPrimitives)} SV primitives of types found: {svp_types}")

        # svPrimitives = [svp for svp in svPrimitives if len(svp.genotypeMeasurements.supporting_reads_start) > 0]

        log.info("Adding local depth of coverage information to the SV primitives...")
        # add coverages to svPrimitives
        # covtrees_dict:dict[str, IntervalTree] = covtree.covtree(path_db=covtrees_db)
        # add_depth_to_svPrimitives_genotypes_inplace(svPrimitives=svPrimitives, _covtree=covtrees_dict)

        svPatterns = svPrimitives_to_svPatterns(
            SVprimitives=svPrimitives, max_del_size=max_del_size
        )
        svp_types = set([svp.get_sv_type() for svp in svPatterns])
        log.warning(
            f"After svPrimitives_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
        )

        svPatterns = add_consensus_sequence_and_size_distortions_to_svPatterns(
            consensus_results_path=input_consensus_container_results,
            svPatterns=svPatterns,
            distance_scale=distance_scale,
            falloff=falloff,
        )
        svp_types = set([svp.get_sv_type() for svp in svPatterns])
        log.warning(
            f"After add_consensus_sequence_and_size_distortions_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
        )

        svPatterns = add_reference_sequence_to_svPatterns(
            svPatterns=svPatterns, reference_sequence=input_reference
        )
        svp_types = set([svp.get_sv_type() for svp in svPatterns])
        log.warning(
            f"After add_reference_sequence_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
        )

        return svPatterns


def get_parser():
    parser = argparse.ArgumentParser(
        description="Align consensus sequences to reference and extract SV primitives"
    )
    parser.add_argument("--samplename", type=str, required=True, help="Sample name.")
    parser.add_argument(
        "--consensus",
        type=Path,
        required=True,
        help="Path to consensus sequences file (consensus_containers.txt).",
    )
    parser.add_argument(
        "--reference", type=Path, required=True, help="Path to reference genome file."
    )
    parser.add_argument(
        "--trf",
        type=Path,
        required=True,
        help="Path to tandem repeat finder output file.",
    )
    parser.add_argument(
        "--output-consensus-fasta",
        type=Path,
        required=True,
        help="Output path for consensus FASTA file.",
    )
    parser.add_argument(
        "--output-consensus-alignments",
        type=Path,
        required=True,
        help="Output path for consensus alignments BAM file.",
    )
    parser.add_argument(
        "--output-svpatterns",
        type=Path,
        required=True,
        help="Output path for SV patterns database file.",
    )
    parser.add_argument(
        "--min-signal-size",
        type=int,
        default=12,
        help="Minimum size of SV signal to consider (default: 12). Increase with low per-base error rates.",
    )
    parser.add_argument(
        "--min-bnd-size",
        type=int,
        default=200,
        help="Minimum size of BND signal to consider (default: 200).",
    )
    parser.add_argument(
        "--max-del-size",
        type=int,
        default=100_000,
        help="Maximum size of deletions to consider for SVpattern generation (default: 100000). They are treated as translocations otherwise. Warning: larger deletions are often spurious and can lead to memory issues.",
    )
    parser.add_argument(
        "--distance-scale",
        type=float,
        default=5000.0,
        help="Distance scale for size distortion calculation (default: 5000.0). Controls the width of the decay function. Make much smaller if no tandem repeats are present.",
    )
    parser.add_argument(
        "--falloff",
        type=float,
        default=1.0,
        help="Falloff for size distortion calculation (default: 1.0). This controls the shape of the distance decay function. A lower value means a slower decay, which is better at noisy loci.",
    )
    parser.add_argument(
        "--sqlite-timeout",
        type=int,
        default=30,
        help="SQLite timeout in seconds (default: 30).",
    )
    parser.add_argument(
        "--tmp-dir-path",
        type=Path,
        required=False,
        default=None,
        help="Path to temporary directory for intermediate files.",
    )
    parser.add_argument(
        "--dont-merge-horizontally",
        action="store_true",
        default=False,
        help="Don't merge SV signals horizontally.",
    )
    parser.add_argument(
        "--enable-profiling",
        action="store_true",
        default=False,
        help="Enable performance profiling for the first few consensus calls.",
    )
    parser.add_argument(
        "--profiling-output-dir",
        type=Path,
        required=False,
        default=None,
        help="Directory to save detailed profiling results.",
    )
    parser.add_argument(
        "--threads", type=int, default=4, help="Number of threads to use (default: 4)."
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO).",
    )
    return parser


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    """Main entry point for the consensus script."""
    parser = get_parser()
    args = parser.parse_args()

    # Configure logging level
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    svPatterns = svPatterns_from_consensus_sequences(
        samplename=args.samplename,
        input_consensus_container_results=args.consensus,
        input_reference=args.reference,
        input_trf=args.trf,
        # covtrees_db=args.coverage_db,
        output_consensus_fasta=args.output_consensus_fasta,
        output_consensus_to_reference_alignments=args.output_consensus_alignments,
        threads=args.threads,
        min_signal_size=args.min_signal_size,
        min_bnd_size=args.min_bnd_size,
        max_del_size=args.max_del_size,
        distance_scale=args.distance_scale,
        falloff=args.falloff,
        tmp_dir_path=args.tmp_dir_path,
        dont_merge_horizontally=args.dont_merge_horizontally,
        enable_profiling=args.enable_profiling,
        profiling_output_dir=args.profiling_output_dir,
    )

    SVpatterns.create_svPatterns_db(database=args.output_svpatterns)
    SVpatterns.write_svPatterns_to_db(
        database=args.output_svpatterns, data=svPatterns, timeout=args.sqlite_timeout
    )
    return


if __name__ == "__main__":
    main()
