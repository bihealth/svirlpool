import argparse
import json
import logging
import multiprocessing as mp
import pickle
import shlex
import subprocess
import tempfile
import time
from collections.abc import Generator
from copy import deepcopy
from io import TextIOWrapper
from pathlib import Path, PosixPath

import attrs
import cattrs
from Bio import SeqIO
from Bio.SeqRecord import Seq, SeqRecord
from intervaltree import IntervalTree
from pysam import AlignedSegment, AlignmentFile, AlignmentHeader
from tqdm import tqdm

from ..util import datatypes, util
from . import SVpatterns, SVprimitives, consensus_align_lib, consensus_class

log = logging.getLogger(__name__)


def parse_crs_container_results(
    path: Path | str,
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
    dict_results: dict[
        str, list[consensus_class.ConsensusAlignment]
    ] = {}  # dict consensusID:list[ConsensusAlignment]
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
    padded_sequences_paths: dict[int, Path],
    input_reference: Path | str,
    path_bamout: Path | str,
    threads: int,
    tmp_dir_path: Path | None = None,
) -> None:
    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tmp_dir:
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
    return core_intervals


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
    _current_consensusID = SVprimitives[0].consensusID  # FIXME: unused?

    # Group SVprimitives by consensusID
    grouped_svprimitives = {}
    for svp in SVprimitives:
        if svp.consensusID not in grouped_svprimitives:
            grouped_svprimitives[svp.consensusID] = []
        grouped_svprimitives[svp.consensusID].append(svp)

    # Process each group of SVprimitives
    svpatterns = []
    for _consensusID, group in grouped_svprimitives.items():
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
                if size_distortions is not None and len(size_distortions) > 0:
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
                    for _svp in svps:  # FIXME: unused?
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


def write_tmp_core_sequences_files_for_parallel_processing(
        input_consensus_container_results: Path | str,
        core_sequences_paths: dict[int, Path],
        consensus_paths: dict[int, Path],
        padded_sequences_paths: dict[int, Path]) -> dict[int, set[str]]:
    """Writes temporary core sequence fasta, padded sequence fasta, and consensus json tsv files for parallel processing."""
    # check if padded_sequences_paths and core_sequences_paths have the same keys
    if set(core_sequences_paths.keys()) != set(padded_sequences_paths.keys()):
        raise ValueError("core_sequences_paths and padded_sequences_paths must have the same keys.")
    # write the core sequences to fasta files
    # 1) open all fasta writers
    n_files = len(core_sequences_paths)
    core_writers:dict[int,TextIOWrapper] = {}
    for i in range(n_files):
        core_writers[i] = open(core_sequences_paths[i], "w")
    padded_writers:dict[int,TextIOWrapper] = {}
    for i in range(n_files):
        padded_writers[i] = open(padded_sequences_paths[i], "w")
    consensus_writers:dict[int,TextIOWrapper] = {}
    for i in range(n_files):
        consensus_writers[i] = open(consensus_paths[i], "w")
    # 2) iterate over all consensus sequences and write them in a round-robin fashion
    writer_idx:int = 0
    index_consensusIDs:dict[int, set[str]] = {i:set() for i in range(n_files)} # to find the correct fasta file for each consensusID

    for crs_container_result in parse_crs_container_results(
            input_consensus_container_results
        ):
            log.debug(
                f"Processing CrsContainerResult with {len(crs_container_result.consensus_dicts)} consensuses..."
            )
            for consensus in crs_container_result.consensus_dicts.values():
                if consensus.consensus_padding is None:
                    raise ValueError(
                        f"Consensus {consensus.ID} has no consensus_padding. This should not happen. Please check the input data."
                    )
                core_seqRec = SeqRecord(
                    seq=Seq(consensus.consensus_sequence),
                    id=consensus.ID,
                    name=consensus.ID,
                    description=""
                )
                log.debug(
                    f"Writing consensus {consensus.ID} with length {len(consensus.consensus_sequence)} to fasta file"
                )
                SeqIO.write(core_seqRec, core_writers[writer_idx], "fasta")
                
                padded_seqRec = SeqRecord(
                    seq=Seq(consensus.consensus_padding.sequence),
                    id=consensus.ID,
                    name=consensus.ID,
                    description="",
                )
                log.debug(
                    f"Writing padded consensus {consensus.ID} with length {len(consensus.consensus_padding.sequence)} to fasta file"
                )
                SeqIO.write(padded_seqRec, padded_writers[writer_idx], "fasta")
                # serialize consensus and write to tsv (consenus.ID \t serialized_consensus)
                
                consensus.consensus_sequence = ""  # remove the sequences to save space and time
                consensus.consensus_padding.sequence = ""  # remove the sequences to save space and time
                print(str(consensus.ID),json.dumps(consensus.unstructure()), sep='\t', file=consensus_writers[writer_idx])
                
                # core interval would be consensus.consensus_padding.consensus_interval_on_sequence_with_padding
                index_consensusIDs[writer_idx].add(consensus.ID)
            writer_idx = (writer_idx + 1) % n_files # round-robin
    # close all fasta writers
    for i in range(n_files):
        core_writers[i].close()
        padded_writers[i].close()
        consensus_writers[i].close()
    return index_consensusIDs


def index_tmp_fastas(paths: list[Path]) -> None:
    """Indexes temporary fasta files using samtools faidx."""
    for path in paths:
        cmd_index = ["samtools", "faidx", str(path)]
        subprocess.check_call(cmd_index)


def write_partitioned_alignments(
    input_bam: Path,
    output_paths: dict[int, Path],
    index_consensusIDs: dict[int, set[str]],
    samplename: str,
) -> None:
    """Parse BAM file and write partitioned datatypes.Alignment objects to TSV files.
    
    Args:
        input_bam: Path to input BAM file
        output_paths: Dictionary mapping partition index to output TSV path
        index_consensusIDs: Dictionary mapping partition index to set of consensusIDs
    """
    # Generate a dict consensusID:index for faster lookup
    consensusID_to_index: dict[str, int] = {
        consensusID: i 
        for i, consensusIDs in index_consensusIDs.items() 
        for consensusID in consensusIDs
    }
    
    # Open all output TSV files
    writers: dict[int, TextIOWrapper] = {}
    for i in output_paths.keys():
        writers[i] = open(output_paths[i], 'w')
    
    # Parse BAM file and write alignments to partitioned files
    with AlignmentFile(str(input_bam), "rb") as bam_in:
        for pysam_aln in bam_in:
            consensusID = pysam_aln.query_name
            
            if consensusID not in consensusID_to_index:
                raise ValueError(
                    f"consensusID {consensusID} not found in index_consensusIDs. "
                    "This should not happen. Please check the input data."
                )
            
            # Convert pysam alignment to datatypes.Alignment
            alignment = datatypes.Alignment.from_pysam(pysam_aln, samplename)
            
            # Write to appropriate partition file: consensusID\tserialized_alignment
            idx = consensusID_to_index[consensusID]
            print(
                consensusID,
                json.dumps(alignment.unstructure()),
                sep='\t',
                file=writers[idx]
            )
    
    # Close all writers
    for writer in writers.values():
        writer.close()
    

def build_alignments_index(alignment_paths: dict[int, Path]) -> dict[str, list[tuple[int, int]]]:
    """Build an index mapping consensusID to list of (file_index, byte_offset) for alignments.
    
    Args:
        alignment_paths: Dictionary mapping file indices to alignment TSV paths
        
    Returns:
        Dictionary mapping consensusID to list of (file_index, byte_offset) tuples,
        one entry per alignment of that consensusID
    """
    index: dict[str, list[tuple[int, int]]] = {}
    
    for file_idx, path in alignment_paths.items():
        with open(path, 'r') as f:
            offset = 0
            for line in f:
                consensusID = line.split('\t', 1)[0]
                if consensusID not in index:
                    index[consensusID] = []
                index[consensusID].append((file_idx, offset))
                offset = f.tell()
    
    return index


def build_alignments_index_partitioned(alignment_paths: dict[int, Path]) -> dict[int, dict[str, list[int]]]:
    """Build partitioned alignment indices directly, partitioned by file.
    
    Args:
        alignment_paths: Dictionary mapping file indices to alignment TSV paths
        
    Returns:
        Dictionary mapping partition_idx -> consensusID -> list of byte_offsets
        Each consensusID can have multiple alignments (multiple offsets) in the same file
    """
    partitioned_index: dict[int, dict[str, list[int]]] = {}
    
    for file_idx, path in alignment_paths.items():
        partitioned_index[file_idx] = {}
        with open(path, 'r') as f:
            offset = 0
            for line in f:
                consensusID = line.split('\t', 1)[0]
                if consensusID not in partitioned_index[file_idx]:
                    partitioned_index[file_idx][consensusID] = []
                partitioned_index[file_idx][consensusID].append(offset)
                offset = f.tell()
    
    return partitioned_index


def build_consensus_tsv_index(consensus_paths: dict[int, Path]) -> dict[str, tuple[int, int]]:
    """Build an index mapping consensusID to (file_index, byte_offset)."""
    index: dict[str, tuple[int, int]] = {}
    
    for file_idx, path in consensus_paths.items():
        with open(path, 'r') as f:
            offset = 0
            for line in f:
                consensusID = line.split('\t', 1)[0]
                index[consensusID] = (file_idx, offset)
                offset = f.tell()
    
    return index


def build_consensus_tsv_index_partitioned(consensus_paths: dict[int, Path]) -> dict[int, dict[str, int]]:
    """Build partitioned indices mapping consensusID to byte_offset, partitioned by file.
    
    Args:
        consensus_paths: Dictionary mapping file indices to consensus TSV paths
        
    Returns:
        Dictionary mapping partition_idx -> consensusID -> byte_offset
    """
    partitioned_index: dict[int, dict[str, int]] = {}
    
    for file_idx, path in consensus_paths.items():
        partitioned_index[file_idx] = {}
        with open(path, 'r') as f:
            offset = 0
            for line in f:
                consensusID = line.split('\t', 1)[0]
                partitioned_index[file_idx][consensusID] = offset
                offset = f.tell()
    
    return partitioned_index


def _process_alignment_file_for_core_intervals(
    file_idx: int,
    padded_alignments_path: Path,
    index_consensusIDs: set[str],
    consensus_json_index_partition: dict[str, int],
    consensus_json_path: Path,
    batch_size: int = 100
) -> dict[str, list[tuple[int, str, int, int]]]:
    """Process a single alignment file to extract core intervals on reference.
    
    This function is designed to be run in parallel for each alignment file.
    
    Args:
        file_idx: Index of the file being processed
        padded_alignments_path: Path to the serialized alignments TSV file
        index_consensusIDs: Set of consensus IDs expected in this file
        consensus_json_index_partition: Partitioned index mapping consensusID to byte_offset for this partition only
        consensus_json_path: Path to the consensus TSV file for this partition
        batch_size: Number of alignments to process per batch
        
    Returns:
        Dictionary mapping consensusID to list of (alignment_idx, ref_name, core_start_on_ref, core_end_on_ref)
        The alignment_idx corresponds to the order in the alignment file for this consensusID
    """
    core_intervals_on_reference: dict[str, list[tuple[int, str, int, int]]] = {}
    alignment_counters: dict[str, int] = {}  # Track alignment index per consensusID
    
    log.debug(f"Worker {file_idx}: Processing {padded_alignments_path}")
    
    alignment_batch = []
    consensusID_batch = []
    alignment_idx_batch = []  # Track indices
    
    # Read serialized alignments from TSV file
    with open(padded_alignments_path, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t', 1)
            if len(fields) < 2:
                continue
            
            consensusID = fields[0]
            if consensusID not in index_consensusIDs:
                raise ValueError(
                    f"consensusID {consensusID} found in alignment file but not in "
                    f"index_consensusIDs for index {file_idx}. This should not happen."
                )
            
            # Track alignment index for this consensusID
            if consensusID not in alignment_counters:
                alignment_counters[consensusID] = 0
            current_alignment_idx = alignment_counters[consensusID]
            alignment_counters[consensusID] += 1
            
            # Deserialize alignment
            alignment = cattrs.structure(
                json.loads(fields[1]),
                datatypes.Alignment
            )
            
            # Add to batch
            alignment_batch.append(alignment)
            consensusID_batch.append(consensusID)
            alignment_idx_batch.append(current_alignment_idx)
            
            # Process batch when it reaches batch_size
            if len(alignment_batch) >= batch_size:
                # Get all consensus objects for this batch at once
                # Reconstruct the index format expected by get_consensus_batch_by_ids
                batch_index = {
                    cid: (file_idx, consensus_json_index_partition[cid])
                    for cid in consensusID_batch
                    if cid in consensus_json_index_partition
                }
                consensus_objs = get_consensus_batch_by_ids(
                    consensusIDs=consensusID_batch,
                    index=batch_index,
                    consensus_paths={file_idx: consensus_json_path}
                )
                
                # Process each alignment in the batch
                for alignment_item, consensusID_item, aln_idx in zip(
                    alignment_batch, consensusID_batch, alignment_idx_batch
                ):
                    consensus_obj = consensus_objs[consensusID_item]
                    
                    # Convert to pysam for processing
                    pysam_aln = alignment_item.to_pysam()
                    
                    # Get core intervals on padded sequence
                    core_interval_on_reference: tuple[str, int, int] = (
                        consensus_class.get_consensus_core_alignment_interval_on_reference(
                            consensus=consensus_obj,
                            alignment=pysam_aln
                        )
                    )
                    
                    # Get core interval on reference via alignment
                    traced_back_ref_start, traced_back_ref_end = util.get_interval_on_ref_in_region(
                        a=pysam_aln,
                        start=core_interval_on_reference[1],
                        end=core_interval_on_reference[2]
                    )
                    
                    if consensusID_item not in core_intervals_on_reference:
                        core_intervals_on_reference[consensusID_item] = []
                    
                    # Store with alignment index for matching
                    core_intervals_on_reference[consensusID_item].append(
                        (
                            aln_idx,  # Alignment index within this consensusID
                            str(pysam_aln.reference_name),
                            min(traced_back_ref_start, traced_back_ref_end),
                            max(traced_back_ref_start, traced_back_ref_end)
                        )
                    )
                
                # Clear the batch
                alignment_batch = []
                consensusID_batch = []
        
        # Process remaining alignments in the last batch
        if alignment_batch:
            # Reconstruct the index format expected by get_consensus_batch_by_ids
            batch_index = {
                cid: (file_idx, consensus_json_index_partition[cid])
                for cid in consensusID_batch
                if cid in consensus_json_index_partition
            }
            consensus_objs = get_consensus_batch_by_ids(
                consensusIDs=consensusID_batch,
                index=batch_index,
                consensus_paths={file_idx: consensus_json_path}
            )
            
            for alignment_item, consensusID_item, aln_idx in zip(
                alignment_batch, consensusID_batch, alignment_idx_batch
            ):
                consensus_obj = consensus_objs[consensusID_item]
                pysam_aln = alignment_item.to_pysam()
                
                core_interval_on_reference = (
                    consensus_class.get_consensus_core_alignment_interval_on_reference(
                        consensus=consensus_obj,
                        alignment=pysam_aln
                    )
                )
                
                traced_back_ref_start, traced_back_ref_end = util.get_interval_on_ref_in_region(
                    a=pysam_aln,
                    start=core_interval_on_reference[1],
                    end=core_interval_on_reference[2]
                )
                
                if consensusID_item not in core_intervals_on_reference:
                    core_intervals_on_reference[consensusID_item] = []
                
                core_intervals_on_reference[consensusID_item].append(
                    (
                        aln_idx,
                        str(pysam_aln.reference_name),
                        min(traced_back_ref_start, traced_back_ref_end),
                        max(traced_back_ref_start, traced_back_ref_end)
                    )
                )
    
    log.debug(f"Worker {file_idx}: Processed {len(core_intervals_on_reference)} consensus sequences")
    return core_intervals_on_reference


def generate_core_intervals_on_reference(
    padded_alignments_paths: dict[int, Path],
    index_consensusIDs: dict[int, set[str]],
    consensus_json_index_partitioned: dict[int, dict[str, int]],
    consensus_json_paths: dict[int, Path],
    threads: int,
    batch_size: int = 100
) -> dict[int, dict[str, list[tuple[int, str, int, int]]]]:
    """Generate core intervals on reference for all alignment files in parallel.
    
    Args:
        padded_alignments_paths: Dictionary mapping file indices to serialized alignment TSV paths
        index_consensusIDs: Dictionary mapping file indices to sets of consensus IDs
        consensus_json_index_partitioned: Partitioned index mapping partition_idx -> consensusID -> byte_offset
        consensus_json_paths: Dictionary mapping file indices to consensus TSV paths
        threads: Number of parallel workers to use
        batch_size: Number of alignments to process per batch
        
    Returns:
        Partitioned dictionary mapping file_index -> consensusID -> list of (alignment_idx, ref_name, core_start_on_ref, core_end_on_ref)
    """
    log.info(f"Generating core intervals on reference using {threads} parallel workers...")
    
    # Prepare arguments for parallel processing
    process_args = []
    for i in range(threads):
        process_args.append((
            i,
            padded_alignments_paths[i],
            index_consensusIDs[i],
            consensus_json_index_partitioned[i],
            consensus_json_paths[i],
            batch_size
        ))
    
    # Process files in parallel with progress bar tracking completed files
    with mp.Pool(processes=threads) as pool:
        results = list(tqdm(
            pool.starmap(_process_alignment_file_for_core_intervals, process_args),
            total=threads,
            desc="Processing alignment files"
        ))
    
    # Keep partitioned structure: file_idx -> consensusID -> intervals
    partitioned_core_intervals: dict[int, dict[str, list[tuple[int, str, int, int]]]] = {}
    for i, result_dict in enumerate(results):
        partitioned_core_intervals[i] = result_dict
    
    total_consensuses = sum(len(intervals) for intervals in partitioned_core_intervals.values())
    log.info(f"Generated core intervals for {total_consensuses} consensus sequences across {threads} files.")
    
    return partitioned_core_intervals


def _process_partition_for_trf_overlaps(
    partition_idx: int,
    core_intervals: dict[str, list[tuple[int, str, int, int]]],
    input_trf: Path,
    reference: Path,
    tmp_dir: Path
) -> dict[str, list[tuple[str, int, int, int]]]:
    """Process a single partition to find TRF overlaps for core intervals using bedtools.
    
    Args:
        partition_idx: Index of the partition being processed
        core_intervals: Dictionary mapping consensusID to list of (alignment_idx, ref_name, core_start, core_end)
        input_trf: Path to TRF bed file
        reference: Path to reference genome (for creating genome file)
        tmp_dir: Temporary directory for intermediate files
        
    Returns:
        Dictionary mapping consensusID to list of TRF intervals (chrom, start, end, repeat_id)
        All TRF overlaps for a consensusID are collected in a single flat list
        The repeat_id is the 0-indexed line number from the TRF bed file
    """
    log.debug(f"Worker {partition_idx}: Finding TRF overlaps for {len(core_intervals)} consensus sequences")
    
    # Write core intervals to BED file with consensusID as 4th column
    bed_file = tmp_dir / f"core_intervals_{partition_idx}.bed"
    with open(bed_file, 'w') as f:
        for consensusID, intervals in core_intervals.items():
            for aln_idx, ref_name, core_start, core_end in intervals:
                print(ref_name, core_start, core_end, consensusID, sep='\t', file=f)
    
    # Sort the BED file
    sorted_bed_file = tmp_dir / f"core_intervals_{partition_idx}.sorted.bed"
    genome_file = tmp_dir / f"genome_{partition_idx}.txt"
    util.genome_file_for_bedtools(reference=reference, output=genome_file)
    
    cmd_sort = f"bedtools sort -i {bed_file} -g {genome_file}"
    with open(sorted_bed_file, 'w') as f:
        subprocess.check_call(shlex.split(cmd_sort), stdout=f)
    
    # Add repeat_id to TRF bed file (0-indexed line number)
    trf_with_id = tmp_dir / f"trf_with_id_{partition_idx}.bed"
    with open(input_trf, 'r') as fin, open(trf_with_id, 'w') as fout:
        for repeat_id, line in enumerate(fin):
            fields = line.rstrip().split('\t')
            if len(fields) >= 3:
                print(fields[0], fields[1], fields[2], repeat_id, sep='\t', file=fout)
    
    # Use bedtools intersect to find overlapping TRFs
    intersect_output = tmp_dir / f"trf_overlaps_{partition_idx}.bed"
    cmd_intersect = f"bedtools intersect -a {sorted_bed_file} -b {trf_with_id} -wa -wb"
    with open(intersect_output, 'w') as f:
        subprocess.check_call(shlex.split(cmd_intersect), stdout=f)
    
    # Parse the intersect output and build the result structure
    # Format: chrom start end consensusID chrom_trf start_trf end_trf repeat_id ...
    result: dict[str, list[tuple[str, int, int, int]]] = {}
    
    # Initialize result structure with empty lists for each consensusID
    for consensusID in core_intervals.keys():
        result[consensusID] = []
    
    # Parse intersect output
    with open(intersect_output, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if len(fields) < 8:
                continue
            
            # Extract core interval info
            consensusID = fields[3]
            
            # Extract TRF interval info (now with repeat_id as 4th column)
            trf_chrom = fields[4]
            trf_start = int(fields[5])
            trf_end = int(fields[6])
            repeat_id = int(fields[7])
            
            # Add TRF overlap to consensusID's list
            if consensusID in result:
                trf_tuple = (trf_chrom, trf_start, trf_end, repeat_id)
                # Avoid duplicates (same TRF might overlap multiple alignments)
                if trf_tuple not in result[consensusID]:
                    result[consensusID].append(trf_tuple)
    
    log.debug(f"Worker {partition_idx}: Processed {len(result)} consensus sequences")
    return result


def find_trf_overlaps_for_core_intervals(
    core_intervals_on_reference: dict[int, dict[str, list[tuple[int, str, int, int]]]],
    input_trf: Path,
    reference: Path,
    threads: int,
    tmp_dir: Path
) -> dict[int, dict[str, list[tuple[str, int, int, int]]]]:
    """Find TRF overlaps for core intervals on reference in parallel using bedtools.
    
    This function processes each partition in parallel to find overlapping TRF intervals
    for each core interval on the reference.
    
    Args:
        core_intervals_on_reference: Partitioned dict mapping partition_idx -> consensusID -> 
                                     list of (alignment_idx, ref_name, core_start, core_end)
        input_trf: Path to TRF bed file
        reference: Path to reference genome
        threads: Number of parallel workers to use
        tmp_dir: Temporary directory for intermediate files
        
    Returns:
        Partitioned dict mapping partition_idx -> consensusID -> list of TRF intervals
        For each consensusID, all overlapping TRF intervals are collected as (chrom, start, end, repeat_id) tuples
        The repeat_id is the 0-indexed line number from the TRF bed file
    """
    log.info(f"Finding TRF overlaps for core intervals using {threads} parallel workers...")
    
    # Prepare arguments for parallel processing
    process_args = []
    for partition_idx, core_intervals in core_intervals_on_reference.items():
        process_args.append((
            partition_idx,
            core_intervals,
            input_trf,
            reference,
            tmp_dir
        ))
    
    # Process partitions in parallel
    with mp.Pool(processes=threads) as pool:
        results = list(tqdm(
            pool.starmap(_process_partition_for_trf_overlaps, process_args),
            total=len(process_args),
            desc="Finding TRF overlaps"
        ))
    
    # Rebuild partitioned structure maintaining partition indices
    trf_overlaps_partitioned: dict[int, dict[str, list[tuple[str, int, int, int]]]] = {}
    for partition_idx in core_intervals_on_reference.keys():
        trf_overlaps_partitioned[partition_idx] = results[partition_idx]
    
    total_consensuses = sum(len(overlaps) for overlaps in trf_overlaps_partitioned.values())
    log.info(f"Found TRF overlaps for {total_consensuses} consensus sequences across {threads} partitions.")
    
    return trf_overlaps_partitioned


def build_fasta_index(fasta_paths: dict[int, Path]) -> dict[str, tuple[int, int, int]]:
    """Build an index mapping consensusID to (file_index, byte_offset, sequence_length) for FASTA files.
    
    Args:
        fasta_paths: Dictionary mapping file indices to FASTA file paths
        
    Returns:
        Dictionary mapping consensusID to (file_index, byte_offset, sequence_length)
    """
    index: dict[str, tuple[int, int, int]] = {}
    
    for file_idx, path in fasta_paths.items():
        with open(path, 'r') as f:
            current_id = None
            header_offset = 0
            seq_length = 0
            
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        index[current_id] = (file_idx, header_offset, seq_length)
                    
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get ID after '>'
                    header_offset = f.tell() - len(line) - 1  # Position of header line
                    seq_length = 0
                else:
                    # Accumulate sequence length
                    seq_length += len(line)
            
            # Save last sequence
            if current_id is not None:
                index[current_id] = (file_idx, header_offset, seq_length)
    
    return index


def build_fasta_index_partitioned(fasta_paths: dict[int, Path]) -> dict[int, dict[str, tuple[int, int]]]:
    """Build partitioned FASTA indices directly, partitioned by file.
    
    Args:
        fasta_paths: Dictionary mapping file indices to FASTA file paths
        
    Returns:
        Dictionary mapping partition_idx -> consensusID -> (byte_offset, sequence_length)
        Note: file_index is not stored since it equals partition_idx
    """
    partitioned_index: dict[int, dict[str, tuple[int, int]]] = {}
    
    for file_idx, path in fasta_paths.items():
        partitioned_index[file_idx] = {}
        with open(path, 'r') as f:
            current_id = None
            header_offset = 0
            seq_length = 0
            
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        partitioned_index[file_idx][current_id] = (header_offset, seq_length)
                    
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get ID after '>'
                    header_offset = f.tell() - len(line) - 1  # Position of header line
                    seq_length = 0
                else:
                    # Accumulate sequence length
                    seq_length += len(line)
            
            # Save last sequence
            if current_id is not None:
                partitioned_index[file_idx][current_id] = (header_offset, seq_length)
    
    return partitioned_index


def get_sequence_by_id(
    consensusID: str,
    index: dict[str, tuple[int, int, int]],
    fasta_paths: dict[int, Path]
) -> str:
    """Retrieve a sequence from indexed FASTA files by consensusID.
    
    Args:
        consensusID: The consensus ID to retrieve
        index: Index mapping consensusID to (file_index, byte_offset, sequence_length)
        fasta_paths: Dictionary mapping file indices to FASTA file paths
        
    Returns:
        The sequence as a string
    """
    if consensusID not in index:
        raise KeyError(f"consensusID {consensusID} not found in FASTA index")
    
    file_idx, offset, _seq_length = index[consensusID]
    
    with open(fasta_paths[file_idx], 'r') as f:
        f.seek(offset)
        
        # Skip header line
        f.readline()
        
        # Read sequence lines until next header or EOF
        sequence_parts = []
        while True:
            pos = f.tell()
            line = f.readline()
            
            if not line or line.startswith('>'):
                break
            
            sequence_parts.append(line.rstrip('\n'))
        
        return ''.join(sequence_parts)


def get_sequences_batch_by_ids(
    consensusIDs: list[str],
    index: dict[str, tuple[int, int, int]],
    fasta_paths: dict[int, Path]
) -> dict[str, str]:
    """Retrieve multiple sequences from indexed FASTA files in a batched manner.
    
    This is more efficient than calling get_sequence_by_id repeatedly
    because it minimizes file open/close operations.
    
    Args:
        consensusIDs: List of consensus IDs to retrieve
        index: Index mapping consensusID to (file_index, byte_offset, sequence_length)
        fasta_paths: Dictionary mapping file indices to FASTA file paths
        
    Returns:
        Dictionary mapping consensusID to sequence string
    """
    # Group consensusIDs by file_index to minimize file operations
    ids_by_file: dict[int, list[tuple[str, int, int]]] = {}
    
    for consensusID in consensusIDs:
        if consensusID not in index:
            log.warning(f"consensusID {consensusID} not found in FASTA index")
            continue
        
        file_idx, offset, seq_length = index[consensusID]
        if file_idx not in ids_by_file:
            ids_by_file[file_idx] = []
        ids_by_file[file_idx].append((consensusID, offset, seq_length))
    
    # Sort by offset for sequential reads (faster than random access)
    for file_idx in ids_by_file:
        ids_by_file[file_idx].sort(key=lambda x: x[1])
    
    # Retrieve sequences
    results: dict[str, str] = {}
    
    for file_idx, id_offset_pairs in ids_by_file.items():
        with open(fasta_paths[file_idx], 'r') as f:
            for consensusID, offset, _seq_length in id_offset_pairs:
                f.seek(offset)
                
                # Skip header line
                f.readline()
                
                # Read sequence lines until next header or EOF
                sequence_parts = []
                while True:
                    line = f.readline()
                    
                    if not line or line.startswith('>'):
                        break
                    
                    sequence_parts.append(line.rstrip('\n'))
                
                results[consensusID] = ''.join(sequence_parts)
    
    return results


def get_sequences_batch_from_partition(
    consensusIDs: list[str],
    index_partition: dict[str, tuple[int, int]],
    fasta_path: Path
) -> dict[str, str]:
    """Retrieve multiple sequences from a single partition's FASTA file.
    
    More efficient than get_sequences_batch_by_ids when all sequences are in one file.
    
    Args:
        consensusIDs: List of consensus IDs to retrieve
        index_partition: Partitioned index mapping consensusID to (byte_offset, sequence_length)
        fasta_path: Path to the FASTA file for this partition
        
    Returns:
        Dictionary mapping consensusID to sequence string
    """
    # Collect all IDs with their offsets
    id_offset_pairs: list[tuple[str, int, int]] = []
    
    for consensusID in consensusIDs:
        if consensusID not in index_partition:
            log.warning(f"consensusID {consensusID} not found in partition FASTA index")
            continue
        
        offset, seq_length = index_partition[consensusID]
        id_offset_pairs.append((consensusID, offset, seq_length))
    
    # Sort by offset for sequential reads (faster than random access)
    id_offset_pairs.sort(key=lambda x: x[1])
    
    # Retrieve sequences
    results: dict[str, str] = {}
    
    with open(fasta_path, 'r') as f:
        for consensusID, offset, _seq_length in id_offset_pairs:
            f.seek(offset)
            
            # Skip header line
            f.readline()
            
            # Read sequence lines until next header or EOF
            sequence_parts = []
            while True:
                line = f.readline()
                
                if not line or line.startswith('>'):
                    break
                
                sequence_parts.append(line.rstrip('\n'))
            
            results[consensusID] = ''.join(sequence_parts)
    
    return results


def get_consensus_batch_by_ids(
    consensusIDs: list[str],
    index: dict[str, tuple[int, int]],
    consensus_paths: dict[int, Path]
) -> dict[str, consensus_class.Consensus]:
    """Retrieve multiple consensus objects by ID in a batched manner.
    
    This is more efficient than calling get_consensus_by_id repeatedly
    because it minimizes file open/close operations.
    """
    # Group consensusIDs by file_index to minimize file operations
    ids_by_file: dict[int, list[tuple[str, int]]] = {}
    
    for consensusID in consensusIDs:
        if consensusID not in index:
            log.warning(f"consensusID {consensusID} not found in index")
            continue
        
        file_idx, offset = index[consensusID]
        if file_idx not in ids_by_file:
            ids_by_file[file_idx] = []
        ids_by_file[file_idx].append((consensusID, offset))
    
    # Sort by offset for sequential reads (faster than random access)
    for file_idx in ids_by_file:
        ids_by_file[file_idx].sort(key=lambda x: x[1])
    
    # Retrieve consensus objects
    results: dict[str, consensus_class.Consensus] = {}
    
    for file_idx, id_offset_pairs in ids_by_file.items():
        with open(consensus_paths[file_idx], 'r') as f:
            for consensusID, offset in id_offset_pairs:
                f.seek(offset)
                line = f.readline()
                _, serialized_data = line.strip().split('\t', 1)
                results[consensusID] = cattrs.structure(
                    json.loads(serialized_data), 
                    consensus_class.Consensus
                )
    
    return results


def get_consensusIDs_by_partition(
    index_alignments_partitioned: dict[int, dict[str, list[int]]]
    ) -> dict[int, set[str]]:
    """Extract consensusIDs from partitioned alignment index.
    
    Args:
        index_alignments_partitioned: Partitioned index mapping partition_idx -> consensusID -> list of offsets
        
    Returns:
        Dictionary mapping partition_idx to set of consensusIDs in that partition
    """
    consensus_ids: dict[int, set[str]] = {}
    for partition_idx, partition_dict in index_alignments_partitioned.items():
        consensus_ids[partition_idx] = set(partition_dict.keys())
    return consensus_ids


def get_alignments_for_consensus_from_partition(
    consensusID: str,
    index_partition: dict[str, list[int]],
    alignments_path: Path
) -> list[datatypes.Alignment]:
    """Retrieve all alignments for a consensusID from a single partition file.
    
    Args:
        consensusID: The consensus ID to retrieve alignments for
        index_partition: Partitioned index mapping consensusID to list of byte_offsets
        alignments_path: Path to the alignments TSV file for this partition
        
    Returns:
        List of Alignment objects for this consensusID
    """
    if consensusID not in index_partition:
        return []
    
    alignments: list[datatypes.Alignment] = []
    offsets = index_partition[consensusID]
    
    with open(alignments_path, 'r') as f:
        for offset in offsets:
            f.seek(offset)
            line = f.readline()
            fields = line.rstrip().split('\t', 1)
            if len(fields) < 2:
                log.warning(f"Malformed alignment line at offset {offset} in partition file {alignments_path}")
                continue
            
            # Deserialize alignment
            alignment = cattrs.structure(
                json.loads(fields[1]),
                datatypes.Alignment
            )
            alignments.append(alignment)
    
    return alignments


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
        # if output_consensus_to_reference_alignments is None:
        #     output_consensus_to_reference_alignments = (
        #         Path(tmp_dir) / "consensus_to_reference.bam"
        #     )

        output_consensus_to_reference_alignments = Path(
            output_consensus_to_reference_alignments
        )

        # the following part will be multi-threaded. All extracted core sequences will be written to {threads} tmp files for parallel processig
        core_sequences_paths:dict[int,Path] = {i: Path(tmp_dir) / f"core_sequences.{i}.fasta" for i in range(threads)}
        padded_sequences_paths:dict[int,Path] = {i: Path(tmp_dir) / f"padded_sequences.{i}.fasta" for i in range(threads)}
        consensus_json_paths:dict[int,Path] = {i: Path(tmp_dir) / f"consensus.{i}.tsv" for i in range(threads)}
        
        # write all core sequences and all padded sequences to fasta files for parallel processing
        index_consensusIDs = write_tmp_core_sequences_files_for_parallel_processing(
            input_consensus_container_results=input_consensus_container_results,
            core_sequences_paths=core_sequences_paths,
            padded_sequences_paths=padded_sequences_paths,
            consensus_paths=consensus_json_paths
        )
        
        # Build indices for FASTA files (samtools faidx for external tools)
        log.info("Building samtools FASTA indices...")
        index_tmp_fastas(paths=list(core_sequences_paths.values()) + list(padded_sequences_paths.values()))
        
        # Build partitioned in-memory indices directly (more efficient than building global then partitioning)
        log.info("Building partitioned in-memory indices for core and padded sequences...")
        core_sequences_index_partitioned = build_fasta_index_partitioned(core_sequences_paths)
        padded_sequences_index_partitioned = build_fasta_index_partitioned(padded_sequences_paths)
        consensus_json_index_partitioned = build_consensus_tsv_index_partitioned(consensus_json_paths)
        
        total_core_seqs = sum(len(partition) for partition in core_sequences_index_partitioned.values())
        total_padded_seqs = sum(len(partition) for partition in padded_sequences_index_partitioned.values())
        log.info(f"Indexed {total_core_seqs} core sequences and {total_padded_seqs} padded sequences across {len(core_sequences_index_partitioned)} partitions.")
        
        tmp_padded_alignments_path : Path = Path(tmp_dir) / "padded_alignments_all.bam"
        
        util.align_reads_with_minimap(
            reference=input_reference,
            bamout=tmp_padded_alignments_path,
            reads=[str(padded_sequences_paths[i]) for i in range(len(padded_sequences_paths))],
            tech="map-ont",
            threads=threads,
            aln_args=" --secondary=no ")
        
        # Parse alignments and write partitioned datatypes.Alignment objects to TSV files
        padded_alignments_paths: dict[int, Path] = {
            i: Path(tmp_dir) / f"padded_alignments.{i}.tsv" 
            for i in range(threads)
        }
        
        write_partitioned_alignments(
            input_bam=tmp_padded_alignments_path,
            output_paths=padded_alignments_paths,
            index_consensusIDs=index_consensusIDs,
            samplename=samplename
        )
        
        # Build partitioned index for random access to alignments by consensusID
        log.info("Building partitioned alignment index for random access...")
        index_alignments_partitioned: dict[int, dict[str, list[int]]] = build_alignments_index_partitioned(padded_alignments_paths)
        total_consensus_with_alignments = sum(len(partition) for partition in index_alignments_partitioned.values())
        log.info(f"Indexed {total_consensus_with_alignments} consensus sequences with alignments across {len(index_alignments_partitioned)} partitions.")
        
        # now iterate the alignments in each padded_alignments_paths
        # for each alignment get the corresponding consensus object from the consensus_json_paths
        # and get its core interval on the reference coordinate system via consensus_class.get_consensus_core_alignment_interval_on_reference
        # this should be done in parallel over each alignments file        
        core_intervals_on_reference: dict[int, dict[str, list[tuple[int, str, int, int]]]] = generate_core_intervals_on_reference(
            padded_alignments_paths=padded_alignments_paths,
            index_consensusIDs=index_consensusIDs,
            consensus_json_index_partitioned=consensus_json_index_partitioned,
            consensus_json_paths=consensus_json_paths,
            threads=threads,
            batch_size=100
        )
        # now we need to find the intersecting trf intervals for each core interval on reference
        # Find TRF overlaps for all core intervals in parallel using bedtools
        trf_overlaps_on_reference: dict[int, dict[str, list[tuple[str, int, int, int]]]] = find_trf_overlaps_for_core_intervals(
            core_intervals_on_reference=core_intervals_on_reference,
            input_trf=Path(input_trf),
            reference=Path(input_reference),
            threads=threads,
            tmp_dir=Path(tmp_dir)
        )
        
        consensusIDs_by_partition = get_consensusIDs_by_partition(index_alignments_partitioned=index_alignments_partitioned)

        # now we can parse SV signals from the consensus alignments
        # fully annotate them
        # and parse them to SV patterns (per consensus and all its alignments)
        partition_index = 0
        # this example is for a single partition index for development purposes and later parallelization
        # get all consensusIDs per partition index from index_alignments
        
        core_intervals_on_reference_partition = core_intervals_on_reference[partition_index]
        trf_overlaps_on_reference_partition = trf_overlaps_on_reference[partition_index]
        
        # Use partitioned indices for this partition to minimize memory footprint
        core_sequences_index_partition = core_sequences_index_partitioned[partition_index]
        padded_sequences_index_partition = padded_sequences_index_partitioned[partition_index]
        consensus_json_index_partition = consensus_json_index_partitioned[partition_index]
        alignments_index_partition = index_alignments_partitioned[partition_index]
        
        batchsize = 100
        for i in range(0, len(consensusIDs_by_partition[partition_index]), batchsize):
            consensusIDs_batch = list(consensusIDs_by_partition[partition_index])[i:i+batchsize]
            # get all consensus objects for this batch at once using partitioned index
            # Need to reconstruct the full index format for get_consensus_batch_by_ids
            consensus_index_for_batch = {
                consensusID: (partition_index, consensus_json_index_partition[consensusID])
                for consensusID in consensusIDs_batch
                if consensusID in consensus_json_index_partition
            }
            consensus_objs = get_consensus_batch_by_ids(
                consensus_paths=consensus_json_paths,
                consensusIDs=consensusIDs_batch,
                index=consensus_index_for_batch
            )
            # get all core sequences for this batch at once using partitioned index
            core_sequences_batch = get_sequences_batch_from_partition(
                consensusIDs=consensusIDs_batch,
                index_partition=core_sequences_index_partition,
                fasta_path=core_sequences_paths[partition_index]
            )
            # get all padded sequences for this batch at once using partitioned index
            padded_sequences_batch = get_sequences_batch_from_partition(
                consensusIDs=consensusIDs_batch,
                index_partition=padded_sequences_index_partition,
                fasta_path=padded_sequences_paths[partition_index]
            )
            # process each consensusID in the batch
            for consensusID in consensusIDs_batch:
                consensus_obj = consensus_objs[consensusID]
                core_sequence = core_sequences_batch[consensusID]
                padded_sequence = padded_sequences_batch[consensusID] # might not be necessary. will see
                trf_overlaps = trf_overlaps_on_reference_partition.get(consensusID, [])
                alignments = get_alignments_for_consensus_from_partition(
                    consensusID=consensusID,
                    index_partition=alignments_index_partition,
                    alignments_path=padded_alignments_paths[partition_index]
                )
                
                # Get the core intervals for this consensus and create lookup by alignment index
                core_intervals_with_idx = core_intervals_on_reference_partition.get(consensusID, [])
                core_interval_by_idx = {
                    aln_idx: (ref_name, core_start, core_end)
                    for aln_idx, ref_name, core_start, core_end in core_intervals_with_idx
                }
                
                mergedSVs: list[datatypes.MergedSVSignal] = []
                svPrimitives: list[SVprimitives.SVprimitive] = []
                for alignment_idx, alignment in enumerate(alignments):
                    # Match by alignment index
                    if alignment_idx not in core_interval_by_idx:
                        log.warning(
                            f"No core interval found for alignment {alignment_idx} of {consensusID}"
                        )
                        continue
                                        
                    # Filter trf intervals to those overlapping with this alignment's reference
                    filtered_trf_overlaps = [
                        (trf[1], trf[2], trf[3]) 
                        for trf in trf_overlaps 
                        if trf[0] == alignment.reference_name
                    ]
                    if consensus_obj.consensus_padding is None:
                        raise ValueError(f"Consensus {consensusID} has no consensus_padding information!")
                    log.debug(f"Parsing SV signals from consensus {consensusID} alignment to {alignment.reference_name}\nCore interval on reference: {str(consensus_obj.consensus_padding.consensus_interval_on_reference)}\nCore interval on sequence with padding: {consensus_obj.consensus_padding.consensus_interval_on_sequence_with_padding}\nNumber of TRF overlaps for this reference: {len(filtered_trf_overlaps)}")
                    mergedSVs.extend(consensus_align_lib.parse_sv_signals_from_consensus(
                        samplename=samplename,
                        consensus_sequence=core_sequence,
                        interval_core=core_interval_by_idx[alignment_idx][1:3],
                        pysam_alignment=alignment.to_pysam(),
                        trf_intervals=filtered_trf_overlaps,
                        reference_name=alignment.reference_name,
                        min_bnd_size=min_bnd_size,
                        min_signal_size=min_signal_size)
                    )
                    # now SV signals are parsed for this alignment
                    svPrimitives.extend(SVprimitives.generate_SVprimitives(
                        samplename=samplename,
                        mergedSVs=mergedSVs,
                        consensus=consensus_obj,
                        consensus_alignment=alignment,
                        alignmentID=alignment_idx,
                        core_interval=core_interval_by_idx[alignment_idx])
                    )
                # now the sv patterns can be created from the svPrimitives
                
                    
                    
                    
                # now all signals can be extracted, annotated and the sv patterns generated
                # use consensus_align_lib.parse_sv_signals_from_consensus to get the signals (mergedSVs)
                # consensusAlignment.proto_svs = merged_svs
                
        
        
        
        
        
        
        
        # consensus_align_lib.parse_sv_signals_from_consensus
        
        log.debug(
            f"loading alignments from {output_consensus_to_reference_alignments}..."
        )
        # consensusID (read name of aligned consensus) : list[Alignment]
        consensus_alignments: dict[str, list[datatypes.Alignment]] = load_alignments(
            path_alignments=output_consensus_to_reference_alignments, parse_DNA=False
        )
        n_alignments = sum(
            len(alignments) for alignments in consensus_alignments.values()
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
            f"{sum(len(alignments) for alignments in dict_alignments.values())} consensus alignments converted to ConsensusAlignment objects."
        )

        log.info("Caching pysam alignments for tracing back core intervals...")
        pysam_cache: dict[int, AlignedSegment] = {}
        # fill the pysam cache
        for _consensusID, consensusAlignments in tqdm(dict_alignments.items()):
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
        consensus_alignments = {}  # has no real use anymore, free memory

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
        _consensusIDs = list(dict_alignments.keys())  # FIXME: unused?
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
        svp_types = {svp.sv_type for svp in svPrimitives}
        log.warning(f"{len(svPrimitives)} SV primitives of types found: {svp_types}")

        # svPrimitives = [svp for svp in svPrimitives if len(svp.genotypeMeasurements.supporting_reads_start) > 0]

        log.info("Adding local depth of coverage information to the SV primitives...")
        # add coverages to svPrimitives
        # covtrees_dict:dict[str, IntervalTree] = covtree.covtree(path_db=covtrees_db)
        # add_depth_to_svPrimitives_genotypes_inplace(svPrimitives=svPrimitives, _covtree=covtrees_dict)

        svPatterns = svPrimitives_to_svPatterns(
            SVprimitives=svPrimitives, max_del_size=max_del_size
        )
        svp_types = {svp.get_sv_type() for svp in svPatterns}
        log.warning(
            f"After svPrimitives_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
        )

        svPatterns = add_consensus_sequence_and_size_distortions_to_svPatterns(
            consensus_results_path=input_consensus_container_results,
            svPatterns=svPatterns,
            distance_scale=distance_scale,
            falloff=falloff,
        )
        svp_types = {svp.get_sv_type() for svp in svPatterns}
        log.warning(
            f"After add_consensus_sequence_and_size_distortions_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
        )

        svPatterns = add_reference_sequence_to_svPatterns(
            svPatterns=svPatterns, reference_sequence=input_reference
        )
        svp_types = {svp.get_sv_type() for svp in svPatterns}
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
