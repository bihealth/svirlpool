import argparse
import gzip
import json
import logging
import multiprocessing as mp
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Generator
from copy import deepcopy
from io import TextIOWrapper
from pathlib import Path

import attrs
import cattrs
from Bio import SeqIO
from Bio.SeqRecord import Seq, SeqRecord
from intervaltree import IntervalTree
from pysam import AlignedSegment, AlignmentFile, AlignmentHeader
from tqdm import tqdm

from ..util import datatypes, util
from . import SVpatterns, SVprimitives, consensus_align_lib, consensus_class
from .SVpatterns import converter as svpattern_converter

log = logging.getLogger(__name__)


@attrs.define
class SVPatternProcessingParams:
    """Parameters for processing consensus objects to SV patterns in parallel."""

    svPatterns_path: Path
    consensusIDs: list[str]
    consensus_json_index_partition: dict[str, int]
    consensus_json_paths: dict[int, Path]
    partition_index: int
    padded_alignments_path: Path
    core_intervals_on_reference: dict[str, list[tuple[int, str, int, int, int, int]]]
    alignments_index_partition: dict[str, list[int]]
    trf_overlaps_on_reference_partition: dict[str, list[tuple[str, int, int, int]]]
    input_reference: Path
    samplename: str
    min_bnd_size: int
    min_signal_size: int
    max_del_size: int
    distance_scale: float
    falloff: float
    max_fourrelations_gap_size: int
    signal_loss_log_dir: Path | None = None
    logfile: Path | None = None
    batchsize: int = 100
    log_level: str = "INFO"


def parse_crs_container_results(
    path: Path | str,
) -> Generator[consensus_class.CrsContainerResult, None, None]:
    with open(path, "r") as f:
        for line in f:
            crs_container_result = json.loads(line)
            yield cattrs.structure(
                crs_container_result, consensus_class.CrsContainerResult
            )


def trf_to_interval_tree(input_trf: Path) -> dict[str, IntervalTree]:
    """produces a dict chrom:IntervalTree"""
    # load the trf intervals to an intervaltree and set their index as their names
    trf_intervals = {}
    open_func = gzip.open if str(input_trf).endswith(".gz") else open
    mode = "rt" if str(input_trf).endswith(".gz") else "r"
    with open_func(input_trf, mode) as f:
        for i, line in enumerate(f):
            chrom, start, end = line.strip().split()[0:3]
            if chrom not in trf_intervals:
                trf_intervals[chrom] = IntervalTree()
            trf_intervals[chrom][int(start) : int(end)] = i
    return trf_intervals


def svPrimitives_to_svPatterns(
    SVprimitives: list[SVprimitives.SVprimitive],
    max_del_size: int,
    max_fourrelations_gap_size: int = 500_000,
) -> list[SVpatterns.SVpatternType]:
    """
    Converts a list of SVprimitives to a list of SVpatterns.
    This is done by parsing the SVprimitives and creating SVpatterns from them.
    IMPORTANT: SVprimitives are in the order of their alignment position on the consensus sequence.
    """

    if len(SVprimitives) == 0:
        log.warning("No SVprimitives provided. Returning empty list of SVpatterns.")
        return []

    # Group SVprimitives by consensusID
    grouped_svprimitives = {}
    for svp in SVprimitives:
        if svp.consensusID not in grouped_svprimitives:
            grouped_svprimitives[svp.consensusID] = []
        grouped_svprimitives[svp.consensusID].append(svp)

    # Process each group of SVprimitives
    svpatterns = []
    for group in grouped_svprimitives.values():
        svpatterns.extend(
            SVpatterns.parse_SVprimitives_to_SVpatterns(
                SVprimitives=group,
                max_del_size=max_del_size,
                max_fourrelations_gap_size=max_fourrelations_gap_size,
            )
        )

    return svpatterns


def add_consensus_sequence_and_size_distortions_to_svPatterns(
    consensus_objects: dict[str, consensus_class.Consensus],
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
    for consensusID, consensus in consensus_objects.items():
        if consensusID not in svPatterns_dict:
            continue

        for idx, svp in svPatterns_dict[consensusID]:
            if idx in processed_indices:
                continue

            # Process the SVpattern - create a copy to avoid modifying original
            processed_svp = svp  # reference is sufficient

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
            if isinstance(processed_svp, SVpatterns.SVpatternInsertion) or isinstance(
                processed_svp, SVpatterns.SVpatternInversion
            ):
                processed_svp.set_sequence(
                    processed_svp.get_sequence_from_consensus(consensus=consensus)
                )
            elif isinstance(processed_svp, SVpatterns.SVpatternDeletion):
                continue  # deletions don't have an alt sequence — this is expected, not data loss
            elif isinstance(processed_svp, SVpatterns.SVpatternSingleBreakend):
                processed_svp.set_clipped_sequence(
                    sequence=processed_svp.get_clipped_sequence_from_consensus(
                        consensus=consensus
                    ),
                )
                processed_svp.set_context_sequence(
                    processed_svp.get_context_sequence_from_consensus(
                        consensus=consensus
                    )
                )
            elif isinstance(processed_svp, SVpatterns.SVpatternAdjacency):
                processed_svp.set_all_sequence_contexts(
                    processed_svp.get_sequence_contexts_from_consensus(
                        consensus=consensus
                    )
                )
                processed_svp.set_sequence(
                    processed_svp.get_sequence_from_consensus(consensus=consensus)
                )
            else:
                try:
                    id_str: str = processed_svp._log_id()
                except Exception as e:
                    log.error(f"Error generating log ID for SVpattern: {e}")
                    id_str = f"consensusID={processed_svp.consensusID}|crID={processed_svp.consensusID.split('.')[0]}|type={type(processed_svp)}"
                log.debug(
                    f"DROPPED::add_consensus_sequence_and_size_distortions_to_svPatterns::(SV type not supported), svpattern={id_str}"
                )
                # _loss_logger = get_signal_loss_logger()
                # _loss_logger.log_skipped(
                #     stage="add_consensus_sequence_and_size_distortions_to_svPatterns",
                #     consensusID=consensusID,
                #     reason="unsupported_SVpattern_type_for_alt_sequence",
                #     details={"svpattern_type": str(type(processed_svp).__name__)},
                # )
                log.warning(
                    f"SVpattern type {type(processed_svp)} is not supported. Skipping alt sequence assignment."
                )

            # Add to output and mark as processed
            output_svPatterns.append(processed_svp)
            processed_indices.add(idx)

    # Add any unprocessed SVpatterns (those without matching consensus)
    for i, svp in enumerate(svPatterns):
        if i not in processed_indices and svp is not None:
            output_svPatterns.append(svp)

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
                for region_key in regions.keys():
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
                        if type(x) is SVpatterns.SVpatternDeletion:
                            x.set_sequence(str(record.seq))
                        elif type(x) in (
                            SVpatterns.SVpatternInversion,
                            SVpatterns.SVpatternInversionDuplication,
                            SVpatterns.SVpatternInversionDeletion,
                        ):
                            x.set_deleted_sequence(
                                sequence=str(record.seq), write_complexity=False
                            )
                        else:
                            raise ValueError(
                                f"SVpattern type {type(x)} not supported for reference sequence assignment. This should not happen. Please check the input data."
                            )
                        output.append(x)  # append the modified copy
                        svp = None  # remove reference to original to free memory
                    regions[region_key].clear()
                else:
                    raise ValueError(
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


def write_consensus_files_for_parallel_processing(
    input_consensus_container_results: Path | str,
    consensus_paths: dict[int, Path],
    padded_sequences_paths: dict[int, Path],
) -> dict[int, set[str]]:
    """Writes temporary core sequence fasta, padded sequence fasta, and consensus json tsv files for parallel processing."""
    # write the core sequences to fasta files
    # 1) open all fasta writers
    n_files = len(consensus_paths)
    padded_writers: dict[int, TextIOWrapper] = {}
    for i in range(n_files):
        padded_writers[i] = open(padded_sequences_paths[i], "w")
    consensus_writers: dict[int, TextIOWrapper] = {}
    for i in range(n_files):
        consensus_writers[i] = open(consensus_paths[i], "w")
    # 2) iterate over all consensus sequences and write them in a round-robin fashion
    writer_idx: int = 0
    index_consensusIDs: dict[int, set[str]] = {
        i: set() for i in range(n_files)
    }  # to find the correct fasta file for each consensusID

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

            # consensus.consensus_sequence = ""  # remove the sequences to save space and time
            consensus.consensus_padding.sequence = (
                ""  # remove the sequences to save space and time
            )
            print(
                str(consensus.ID),
                json.dumps(consensus.unstructure()),
                sep="\t",
                file=consensus_writers[writer_idx],
            )

            # core interval would be consensus.consensus_padding.consensus_interval_on_sequence_with_padding
            index_consensusIDs[writer_idx].add(consensus.ID)
        writer_idx = (writer_idx + 1) % n_files  # round-robin
    # close all fasta writers
    for i in range(n_files):
        padded_writers[i].close()
        consensus_writers[i].close()
    return index_consensusIDs


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
        writers[i] = open(output_paths[i], "w")

    # Parse BAM file and write alignments to partitioned files
    with AlignmentFile(str(input_bam), "rb") as bam_in:
        for pysam_aln in bam_in:
            consensusID = pysam_aln.query_name

            if consensusID not in consensusID_to_index:
                raise ValueError(
                    f"consensusID {consensusID} not found in index_consensusIDs. "
                    "This should not happen. Please check the input data."
                )

            if pysam_aln.reference_end is None:
                raise ValueError(
                    f"Alignment for consensusID {consensusID} has no reference_end. "
                    "This should not happen. Please check the input data."
                )

            qstart, qend = util.get_interval_on_read_in_region(
                a=pysam_aln,
                start=pysam_aln.reference_start,
                end=pysam_aln.reference_end,
            )

            # Convert pysam alignment to datatypes.Alignment
            annotations = {"qstart": qstart, "qend": qend, "qmid": (qstart + qend) // 2}
            alignment = datatypes.Alignment.from_pysam(
                aln=pysam_aln, samplename=samplename, annotations=annotations
            )

            # Write to appropriate partition file: consensusID\tserialized_alignment
            # This ensures that all alignments of one consensusID go to the same file
            idx = consensusID_to_index[consensusID]
            print(
                consensusID,
                qstart,
                qend,
                json.dumps(alignment.unstructure()),
                sep="\t",
                file=writers[idx],
            )

    # Close all writers
    for writer in writers.values():
        writer.close()


def build_alignments_index(
    alignment_paths: dict[int, Path],
) -> dict[str, list[tuple[int, int]]]:
    """Build an index mapping consensusID to list of (file_index, byte_offset) for alignments.

    Args:
        alignment_paths: Dictionary mapping file indices to alignment TSV paths

    Returns:
        Dictionary mapping consensusID to list of (file_index, byte_offset) tuples,
        one entry per alignment of that consensusID
    """
    index: dict[str, list[tuple[int, int]]] = {}

    for file_idx, path in alignment_paths.items():
        with open(path, "r") as f:
            while True:
                offset = f.tell()
                line = f.readline()

                if not line:  # EOF
                    break

                consensusID = line.split("\t", 1)[0]
                if consensusID not in index:
                    index[consensusID] = []
                index[consensusID].append((file_idx, offset))

    return index


def build_alignments_index_partitioned(
    alignment_paths: dict[int, Path],
) -> dict[int, dict[str, list[int]]]:
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
        with open(path, "r") as f:
            while True:
                offset = f.tell()
                line = f.readline()

                if not line:  # EOF
                    break

                consensusID = line.split("\t", 1)[0]
                if consensusID not in partitioned_index[file_idx]:
                    partitioned_index[file_idx][consensusID] = []
                partitioned_index[file_idx][consensusID].append(offset)

    return partitioned_index


def build_consensus_tsv_index(
    consensus_paths: dict[int, Path],
) -> dict[str, tuple[int, int]]:
    """Build an index mapping consensusID to (file_index, byte_offset)."""
    index: dict[str, tuple[int, int]] = {}

    for file_idx, path in consensus_paths.items():
        with open(path, "r") as f:
            while True:
                offset = f.tell()
                line = f.readline()

                if not line:  # EOF
                    break

                consensusID = line.split("\t", 1)[0]
                index[consensusID] = (file_idx, offset)

    return index


def build_consensus_tsv_index_partitioned(
    consensus_paths: dict[int, Path],
) -> dict[int, dict[str, int]]:
    """Build partitioned indices mapping consensusID to byte_offset, partitioned by file.

    Args:
        consensus_paths: Dictionary mapping file indices to consensus TSV paths

    Returns:
        Dictionary mapping partition_idx -> consensusID -> byte_offset
    """
    partitioned_index: dict[int, dict[str, int]] = {}

    for file_idx, path in consensus_paths.items():
        partitioned_index[file_idx] = {}
        with open(path, "r") as f:
            while True:
                offset = f.tell()
                line = f.readline()

                if not line:  # EOF
                    break

                consensusID = line.split("\t", 1)[0]
                partitioned_index[file_idx][consensusID] = offset

    return partitioned_index


def _process_alignment_file_for_core_intervals(
    file_idx: int,
    padded_alignments_path: Path,
    index_consensusIDs: set[str],
    consensus_json_index_partition: dict[str, int],
    consensus_json_path: Path,
    batch_size: int = 35,
    log_level: str = "INFO",
) -> dict[str, list[tuple[int, str, int, int, int, int]]]:
    """Process a single alignment file to extract core intervals on reference.

    This function is designed to be run in parallel for each alignment file.

    Args:
        file_idx: Index of the file being processed
        padded_alignments_path: Path to the serialized alignments TSV file
        index_consensusIDs: Set of consensus IDs expected in this file
        consensus_json_index_partition: Partitioned index mapping consensusID to byte_offset for this partition only
        consensus_json_path: Path to the consensus TSV file for this partition
        batch_size: Number of consensusIDs for which alignments are loaded to process per batch
        log_level: Logging level string (e.g. "DEBUG", "INFO") – must be passed explicitly
                   because multiprocessing workers do not inherit the parent's logging config.

    Returns:
        Dictionary mapping consensusID to list of (alignment_idx, ref_name, core_start_on_ref, core_end_on_ref, core_read_start, core_read_end)
        The alignment_idx corresponds to the order in the alignment file for this consensusID
    """
    # Configure logging for this worker.  On Linux mp.Pool uses fork, so workers
    # inherit the parent's named-logger objects with whatever levels were set there.
    # basicConfig(force=True) resets the root logger but does NOT reset named child
    # loggers.  We therefore also call setLevel on the named logger explicitly.
    worker_log_level = getattr(logging, log_level, logging.INFO)
    logging.basicConfig(
        level=worker_log_level,
        format=f"[Worker {file_idx}] %(asctime)s - %(name)s - %(levelname)s - %(message)s",
        stream=sys.stderr,
        force=True,
    )
    worker_log = logging.getLogger(__name__)
    worker_log.setLevel(worker_log_level)

    # core ref chr, core ref start, core ref end, core read start, core read end
    core_intervals_on_reference: dict[str, list[tuple[int, str, int, int, int, int]]] = {}

    worker_log.debug(f"Worker {file_idx}: Processing {padded_alignments_path}")

    # Batch structure: maps consensusID to list of (qstart, qend, alignment)
    current_batch: dict[str, list[tuple[int, int, datatypes.Alignment]]] = {}

    all_seen_consensusIDs = (
        set()
    )  # to track all consensusIDs seen in the file for validation

    def process_batch(
        batch: dict[str, list[tuple[int, int, datatypes.Alignment]]],
    ) -> None:
        """Process a batch of alignments grouped by consensusID."""
        if not batch:
            return

        # Get all unique consensusIDs in this batch
        batch_consensusIDs = list(batch.keys())
        all_seen_consensusIDs.update(batch_consensusIDs)

        # Reconstruct the index format expected by get_consensus_batch_by_ids
        batch_index = {
            cid: (file_idx, consensus_json_index_partition[cid])
            for cid in batch_consensusIDs
            if cid in consensus_json_index_partition
        }
        consensus_objs = get_consensus_batch_by_ids(
            consensusIDs=batch_consensusIDs,
            index=batch_index,
            consensus_paths={file_idx: consensus_json_path},
        )

        # Process each consensusID's alignments
        for consensusID, alignments_list in batch.items():
            consensus_obj = consensus_objs[consensusID]

            # Sort alignments by (qstart, qend) to ensure consistent ordering
            alignments_list.sort(key=lambda x: (x[0], x[1]))

            # DEBUG START
            # --- write alignment(s) using the same container format the test loader expects ---
            if consensus_obj.ID == "11.0":
                import gzip as _gzip
                import json as _json

                path_debug_save = Path(
                    "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.11.alignments.json.gz"
                )
                # convert the pysam AlignedSegment -> datatypes.Alignment -> plain dict
                alns_unstructured = [aln[2].unstructure() for aln in alignments_list]
                # append to existing alignments if the file already exists
                if path_debug_save.exists():
                    with _gzip.open(path_debug_save, "rt") as f:
                        existing = _json.load(f)
                    alns_unstructured = existing["alignments"] + alns_unstructured
                with _gzip.open(path_debug_save, "wt") as f:
                    _json.dump({"alignments": alns_unstructured}, f, indent=2)

                # write the consensus object exactly as the test loader expects
                path_debug_save_consensus = Path(
                    "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/consensus_class/INV.11.consensus.json.gz"
                )
                with _gzip.open(path_debug_save_consensus, "wt") as f:
                    _json.dump(consensus_obj.unstructure(), f, indent=2)
            # DEBUG END


            # Process each alignment with its index
            for aln_idx, (_qstart, _qend, alignment) in enumerate(alignments_list):
                # Convert to pysam for processing
                pysam_aln = alignment.to_pysam()

                # Get core interval on reference (already traced back from consensus coordinates)
                ref_name, core_start_on_ref, core_end_on_ref = (
                    consensus_class.get_consensus_core_alignment_interval_on_reference(
                        consensus=consensus_obj, alignment=pysam_aln
                    )
                )

                # Skip alignments that don't overlap with the core
                # - this is a bit tricky. only one alignment might be overlapping, so that many times
                # - a false filter would be reported.
                if core_start_on_ref == core_end_on_ref:
                    continue

                if consensusID not in core_intervals_on_reference:
                    core_intervals_on_reference[consensusID] = []

                # Store with alignment index for matching
                if consensus_obj.consensus_padding.consensus_interval_on_sequence_with_padding is None:
                    raise ValueError(
                        f"_process_alignment_file_for_core_intervals::process_batch: Consensus {consensusID} has no consensus_interval_on_sequence_with_padding."
                    )
                core_intervals_on_reference[consensusID].append((
                    aln_idx,  # Alignment index within this consensusID (sorted by qstart, qend)
                    ref_name,
                    min(core_start_on_ref, core_end_on_ref),
                    max(core_start_on_ref, core_end_on_ref),
                    consensus_obj.consensus_padding.consensus_interval_on_sequence_with_padding[0],
                    consensus_obj.consensus_padding.consensus_interval_on_sequence_with_padding[1],
                ))

    # Read serialized alignments from TSV file
    # File is sorted by: consensusID, qstart, qend
    with open(padded_alignments_path, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t", 4)
            if len(fields) < 4:
                continue

            consensusID = fields[0]
            qstart = int(fields[1])
            qend = int(fields[2])

            if consensusID not in index_consensusIDs:
                raise ValueError(
                    f"consensusID {consensusID} found in alignment file but not in "
                    f"index_consensusIDs for index {file_idx}. This should not happen."
                )

            # Deserialize alignment
            alignment = cattrs.structure(json.loads(fields[3]), datatypes.Alignment)

            # Check if we need to start a new batch (when we hit a new consensusID)
            if consensusID not in current_batch:
                # If batch is full, process it before adding new consensusID
                if len(current_batch) >= batch_size:
                    process_batch(current_batch)
                    current_batch.clear()

                current_batch[consensusID] = []

            # Add alignment to current batch
            current_batch[consensusID].append((qstart, qend, alignment))

        # Process remaining batch
        if current_batch:
            process_batch(current_batch)

    # # finally compare all all_seen_consensusIDs with core_intervals_on_reference.keys()
    # # and report all that have beed dropped due to missing alignments (should be very few or none)
    # missing_alignments_consensusIDs = set(core_intervals_on_reference.keys()) - all_seen_consensusIDs
    # if missing_alignments_consensusIDs:
    #     log.warning(
    #         f"Worker {file_idx}: Found {len(missing_alignments_consensusIDs)} consensusIDs with core intervals but no alignments. This should not happen. Please check the input data. Missing consensusIDs: {missing_alignments_consensusIDs}"
    #     )

    worker_log.debug(
        f"Worker {file_idx}: Processed {len(core_intervals_on_reference)} consensus sequences:"
    )
    if worker_log.isEnabledFor(logging.DEBUG):
        for consensusID, intervals in core_intervals_on_reference.items():
            for aln_idx, ref_name, core_start_on_ref, core_end_on_ref, core_read_start, core_read_end in intervals:
                worker_log.debug(
                    f"Worker {file_idx}: consensusID {consensusID}, alignment_idx {aln_idx}: core interval on reference/ region={ref_name}:{core_start_on_ref}-{core_end_on_ref}, core interval on read={core_read_start}-{core_read_end}"
                )

    return core_intervals_on_reference


def generate_core_intervals_on_reference(
    padded_alignments_paths: dict[int, Path],
    index_consensusIDs: dict[int, set[str]],
    consensus_json_index_partitioned: dict[int, dict[str, int]],
    consensus_json_paths: dict[int, Path],
    threads: int,
    batch_size: int = 100,
) -> dict[int, dict[str, list[tuple[int, str, int, int, int, int]]]]:
    """Generate core intervals on reference for all alignment files in parallel.

    Args:
        padded_alignments_paths: Dictionary mapping file indices to serialized alignment TSV paths
        index_consensusIDs: Dictionary mapping file indices to sets of consensus IDs
        consensus_json_index_partitioned: Partitioned index mapping partition_idx -> consensusID -> byte_offset
        consensus_json_paths: Dictionary mapping file indices to consensus TSV paths
        threads: Number of parallel workers to use
        batch_size: Number of alignments to process per batch

    Returns:
        Partitioned dictionary mapping file_index -> consensusID -> list of (alignment_idx, ref_name, core_start_on_ref, core_end_on_ref, core_read_start, core_read_end)
    """
    log.info(
        f"Generating core intervals on reference using {threads} parallel workers..."
    )

    # Prepare arguments for parallel processing
    # Pass the effective log level so workers configure their loggers correctly.
    effective_log_level = logging.getLevelName(log.getEffectiveLevel())
    process_args = []
    for i in range(threads):
        process_args.append((
            i,
            padded_alignments_paths[i],
            index_consensusIDs[i],
            consensus_json_index_partitioned[i],
            consensus_json_paths[i],
            batch_size,
            effective_log_level,
        ))

    # Process files in parallel with progress bar tracking completed files
    with mp.Pool(processes=threads) as pool:
        results = list(
            tqdm(
                pool.starmap(_process_alignment_file_for_core_intervals, process_args),
                total=threads,
                desc="Processing alignment files",
            )
        )

    # Keep partitioned structure: file_idx -> consensusID -> intervals
    partitioned_core_intervals: dict[
        int, dict[str, list[tuple[int, str, int, int, int, int]]]
    ] = {}
    for i, result_dict in enumerate(results):
        partitioned_core_intervals[i] = result_dict

    total_consensuses = sum(
        len(intervals) for intervals in partitioned_core_intervals.values()
    )
    log.info(
        f"Generated core intervals for {total_consensuses} consensus sequences across {threads} files."
    )
    # print all core intervals of each consensusID and alignmend index if the log level is DEBUG
    if log.isEnabledFor(logging.DEBUG):
        for partition_idx, consensus_intervals in partitioned_core_intervals.items():
            for consensusID, intervals in consensus_intervals.items():
                for aln_idx, ref_name, core_start_on_ref, core_end_on_ref, core_read_start, core_read_end in intervals:
                    log.debug(
                        f"TRANSFORMED::generate_core_intervals_on_reference: Partition {partition_idx}, consensusID {consensusID}, crID={consensusID.split('.')[0]}, alignment_idx {aln_idx}: core interval on reference/ region={ref_name}:{core_start_on_ref}-{core_end_on_ref}, core interval on read={core_read_start}-{core_read_end}"
                    )
    
    
    # to log every case of data loss, we can compare the consensusIDs for which we generated core intervals with the consensusIDs in the consensus_json_index_partitioned
    # and log any consensusID for which we generated core intervals but is not present in the consensus_json_index_partitioned
    # (which means we won't be able to find the consensus object for that consensusID later.
    # first create a set of all input consensusIDs from the consensus_json_index_partitioned
    all_input_consensusIDs = set()
    for partition in consensus_json_index_partitioned.values():
        all_input_consensusIDs.update(partition.keys())
    # then create a set of all consensusIDs for which we generated core intervals
    generated_core_interval_consensusIDs = set()
    for partition in partitioned_core_intervals.values():
        generated_core_interval_consensusIDs.update(partition.keys())
    # find the difference
    missing_consensusIDs = generated_core_interval_consensusIDs - all_input_consensusIDs
    if missing_consensusIDs:
        for consensusID in missing_consensusIDs:
            log.debug(
                f"DROPPED::generate_core_intervals_on_reference:(couldn't find intervals), consensusID {consensusID}."
            )

    return partitioned_core_intervals


def _process_partition_for_trf_overlaps(
    partition_idx: int,
    core_intervals: dict[str, list[tuple[int, str, int, int, int, int]]],
    input_trf: Path,
    reference: Path,
    tmp_dir: Path,
) -> dict[str, list[tuple[str, int, int, int]]]:
    """Process a single partition to find TRF overlaps for core intervals using bedtools.

    Args:
        partition_idx: Index of the partition being processed
        core_intervals: Dictionary mapping consensusID to list of (alignment_idx, ref_name, core_start_on_ref, core_end_on_ref, core_read_start, core_read_end)
        input_trf: Path to TRF bed file
        reference: Path to reference genome (for creating genome file)
        tmp_dir: Temporary directory for intermediate files

    Returns:
        Dictionary mapping consensusID to list of TRF intervals (chrom, start, end, repeat_id)
        All TRF overlaps for a consensusID are collected in a single flat list
        The repeat_id is the 0-indexed line number from the TRF bed file
    """
    log.debug(
        f"Worker {partition_idx}: Finding TRF overlaps for {len(core_intervals)} consensus sequences"
    )

    # Write core intervals to BED file with consensusID as 4th column
    bed_file = tmp_dir / f"core_intervals_{partition_idx}.bed"
    with open(bed_file, "w") as f:
        for consensusID, intervals in core_intervals.items():
            for (
                _,
                ref_name,
                core_start,
                core_end,
                _core_read_start,
                _core_read_end,
            ) in intervals:  # alignment index not needed here.
                print(ref_name, core_start, core_end, consensusID, sep="\t", file=f)

    # Sort the BED file
    sorted_bed_file = tmp_dir / f"core_intervals_{partition_idx}.sorted.bed"
    genome_file = tmp_dir / f"genome_{partition_idx}.txt"
    util.genome_file_for_bedtools(reference=reference, output=genome_file)

    cmd_sort = f"bedtools sort -i {bed_file} -g {genome_file}"
    with open(sorted_bed_file, "w") as f:
        subprocess.check_call(shlex.split(cmd_sort), stdout=f)

    # Add repeat_id to TRF bed file (0-indexed line number)
    trf_with_id = tmp_dir / f"trf_with_id_{partition_idx}.bed"
    open_func = gzip.open if str(input_trf).endswith(".gz") else open
    mode = "rt" if str(input_trf).endswith(".gz") else "r"
    with open_func(input_trf, mode) as fin, open(trf_with_id, "w") as fout:
        for repeat_id, line in enumerate(fin):
            fields = line.rstrip().split("\t")
            if len(fields) >= 3:
                print(fields[0], fields[1], fields[2], repeat_id, sep="\t", file=fout)

    # Use bedtools intersect to find overlapping TRFs
    intersect_output = tmp_dir / f"trf_overlaps_{partition_idx}.bed"
    cmd_intersect = f"bedtools intersect -a {sorted_bed_file} -b {trf_with_id} -wa -wb"
    with open(intersect_output, "w") as f:
        subprocess.check_call(shlex.split(cmd_intersect), stdout=f)

    # Parse the intersect output and build the result structure
    # Format: chrom start end consensusID chrom_trf start_trf end_trf repeat_id ...
    result: dict[str, list[tuple[str, int, int, int]]] = {}

    # Initialize result structure with empty lists for each consensusID
    for consensusID in core_intervals.keys():
        result[consensusID] = []

    # Parse intersect output
    with open(intersect_output, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
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
    core_intervals_on_reference: dict[int, dict[str, list[tuple[int, str, int, int, int, int]]]],
    input_trf: Path,
    reference: Path,
    threads: int,
    tmp_dir: Path,
) -> dict[int, dict[str, list[tuple[str, int, int, int]]]]:
    """Find TRF overlaps for core intervals on reference in parallel using bedtools.

    This function processes each partition in parallel to find overlapping TRF intervals
    for each core interval on the reference.

    Args:
        core_intervals_on_reference: Partitioned dict mapping partition_idx -> consensusID ->
                                     list of (alignment_idx, ref_name, core_start_on_ref, core_end_on_ref, core_read_start, core_read_end)
        input_trf: Path to TRF bed file
        reference: Path to reference genome
        threads: Number of parallel workers to use
        tmp_dir: Temporary directory for intermediate files

    Returns:
        Partitioned dict mapping partition_idx -> consensusID -> list of TRF intervals
        For each consensusID, all overlapping TRF intervals are collected as (chrom, start, end, repeat_id) tuples
        The repeat_id is the 0-indexed line number from the TRF bed file
    """
    log.info(
        f"Finding TRF overlaps for core intervals using {threads} parallel workers..."
    )

    # Prepare arguments for parallel processing
    process_args = []
    for partition_idx, core_intervals in core_intervals_on_reference.items():
        process_args.append((
            partition_idx,
            core_intervals,
            input_trf,
            reference,
            tmp_dir,
        ))

    # Process partitions in parallel
    with mp.Pool(processes=threads) as pool:
        results = list(
            tqdm(
                pool.starmap(_process_partition_for_trf_overlaps, process_args),
                total=len(process_args),
                desc="Finding TRF overlaps",
            )
        )

    # Rebuild partitioned structure maintaining partition indices
    trf_overlaps_partitioned: dict[int, dict[str, list[tuple[str, int, int, int]]]] = {}
    for partition_idx in core_intervals_on_reference.keys():
        trf_overlaps_partitioned[partition_idx] = results[partition_idx]

    total_consensuses = sum(
        len(overlaps) for overlaps in trf_overlaps_partitioned.values()
    )
    log.info(
        f"Found TRF overlaps for {total_consensuses} consensus sequences across {threads} partitions."
    )

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
        with open(path, "r") as f:
            current_id = None
            header_offset = 0
            seq_length = 0

            while True:
                # Store position before reading line
                pos = f.tell()
                line = f.readline()

                if not line:  # EOF
                    break

                line = line.rstrip("\n")
                if line.startswith(">"):
                    # Save previous sequence if exists
                    if current_id is not None:
                        index[current_id] = (file_idx, header_offset, seq_length)

                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get ID after '>'
                    header_offset = pos  # Position of header line
                    seq_length = 0
                else:
                    # Accumulate sequence length
                    seq_length += len(line)

            # Save last sequence
            if current_id is not None:
                index[current_id] = (file_idx, header_offset, seq_length)

    return index


def build_fasta_index_partitioned(
    fasta_paths: dict[int, Path],
) -> dict[int, dict[str, tuple[int, int]]]:
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
        with open(path, "r") as f:
            current_id = None
            header_offset = 0
            seq_length = 0

            while True:
                # Store position before reading line
                pos = f.tell()
                line = f.readline()

                if not line:  # EOF
                    break

                line = line.rstrip("\n")
                if line.startswith(">"):
                    # Save previous sequence if exists
                    if current_id is not None:
                        partitioned_index[file_idx][current_id] = (
                            header_offset,
                            seq_length,
                        )

                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get ID after '>'
                    header_offset = pos  # Position of header line
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
    fasta_paths: dict[int, Path],
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

    with open(fasta_paths[file_idx], "r") as f:
        f.seek(offset)

        # Skip header line
        f.readline()

        # Read sequence lines until next header or EOF
        sequence_parts = []
        while True:
            line = f.readline()

            if not line or line.startswith(">"):
                break

            sequence_parts.append(line.rstrip("\n"))

        return "".join(sequence_parts)


def get_sequences_batch_by_ids(
    consensusIDs: list[str],
    index: dict[str, tuple[int, int, int]],
    fasta_paths: dict[int, Path],
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
        with open(fasta_paths[file_idx], "r") as f:
            for consensusID, offset, _seq_length in id_offset_pairs:
                f.seek(offset)

                # Skip header line
                f.readline()

                # Read sequence lines until next header or EOF
                sequence_parts = []
                while True:
                    line = f.readline()

                    if not line or line.startswith(">"):
                        break

                    sequence_parts.append(line.rstrip("\n"))

                results[consensusID] = "".join(sequence_parts)

    return results


def get_sequences_batch_from_partition(
    consensusIDs: list[str],
    index_partition: dict[str, tuple[int, int]],
    fasta_path: Path,
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

    with open(fasta_path, "r") as f:
        for consensusID, offset, _seq_length in id_offset_pairs:
            f.seek(offset)

            # Skip header line
            f.readline()

            # Read sequence lines until next header or EOF
            sequence_parts = []
            while True:
                line = f.readline()

                if not line or line.startswith(">"):
                    break

                sequence_parts.append(line.rstrip("\n"))

            results[consensusID] = "".join(sequence_parts)

    return results


def get_consensus_batch_by_ids(
    consensusIDs: list[str],
    index: dict[str, tuple[int, int]],
    consensus_paths: dict[int, Path],
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
        with open(consensus_paths[file_idx], "r") as f:
            for consensusID, offset in id_offset_pairs:
                f.seek(offset)
                line = f.readline()
                _, serialized_data = line.strip().split("\t", 1)
                results[consensusID] = cattrs.structure(
                    json.loads(serialized_data), consensus_class.Consensus
                )

    return results


def get_consensusIDs_by_partition(
    index_alignments_partitioned: dict[int, dict[str, list[int]]],
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
    consensusID: str, index_partition: dict[str, list[int]], alignments_path: Path
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

    with open(alignments_path, "r") as f:
        for offset in offsets:
            f.seek(offset)
            line = f.readline()
            fields = line.rstrip().split("\t")
            if (
                len(fields) < 4
            ):  # expected consensusID, qstart, qend, serialized_alignment
                log.warning(
                    f"Malformed alignment line at offset {offset} in partition file {alignments_path}. Expected: consensusID, qstart, qend, serialized_alignment, but got: {fields}"
                )
                continue

            # Deserialize alignment
            alignment = cattrs.structure(json.loads(fields[-1]), datatypes.Alignment)
            alignments.append(alignment)
    return alignments


def write_sequence_fastas_to_file(
    padded_sequences_paths: dict[int, Path], output_consensus_fasta: Path | str
) -> None:
    # cat all consensus sequences fasta files and write them to output_consensus_fasta
    with open(output_consensus_fasta, "w") as outfile:
        for i in range(len(padded_sequences_paths)):
            with open(padded_sequences_paths[i], "r") as infile:
                shutil.copyfileobj(infile, outfile)
                # Ensure newline between files to preserve FASTA format
                if (
                    i < len(padded_sequences_paths) - 1
                ):  # Don't add extra newline after last file
                    outfile.write("\n")
    log.info(
        f"Concatenated {len(padded_sequences_paths)} FASTA files to {output_consensus_fasta}"
    )


def _process_consensus_objects_to_svPatterns(params: SVPatternProcessingParams):
    """Process consensus objects to SV patterns for a single partition.

    Args:
        params: SVPatternProcessingParams containing all necessary parameters
    """
    # Configure logging for this worker process

    worker_log_level = getattr(logging, params.log_level, logging.INFO)
    worker_format = f"[Worker {params.partition_index}] %(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        level=worker_log_level,
        format=worker_format,
        stream=sys.stderr,
        force=True,
    )

    # Add file handler so worker logs are written to the log file
    if params.logfile is not None:
        file_handler = logging.FileHandler(params.logfile)
        file_handler.setLevel(worker_log_level)
        file_handler.setFormatter(logging.Formatter(worker_format))
        logging.getLogger().addHandler(file_handler)

    # Initialise per-worker signal loss logger
    # if params.signal_loss_log_dir is not None:
    #     _worker_loss_log = params.signal_loss_log_dir / f"signal_loss.partition_{params.partition_index}.tsv"
    #     init_signal_loss_logger(output_path=_worker_loss_log)
    # else:
    #     init_signal_loss_logger()  # fallback to Python logger

    with open(
        params.svPatterns_path, "w"
    ) as svp_outfile:  # tmp output file for this partition
        for i in range(0, len(params.consensusIDs), params.batchsize):
            consensusIDs_batch = list(params.consensusIDs)[i : i + params.batchsize]
            # get all consensus objects for this batch at once using partitioned index
            # Need to reconstruct the full index format for get_consensus_batch_by_ids
            consensus_index_for_batch = {
                consensusID: (
                    params.partition_index,
                    params.consensus_json_index_partition[consensusID],
                )
                for consensusID in consensusIDs_batch
                if consensusID in params.consensus_json_index_partition
            }
            consensus_objs = get_consensus_batch_by_ids(
                consensus_paths=params.consensus_json_paths,
                consensusIDs=consensusIDs_batch,
                index=consensus_index_for_batch,
            )
            svPrimitives: list[SVprimitives.SVprimitive] = []
            # process each consensusID in the batch
            for consensusID in consensusIDs_batch:
                consensus_obj = consensus_objs[consensusID]
                trf_overlaps = params.trf_overlaps_on_reference_partition.get(
                    consensusID, []
                )
                alignments = get_alignments_for_consensus_from_partition(
                    consensusID=consensusID,
                    index_partition=params.alignments_index_partition,
                    alignments_path=params.padded_alignments_path,
                )

                # sort alignments by alignment.annotations.qstart
                for aln in alignments:
                    if aln.annotations is None or "qstart" not in aln.annotations:
                        raise ValueError(
                            f"Alignment for consensus {consensusID} missing 'qstart' in annotations!"
                        )
                alignments.sort(
                    key=lambda aln: aln.annotations["qmid"]
                )  # alignments need to be sorted by mid point of query sequence
                # to have all break ends in the right order for sv pattern generation, since tips can overlap by a few bases

                # Get the core intervals for this consensus and create lookup by alignment index
                core_intervals_with_idx = params.core_intervals_on_reference.get(
                    consensusID, []
                )
                core_interval_by_idx = {
                    aln_idx: (ref_name, core_start, core_end, core_read_start, core_read_end)
                    for aln_idx, ref_name, core_start, core_end, core_read_start, core_read_end in core_intervals_with_idx
                }

                for alignment_idx, alignment in enumerate(
                    alignments
                ):  # alignment_idx order on consensus sequence!
                    # Match by alignment index

                    if alignment_idx not in core_interval_by_idx:
                        log.debug(
                            f"DROPPED::_process_consensus_objects_to_svPatterns:(no core interval for alignment), consensusID {consensusID}, alignment_idx {alignment_idx}; alignment_region={alignment.reference_name}:{alignment.reference_start}-{alignment.reference_end}."
                            f"\ncore intervals on reference: {' '.join([f'region={ci[1]}:{ci[2]}-{ci[3]} (core interval={ci[4]}-{ci[5]})' for ci in core_intervals_with_idx])}"
                        )
                        continue
                    
                    # Filter trf intervals to those overlapping with this alignment's reference
                    filtered_trf_overlaps = [
                        (trf[1], trf[2], trf[3])
                        for trf in trf_overlaps
                        if trf[0] == alignment.reference_name
                    ]
                    if consensus_obj.consensus_padding is None:
                        raise ValueError(
                            f"Consensus {consensusID} has no consensus_padding information!"
                        )
                    log.debug(
                        f"Parsing SV signals from consensus {consensusID} alignment to {alignment.reference_name}\nCore interval on reference: {core_interval_by_idx[alignment_idx][1:3]}\nCore interval on sequence with padding: {consensus_obj.consensus_padding.consensus_interval_on_sequence_with_padding}\nNumber of TRF overlaps for this reference: {len(filtered_trf_overlaps)}"
                    )
                    mergedSVs = consensus_align_lib.parse_sv_signals_from_consensus(
                        samplename=params.samplename,
                        consensus_sequence=consensus_obj.consensus_sequence,
                        interval_core=consensus_obj.consensus_padding.consensus_interval_on_sequence_with_padding,
                        pysam_alignment=alignment.to_pysam(),
                        trf_intervals=filtered_trf_overlaps,
                        reference_name=alignment.reference_name,
                        min_bnd_size=params.min_bnd_size,
                        min_signal_size=params.min_signal_size,
                    )  # sorted in order of consensus sequence

                    log.debug(
                        f"Consensus {consensusID} alignment {alignment_idx}: Found {len(mergedSVs)} merged SV signals"
                    )

                    # now SV signals are parsed for this alignment
                    svPrimitives.extend(
                        SVprimitives.generate_SVprimitives(
                            samplename=params.samplename,
                            mergedSVs=mergedSVs,
                            consensus=consensus_obj,
                            consensus_alignment=alignment,
                            alignmentID=alignment_idx,
                            core_interval=core_interval_by_idx[alignment_idx],
                        )
                    )
                    log.debug(
                        f"Consensus {consensusID} alignment {alignment_idx}: Generated {len(svPrimitives)} SV primitives so far"
                    )
                # now the sv patterns can be created from the svPrimitives
            log.debug(
                f"Batch {i}: Total {len(svPrimitives)} SV primitives from {len(consensusIDs_batch)} consensus sequences"
            )
            svPatterns = svPrimitives_to_svPatterns(
                SVprimitives=svPrimitives,
                max_del_size=params.max_del_size,
                max_fourrelations_gap_size=params.max_fourrelations_gap_size,
            )
            svp_types = {svp.get_sv_type() for svp in svPatterns}
            log.warning(
                f"After svPrimitives_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
            )
            # TODO: rewrite to not use input_consensus_container_results but the ready-made data files
            svPatterns = add_consensus_sequence_and_size_distortions_to_svPatterns(
                consensus_objects=consensus_objs,
                svPatterns=svPatterns,
                distance_scale=params.distance_scale,
                falloff=params.falloff,
            )
            svp_types = {svp.get_sv_type() for svp in svPatterns}
            log.warning(
                f"After add_consensus_sequence_and_size_distortions_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
            )

            # todo: re-write necessary?
            svPatterns = add_reference_sequence_to_svPatterns(
                svPatterns=svPatterns, reference_sequence=str(params.input_reference)
            )
            svp_types = {svp.get_sv_type() for svp in svPatterns}
            log.warning(
                f"After add_reference_sequence_to_svPatterns. {len(svPatterns)} SV patterns of types found: {svp_types}"
            )
            # write all svPatterns to the tmp_svPrimitives_path
            for svp in svPatterns:
                print(
                    json.dumps(svpattern_converter.unstructure(svp)), file=svp_outfile
                )
            log.debug(
                f"Processed batch {i} of consensus sequences for partition {params.partition_index}"
            )

    # Flush and close the per-worker signal loss logger before the worker exits.
    # Without this, mp.Pool workers may terminate (via os._exit) without
    # flushing Python I/O buffers, losing all buffered log entries.
    # _worker_logger = get_signal_loss_logger()
    # _worker_logger.flush()
    # _worker_logger.close()


def sort_partitioned_alignments(
    padded_alignments_paths: dict[int, Path],
) -> None:
    """Sort partitioned alignment TSV files by consensusID, qstart, qend.
    Those are the first three columns in the tsv file."""
    cmd_sort = "sort -k1,1 -k2,2n -k3,3n "  # sort by consensusID, qstart, qend
    for partition_idx, path in padded_alignments_paths.items():
        sorted_path = path.with_suffix(".sorted.tsv")
        with open(sorted_path, "w") as outfile:
            with open(path, "r") as infile:
                subprocess.check_call(
                    shlex.split(cmd_sort), stdin=infile, stdout=outfile
                )
        # replace original file with sorted file
        shutil.move(sorted_path, path)
        log.info(f"Sorted alignments in partition {partition_idx} and updated file.")


def svPatterns_from_consensus_sequences(
    samplename: str,
    input_consensus_container_results: Path | str,
    input_reference: Path | str,
    input_trf: Path | str,
    output_consensus_fasta: Path | str,
    output_consensus_to_reference_alignments: Path | str,
    output_svPatterns: Path | str,
    threads: int,
    min_signal_size: int,
    min_bnd_size: int,
    max_del_size: int,
    distance_scale: float,
    falloff: float,
    max_fourrelations_gap_size: int = 500_000,
    tmp_dir_path: Path | str | None = None,
    dont_merge_horizontally: bool = False,
    batchsize: int = 100,
    signal_loss_log: Path | str | None = None,
    logfile: Path | str | None = None,
    log_level: str = "INFO",
) -> None:
    # write all consensus sequences to a fasta file and align to the target reference
    # write all alignments to a file
    if dont_merge_horizontally:
        log.warning(
            "dont_merge_horizontally is set to True. This will not merge SVs within the same consensus alignment. This is not recommended for general use cases."
        )

    # Initialise the signal loss logger for the main process
    _signal_loss_log_dir: Path | None = None
    _logfile: Path | None = Path(logfile) if logfile is not None else None
    # if signal_loss_log is not None:
    #     _signal_loss_log_path = Path(signal_loss_log)
    #     _signal_loss_log_dir = _signal_loss_log_path.parent
    #     init_signal_loss_logger(output_path=_signal_loss_log_path)
    #     log.info(f"Signal loss log will be written to {_signal_loss_log_path}")
    # else:
    #     init_signal_loss_logger()  # fallback to Python logger

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
        padded_sequences_paths: dict[int, Path] = {
            i: Path(tmp_dir) / f"padded_sequences.{i}.fasta" for i in range(threads)
        }
        consensus_json_paths: dict[int, Path] = {
            i: Path(tmp_dir) / f"consensus.{i}.tsv" for i in range(threads)
        }

        # write all core sequences and all padded sequences to fasta files for parallel processing
        index_consensusIDs = write_consensus_files_for_parallel_processing(
            input_consensus_container_results=input_consensus_container_results,
            padded_sequences_paths=padded_sequences_paths,
            consensus_paths=consensus_json_paths,
        )
        write_sequence_fastas_to_file(
            output_consensus_fasta=output_consensus_fasta,
            padded_sequences_paths=padded_sequences_paths,
        )

        # Build partitioned in-memory indices directly (more efficient than building global then partitioning)
        log.info(
            "Building partitioned in-memory indices for core and padded sequences..."
        )
        padded_sequences_index_partitioned = build_fasta_index_partitioned(
            padded_sequences_paths
        )
        consensus_json_index_partitioned = build_consensus_tsv_index_partitioned(
            consensus_json_paths
        )

        total_padded_seqs = sum(
            len(partition) for partition in padded_sequences_index_partitioned.values()
        )
        log.info(
            f"Indexed {total_padded_seqs} padded sequences across {len(padded_sequences_index_partitioned)} partitions."
        )

        util.align_reads_with_minimap(
            reference=input_reference,
            bamout=output_consensus_to_reference_alignments,
            reads=[
                str(padded_sequences_paths[i])
                for i in range(len(padded_sequences_paths))
            ],
            tech="map-ont",
            threads=threads,
            aln_args=" --secondary=no",
        )

        # Parse alignments and write partitioned datatypes.Alignment objects to TSV files
        padded_alignments_paths: dict[int, Path] = {
            i: Path(tmp_dir) / f"padded_alignments.{i}.tsv" for i in range(threads)
        }

        write_partitioned_alignments(
            input_bam=output_consensus_to_reference_alignments,
            output_paths=padded_alignments_paths,
            index_consensusIDs=index_consensusIDs,
            samplename=samplename,
        )
        # sort partitioned alignments by consenusID, qstart, qend
        sort_partitioned_alignments(padded_alignments_paths=padded_alignments_paths)

        # Build partitioned index for random access to alignments by consensusID
        log.info("Building partitioned alignment index for random access...")
        index_alignments_partitioned: dict[int, dict[str, list[int]]] = (
            build_alignments_index_partitioned(padded_alignments_paths)
        )
        total_consensus_with_alignments = sum(
            len(partition) for partition in index_alignments_partitioned.values()
        )
        log.info(
            f"Indexed {total_consensus_with_alignments} consensus sequences with alignments across {len(index_alignments_partitioned)} partitions."
        )

        # now iterate the alignments in each padded_alignments_paths
        # for each alignment get the corresponding consensus object from the consensus_json_paths
        # and get its core interval on the reference coordinate system via consensus_class.get_consensus_core_alignment_interval_on_reference
        # this should be done in parallel over each alignments file
        core_intervals_on_reference_partitioned: dict[
            int, dict[str, list[tuple[int, str, int, int, int, int]]]
        ] = generate_core_intervals_on_reference(
            padded_alignments_paths=padded_alignments_paths,
            index_consensusIDs=index_consensusIDs,
            consensus_json_index_partitioned=consensus_json_index_partitioned,
            consensus_json_paths=consensus_json_paths,
            threads=threads,
            batch_size=batchsize,
        )
        # now we need to find the intersecting trf intervals for each core interval on reference
        # Find TRF overlaps for all core intervals in parallel using bedtools
        trf_overlaps_on_reference_partitioned: dict[
            int, dict[str, list[tuple[str, int, int, int]]]
        ] = find_trf_overlaps_for_core_intervals(
            core_intervals_on_reference=core_intervals_on_reference_partitioned,
            input_trf=Path(input_trf),
            reference=Path(input_reference),
            threads=threads,
            tmp_dir=Path(tmp_dir),
        )

        consensusIDs_by_partition = get_consensusIDs_by_partition(
            index_alignments_partitioned=index_alignments_partitioned
        )

        svPatterns_paths: dict[int, Path] = {
            i: Path(tmp_dir) / f"svPatterns_partition_{i}.tsv" for i in range(threads)
        }

        # Prepare arguments for parallel processing
        log.info(
            f"Processing consensus objects to SV patterns using {threads} parallel workers..."
        )
        process_args = []
        for partition_index in range(threads):
            process_args.append(
                SVPatternProcessingParams(
                    svPatterns_path=svPatterns_paths[partition_index],
                    consensusIDs=list(consensusIDs_by_partition[partition_index]),
                    consensus_json_index_partition=consensus_json_index_partitioned[
                        partition_index
                    ],
                    consensus_json_paths=consensus_json_paths,
                    partition_index=partition_index,
                    padded_alignments_path=padded_alignments_paths[partition_index],
                    core_intervals_on_reference=core_intervals_on_reference_partitioned[
                        partition_index
                    ],
                    alignments_index_partition=index_alignments_partitioned[
                        partition_index
                    ],
                    trf_overlaps_on_reference_partition=trf_overlaps_on_reference_partitioned[
                        partition_index
                    ],
                    input_reference=Path(input_reference),
                    samplename=samplename,
                    min_bnd_size=min_bnd_size,
                    min_signal_size=min_signal_size,
                    max_del_size=max_del_size,
                    distance_scale=distance_scale,
                    falloff=falloff,
                    max_fourrelations_gap_size=max_fourrelations_gap_size,
                    signal_loss_log_dir=_signal_loss_log_dir
                    if _signal_loss_log_dir
                    else Path(tmp_dir),
                    logfile=_logfile,
                    batchsize=100,
                    log_level=log_level,
                )
            )

        # Process partitions in parallel
        with mp.Pool(processes=threads) as pool:
            list(
                tqdm(
                    pool.map(_process_consensus_objects_to_svPatterns, process_args),
                    total=threads,
                    desc="Processing consensus sequences to SV patterns",
                )
            )

        log.info(
            f"Completed parallel processing of SV patterns across {threads} partitions."
        )

        # Create a temporary file for concatenated JSON lines
        tmp_svpatterns_json = Path(tmp_dir) / "svpatterns_all.json"

        # open each svPatterns file sequentially, read lines and write to tmp file
        with open(tmp_svpatterns_json, "w") as outfile:
            for partition_idx in range(threads):
                svp_path = svPatterns_paths[partition_idx]
                if not svp_path.exists():
                    log.warning(
                        f"SV patterns file for partition {partition_idx} does not exist at {svp_path}"
                    )
                    continue
                with open(svp_path, "r") as infile:
                    shutil.copyfileobj(infile, outfile)

        log.info(f"Concatenated {threads} SV pattern partition files")

        # Create the database and populate it from the temporary JSON file
        log.info(f"Creating SV patterns database at {output_svPatterns}")
        SVpatterns.create_svPatterns_db(database=Path(output_svPatterns))
        SVpatterns.write_svPatterns_to_db(
            database=Path(output_svPatterns),
            svPatterns_file=tmp_svpatterns_json,
            timeout=30,
        )
        log.info(f"Successfully wrote SV patterns to database {output_svPatterns}")

        # Merge per-worker signal loss logs into the main signal loss log
        # if signal_loss_log is not None:
        #     _main_loss_log = Path(signal_loss_log)
        #     # Close the main process logger first so we don't have two
        #     # open handles to the same file when we append worker data.
        #     _main_logger = get_signal_loss_logger()
        #     _main_entries = _main_logger.get_entry_count()
        #     _main_logger.flush()
        #     _main_logger.close()

        #     _loss_log_dir = _signal_loss_log_dir if _signal_loss_log_dir else Path(tmp_dir)
        #     _worker_logs = sorted(_loss_log_dir.glob("signal_loss.partition_*.tsv"))
        #     if _worker_logs:
        #         with open(_main_loss_log, "a") as outfile:
        #             for wlog in _worker_logs:
        #                 with open(wlog, "r") as infile:
        #                     # Skip header line from each worker log
        #                     header = infile.readline()
        #                     shutil.copyfileobj(infile, outfile)
        #         log.info(
        #             f"Merged {len(_worker_logs)} per-worker signal loss logs into {_main_loss_log}"
        #         )
        #     log.info(
        #         f"Signal loss log written to {_main_loss_log} "
        #         f"({_main_entries} entries from main process)"
        #     )


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
        default=50,
        help="Minimum size of BND signal to consider (default: 50).",
    )
    parser.add_argument(
        "--max-del-size",
        type=int,
        default=500_000,
        help="Maximum size of deletions to consider for SVpattern generation (default: 500000). They are treated as translocations otherwise. Warning: larger deletions are often spurious and can lead to memory issues.",
    )
    parser.add_argument(
        "--max-fourrelations-gap-size",
        type=int,
        default=500_000,
        help="Maximum gap size for four-relations analysis in SVpattern generation (default: 500000). Controls detection of complex SV patterns like inversions.",
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
    parser.add_argument(
        "--signal-loss-log",
        type=Path,
        required=False,
        default=None,
        help="Path to write a detailed signal-loss log file (TSV). Records every signal, "
        "primitive, or pattern that was filtered, merged, skipped, or otherwise discarded "
        "during the pipeline. If not provided, loss events are logged via the Python logger only.",
    )
    parser.add_argument(
        "--logfile",
        type=Path,
        required=False,
        default=None,
        help="Path to log file. If not provided, logs will be printed to stdout.",
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

    # Add file handler if log file is specified
    if args.logfile:
        file_handler = logging.FileHandler(args.logfile)
        file_handler.setLevel(getattr(logging, args.log_level))
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

    svPatterns_from_consensus_sequences(
        samplename=args.samplename,
        input_consensus_container_results=args.consensus,
        input_reference=args.reference,
        input_trf=args.trf,
        output_consensus_fasta=args.output_consensus_fasta,
        output_consensus_to_reference_alignments=args.output_consensus_alignments,
        output_svPatterns=args.output_svpatterns,
        threads=args.threads,
        min_signal_size=args.min_signal_size,
        min_bnd_size=args.min_bnd_size,
        max_del_size=args.max_del_size,
        distance_scale=args.distance_scale,
        falloff=args.falloff,
        max_fourrelations_gap_size=args.max_fourrelations_gap_size,
        tmp_dir_path=args.tmp_dir_path,
        dont_merge_horizontally=args.dont_merge_horizontally,
        signal_loss_log=args.signal_loss_log,
        logfile=args.logfile,
        log_level=args.log_level,
    )

    # Database creation is now handled inside svPatterns_from_consensus_sequences
    return


if __name__ == "__main__":
    main()
