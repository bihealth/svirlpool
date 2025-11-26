import logging
import multiprocessing as mp
import time
from collections.abc import Iterator
from datetime import datetime

import attrs
import cattrs
from pysam import AlignedSegment
from tqdm import tqdm

from ..svcalling import genotyping
from ..util import datatypes
from ..util.minimizer_matching import MinimizerIndex
from ..util.util import Direction, get_read_position_on_ref
from . import consensus_class

log = logging.getLogger(__name__)


@attrs.define
class Adjacency:
    chr: str
    ref_pos: int
    gapsize: int
    tp_str: str  # needs to be of form t[p[ or t]p] or [p[t or ]p]t; p is chr:start, t is alt seq
    mate_id: str

    def unstructure(self):
        return cattrs.unstructure(self)


@attrs.define
class SVprimitive(datatypes.MergedSVSignal):  # can be ins,del,bndl,bndr
    samplename: str
    consensusID: str
    alignmentID: int  # ID of the alignment to the reference
    svID: int  # ID of the SV primitive on the consensus alignment
    aln_is_reverse: bool
    consensus_aln_interval: tuple[
        str, int, int
    ]  # chr, start, end of the according consensus sequence to reference alignment
    genotypeMeasurement: genotyping.GenotypeMeasurement | None = (
        None  # samplename: GenotypeMeasurement
    )
    adjacent_bnd: None | Adjacency = None
    similar_sequence_intervals_on_consensus: list[
        tuple[int, int]
    ] = []  # in core sequence coords

    @classmethod
    def from_merged_sv_signal(
        cls,
        merged_sv_signal: datatypes.MergedSVSignal,
        samplename: str,
        consensusID: str,
        alignmentID: int,
        svID: int,
        aln_is_reverse: bool,
        consensus_aln_interval: tuple[str, int, int],
    ):
        return cls(
            # Copy all fields from MergedSVSignal
            ref_start=merged_sv_signal.ref_start,
            ref_end=merged_sv_signal.ref_end,
            read_start=merged_sv_signal.read_start,
            read_end=merged_sv_signal.read_end,
            size=merged_sv_signal.size,
            sv_type=merged_sv_signal.sv_type,
            chr=merged_sv_signal.chr,
            repeatIDs=merged_sv_signal.repeatIDs,
            original_alt_sequences=merged_sv_signal.original_alt_sequences,
            original_ref_sequences=merged_sv_signal.original_ref_sequences,
            # Add new fields
            samplename=samplename,
            consensusID=consensusID,
            alignmentID=alignmentID,
            svID=svID,
            aln_is_reverse=aln_is_reverse,
            consensus_aln_interval=consensus_aln_interval,
        )

    def unstructure(self):
        return super().unstructure()

    def __lt__(self, other: object) -> bool:
        return super().__lt__(other)

    def __hash__(self) -> int:
        return super().__hash__()

    def get_vcfID(self) -> str:
        return f"{datatypes.SV_TYPE_DICT[self.sv_type]}.{self.consensusID}:{str(self.alignmentID)}.{str(self.svID)}"

    def get_vcfID_with_samplename(self) -> str:
        return f"{self.samplename}.{self.get_vcfID()}"

    def get_total_support(self) -> int:
        if self.genotypeMeasurement is None:
            return 0
        return max(
            len(self.genotypeMeasurement.supporting_reads_start),
            (
                len(self.genotypeMeasurement.supporting_reads_end)
                if self.genotypeMeasurement.supporting_reads_end
                else 0
            ),
        )

    # def get_total_coverage(self) -> int:
    #     if self.genotypeMeasurement is None:
    #         return 0
    #     return max(self.genotypeMeasurement.estimated_total_depth_start, (self.genotypeMeasurement.estimated_total_depth_end if self.genotypeMeasurement.estimated_total_depth_end else 0), 0)


def add_adjacencies_to_svPrimitives(
    svPrimitives: list[SVprimitive],
) -> list[SVprimitive]:
    # sort all svPrimitives by consensusID, chr, and read_start
    svPrimitives = sorted(svPrimitives, key=lambda x: (x.consensusID, x.read_start))

    # two BNDs are adjacent if they are consecutive in the list and have the same consensusID
    # given a pair (A,B) of adjacent svPs:
    i = 0
    while i < len(svPrimitives) - 1:
        A: SVprimitive = svPrimitives[i]
        B: SVprimitive = svPrimitives[i + 1]
        # if A.sv_type >= 3 and B.sv_type >= 3:
        #     print(f"Checking {A.get_vcfID()} and {B.get_vcfID()} for adjacency. A.consensusID == B.consensusID: {A.consensusID == B.consensusID}, A.alignmentID != B.alignmentID: {A.alignmentID != B.alignmentID}, A.sv_type >= 3 and B.sv_type >= 3: {A.sv_type >= 3 and B.sv_type >= 3}")
        if (
            A.consensusID == B.consensusID
            and A.alignmentID != B.alignmentID
            and A.sv_type >= 3
            and B.sv_type >= 3
        ):
            tp_str_A = ""
            tp_str_B = ""
            # if A and B are BNDs, add the adjacency to both
            if A.sv_type == 4:  # BNDR
                if B.sv_type == 3:  # t[p[, ]p]t
                    tp_str_A = "t[p["
                    tp_str_B = "]p]t"
                if B.sv_type == 4:  # t]p], t]p]
                    tp_str_A = "t]p]"
                    tp_str_B = "t]p]"
            if A.sv_type == 3:  # BNDL
                if B.sv_type == 3:  # [p[t, [p[t
                    tp_str_A = "[p[t"
                    tp_str_B = "[p[t"
                if B.sv_type == 4:  # ]p]t, t[p[
                    tp_str_A = "]p]t"
                    tp_str_B = "t[p["
            A.adjacent_bnd = Adjacency(
                chr=B.chr,
                ref_pos=B.ref_start,
                gapsize=abs(B.read_start - A.read_start - 1),
                tp_str=tp_str_A,
                mate_id=B.get_vcfID_with_samplename(),
            )
            B.adjacent_bnd = Adjacency(
                chr=A.chr,
                ref_pos=A.ref_start,
                gapsize=abs(B.read_start - A.read_start - 1),
                tp_str=tp_str_B,
                mate_id=A.get_vcfID_with_samplename(),
            )
            i += 2
        else:
            i += 1

    return svPrimitives


def add_genotypeMeasurements_to_SVprimitives(
    svps: list[SVprimitive],
    pysam_aln: AlignedSegment,
    intervals_cutread_alignments: list[tuple[int, int, str, bool]],
    core_interval_start: int,
) -> None:
    for svp in svps:
        start_on_consensus = get_read_position_on_ref(
            position=svp.read_start, alignment=pysam_aln, direction=Direction.LEFT
        )
        # estimated_total_depth_start = -1 # will be set in the multisample sv calling step
        # find all reads that support this SV signal. This means to find all overlapping intervals_cutread_alignments
        supporting_reads_start = []
        svp_start = svp.read_start - core_interval_start
        for start, end, readname, _forward in intervals_cutread_alignments:
            start, end = (end, start) if start > end else (start, end)
            if start <= svp_start <= end:
                supporting_reads_start.append(readname)
        if svp.sv_type == 1:  # not deletion
            svp.genotypeMeasurement = genotyping.GenotypeMeasurement(
                start_on_consensus=start_on_consensus,
                end_on_consensus=None,  # None for insertions and breakends
                # estimated_total_depth_start=estimated_total_depth_start,
                # estimated_total_depth_end=None,  # None for insertions and breakends
                supporting_reads_start=supporting_reads_start,
                supporting_reads_end=None,  # None for insertions and breakends
            )
            continue
        # in case of a deletion
        end_on_consensus = get_read_position_on_ref(
            position=svp.read_end, alignment=pysam_aln, direction=Direction.RIGHT
        )
        supporting_reads_end = []
        svp_end = svp.read_end - core_interval_start
        for start, end, readname, _forward in intervals_cutread_alignments:
            if start <= svp_end <= end:
                supporting_reads_end.append(readname)
        svp.genotypeMeasurement = genotyping.GenotypeMeasurement(
            start_on_consensus=start_on_consensus,
            end_on_consensus=end_on_consensus,  # None for insertions and breakends
            # estimated_total_depth_start=estimated_total_depth_start,
            # estimated_total_depth_end=estimated_total_depth_end,  # None for insertions and breakends
            supporting_reads_start=supporting_reads_start,
            supporting_reads_end=supporting_reads_end,  # None for insertions and breakends
        )


def calculate_similar_sequence_intervals_on_consensus(
    svPrimitive: SVprimitive, core_padding_left: int, minimizer_index: MinimizerIndex
) -> list[tuple[int, int]]:
    sv_start_on_core = svPrimitive.read_start - core_padding_left
    if svPrimitive.sv_type == 0:  # INS
        start = sv_start_on_core
        end = sv_start_on_core + svPrimitive.size
        if svPrimitive.size < minimizer_index.w + minimizer_index.k:
            diff = minimizer_index.w + minimizer_index.k - (end - start)
            start -= diff // 2
            end += diff // 2
            query = svPrimitive.get_alt_sequence()
            if len(query) < minimizer_index.k:
                start = max(0, svPrimitive.read_start - 50 - core_padding_left)
                end = svPrimitive.read_end + 50 - core_padding_left
                return [(start, end)]
            else:
                return minimizer_index.find_similar_regions_of_query(query=query)
    elif svPrimitive.sv_type == 1:  # DEL
        radius_size = 30
        query = get_query_from_reference(
            reference=minimizer_index.reference,
            start=sv_start_on_core - radius_size,
            end=sv_start_on_core + radius_size,
        )
        # find similar regions of the query in the minimizer index
        return minimizer_index.find_similar_regions_of_query(query=query)

    else:  # BND
        radius_size = 200
        query = get_query_from_reference(
            reference=minimizer_index.reference,
            start=sv_start_on_core - radius_size,
            end=sv_start_on_core + radius_size,
        )
        return minimizer_index.find_similar_regions_of_query(query=query)


def get_query_from_reference(reference: str, start: int, end: int) -> str:
    # this function tried to get the sequence from the reference by the given interval
    # however, if the start is negative or the end is larger than the reference length, the interval
    # is dajusted to keep its length
    if start > end:
        start, end = end, start  # swap if start is larger than end
    if start < 0:
        start = 0
        end = start + (end - start)
    if end > len(reference):
        end = len(reference)
        start = end - (end - start)
    return reference[start:end]


def generate_SVprimitives(
    crs_container_results_iter: Iterator[consensus_class.CrsContainerResult],
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    pysam_cache: dict[int, AlignedSegment],
    intervals_core: dict[str, tuple[int, int]],
    samplename: str,
) -> list[SVprimitive]:
    """
    Generate SV primitives from consensus results.

    Args:
        use_minimizer_optimization: If True, use the faster minimizer-based approach
    """
    print_performance_times: bool = False
    result = []

    for crs_container_result in tqdm(crs_container_results_iter):
        
        
def container_result_to_svprimitives(
        consensus:consensus_class.Consensus,
        consensus_alignments: list[consensus_class.ConsensusAlignment],
        ) -> list[SVprimitive]:
    
    results : list[SVprimitive] = []

    for i, consensusAlignment in enumerate(consensus_alignments):
        svp_cache: list[SVprimitive] = []
        # trace back the core interval of the consensus sequence on the consensus to reference alignment
        pysam_alignment = pysam_cache[consensusAlignment.uid]
        ref_start = get_read_position_on_ref(
            alignment=pysam_alignment,
            position=intervals_core[consensus.ID][0],
            direction=Direction.LEFT,
        )
        ref_end = get_read_position_on_ref(
            alignment=pysam_alignment,
            position=intervals_core[consensus.ID][1],
            direction=Direction.RIGHT,
        )
        ref_start, ref_end = (
            (ref_end, ref_start)
            if ref_start > ref_end
            else (ref_start, ref_end)
        )  # ensure ref_start <= ref_end
        for j, svCandidate in enumerate(consensusAlignment.proto_svs):
            consensus_aln_interval = (
                pysam_alignment.reference_name,
                ref_start,
                ref_end,
            )
            svp_cache.append(
                SVprimitive.from_merged_sv_signal(
                    merged_sv_signal=svCandidate,
                    samplename=samplename,
                    consensusID=consensus.ID,
                    alignmentID=i,
                    svID=j,
                    aln_is_reverse=consensusAlignment.alignment.is_reverse(),
                    consensus_aln_interval=consensus_aln_interval,
                )
            )
        add_genotypeMeasurements_to_SVprimitives(
            svps=svp_cache,
            pysam_aln=pysam_alignment,
            intervals_cutread_alignments=consensus.intervals_cutread_alignments,
            core_interval_start=consensus.consensus_padding.padding_size_left,
        )

        # filter svp_cache for all svps that have no supporting reads
        svp_cache = [
            svp
            for svp in svp_cache
            if len(svp.genotypeMeasurement.supporting_reads_start) > 0
        ]  # this should really filter all empty svPrimitives, but they still occur downstream. why?

        for idx, svPrimitive in enumerate(svp_cache):
            debug_performance__len_consensus_seuquence = len(
                consensus.consensus_sequence
            )
            debug_performance__sum_svPrimitive_sizes = sum(
                (
                    svp.size
                    if svp.sv_type == 0
                    else (50 if svp.sv_type == 1 else 200)
                )
                for svp in svp_cache
            )
            debug_performance__num_svPrimitives = len(svp_cache)

            # Start timing
            start_time = time.perf_counter()

            # Use optimized or original method based on flag
            svPrimitive.similar_sequence_intervals_on_consensus = calculate_similar_sequence_intervals_on_consensus(
                svPrimitive=svPrimitive,
                core_padding_left=consensus.consensus_padding.padding_size_left,
                minimizer_index=minimizer_index,
            )

            if (
                not svPrimitive.similar_sequence_intervals_on_consensus
                or len(svPrimitive.similar_sequence_intervals_on_consensus) == 0
            ):
                # make the location of the SV primitive on the consensus sequence
                start = (
                    svPrimitive.read_start
                    - 50
                    - consensus.consensus_padding.padding_size_left
                )
                end = (
                    svPrimitive.read_end
                    + 50
                    - consensus.consensus_padding.padding_size_left
                )
                svPrimitive.similar_sequence_intervals_on_consensus = [
                    (start, end)
                ]

            # End timing and log performance data
            end_time = time.perf_counter()
            execution_time_ms = (end_time - start_time) * 1000.0

            # Get performance logger
            perf_logger = logging.getLogger("performance")
            perf_logger.info(
                f"{consensus.ID},{idx},{debug_performance__len_consensus_seuquence},{debug_performance__num_svPrimitives},{debug_performance__sum_svPrimitive_sizes},{svPrimitive.sv_type},{svPrimitive.size},{execution_time_ms:.3f}"
            )
            results.append(svPrimitive)
    return results



def _process_consensus_batch_serialized(
    serialized_batch_data: dict,
) -> list[SVprimitive]:
    """Process a batch of serialized consensus data to generate SVprimitives.

    Args:
        serialized_batch_data: Dictionary containing serialized batch data

    Returns:
        List of SVprimitives
    """
    from .consensus_align import (
        deserialize_pysam_AlignedSegment,
    )  # avoid circular imports

    # Deserialize the batch data
    consensus_batch = []
    for serialized_consensus in serialized_batch_data["consensus_batch"]:
        consensus = cattrs.structure(serialized_consensus, consensus_class.Consensus)
        consensus_batch.append(consensus)

    dict_alignments_subset = {}
    for consensus_id, serialized_alignments in serialized_batch_data[
        "dict_alignments"
    ].items():
        dict_alignments_subset[consensus_id] = []
        for serialized_alignment in serialized_alignments:
            consensus_alignment = cattrs.structure(
                serialized_alignment, consensus_class.ConsensusAlignment
            )
            dict_alignments_subset[consensus_id].append(consensus_alignment)

    intervals_core_subset = serialized_batch_data["intervals_core"]
    samplename = serialized_batch_data["samplename"]

    # Deserialize pysam alignments
    pysam_cache_subset = {}
    for uid_str, serialized_alignment_data in serialized_batch_data[
        "pysam_cache"
    ].items():
        uid = int(uid_str)
        pysam_cache_subset[uid] = deserialize_pysam_AlignedSegment(
            serialized_alignment_data
        )

    # Process using the existing single-process logic
    return process_consensus_batch(
        consensus_batch=consensus_batch,
        dict_alignments_subset=dict_alignments_subset,
        alignment_cache_subset={},  # Not needed since we have pysam_cache_subset
        intervals_core_subset=intervals_core_subset,
        samplename=samplename,
        pysam_cache_subset=pysam_cache_subset,
    )


def _create_serialized_consensus_batch(
    consensus_batch: list,
    dict_alignments: dict,
    pysam_cache: dict,
    intervals_core: dict,
    samplename: str,
) -> dict:
    """Create a serialized batch for the given consensus sequences.

    Args:
        consensus_batch: List of consensus objects for this batch
        dict_alignments: Full dictionary of alignments
        pysam_cache: Full pysam cache
        intervals_core: Full dictionary of core intervals
        samplename: Sample name

    Returns:
        Serialized batch data ready for processing
    """
    from .consensus_align import (
        serialize_pysam_AlignedSegment,
    )  # avoid circular imports

    # Get consensus IDs for this batch
    consensus_ids = {consensus.ID for consensus in consensus_batch}

    # Create sub-containers for this batch
    batch_dict_alignments = {}
    batch_intervals_core = {}
    batch_pysam_cache = {}

    # Collect all UIDs needed for this batch
    required_uids = set()

    for consensus_id in consensus_ids:
        if consensus_id in dict_alignments:
            batch_dict_alignments[consensus_id] = dict_alignments[consensus_id]
            for consensus_alignment in dict_alignments[consensus_id]:
                required_uids.add(consensus_alignment.uid)

        if consensus_id in intervals_core:
            batch_intervals_core[consensus_id] = intervals_core[consensus_id]

    # Create pysam cache subset for this batch
    for uid in required_uids:
        if uid in pysam_cache:
            batch_pysam_cache[uid] = serialize_pysam_AlignedSegment(pysam_cache[uid])

    # Serialize the batch data
    serialized_batch = {
        "consensus_batch": [consensus.unstructure() for consensus in consensus_batch],
        "dict_alignments": {
            consensus_id: [
                consensus_alignment.unstructure() for consensus_alignment in alignments
            ]
            for consensus_id, alignments in batch_dict_alignments.items()
        },
        "intervals_core": batch_intervals_core,
        "pysam_cache": {
            str(uid): serialized_alignment
            for uid, serialized_alignment in batch_pysam_cache.items()
        },
        "samplename": samplename,
    }

    return serialized_batch


def process_consensus_batch(
    consensus_batch: list,
    dict_alignments_subset: dict[str, list],
    alignment_cache_subset: dict[int, datatypes.Alignment],
    intervals_core_subset: dict[str, tuple[int, int]],
    samplename: str,
    pysam_cache_subset: dict[int, AlignedSegment] = None,
) -> list[SVprimitive]:
    """Process a batch of consensus sequences in a single worker process."""
    # Use provided pysam cache or convert from alignment cache
    if pysam_cache_subset is None:
        pysam_cache_subset = {
            alignment_hash: alignment.to_pysam()
            for alignment_hash, alignment in alignment_cache_subset.items()
        }

    print_performance_times: bool = False
    result = []

    for consensus in consensus_batch:
        cache: list[SVprimitive] = []
        if print_performance_times:
            start_idx_build_time = datetime.now()
            minimizer_index = MinimizerIndex(consensus.consensus_sequence, k=15, w=10)
            end_idx_build_time = datetime.now()
            log.info(
                f"Minimizer index built in {(end_idx_build_time - start_idx_build_time).total_seconds() * 1000:.3f} ms"
            )
        else:
            minimizer_index = MinimizerIndex(consensus.consensus_sequence, k=15, w=10)

        start_querying_time = datetime.now()
        for i, consensusAlignment in enumerate(
            dict_alignments_subset.get(consensus.ID, [])
        ):
            svp_cache: list[SVprimitive] = []
            # trace back the core interval of the consensus sequence on the consensus to reference alignment
            pysam_alignment = pysam_cache_subset[consensusAlignment.uid]
            ref_start = get_read_position_on_ref(
                alignment=pysam_alignment,
                position=intervals_core_subset[consensus.ID][0],
                direction=Direction.LEFT,
            )
            ref_end = get_read_position_on_ref(
                alignment=pysam_alignment,
                position=intervals_core_subset[consensus.ID][1],
                direction=Direction.RIGHT,
            )
            ref_start, ref_end = (
                (ref_end, ref_start) if ref_start > ref_end else (ref_start, ref_end)
            )  # ensure ref_start <= ref_end

            for j, svCandidate in enumerate(consensusAlignment.proto_svs):
                consensus_aln_interval = (
                    pysam_alignment.reference_name,
                    ref_start,
                    ref_end,
                )
                svp_cache.append(
                    SVprimitive.from_merged_sv_signal(
                        merged_sv_signal=svCandidate,
                        samplename=samplename,
                        consensusID=consensus.ID,
                        alignmentID=i,
                        svID=j,
                        aln_is_reverse=consensusAlignment.alignment.is_reverse(),
                        consensus_aln_interval=consensus_aln_interval,
                    )
                )

            add_genotypeMeasurements_to_SVprimitives(
                svps=svp_cache,
                pysam_aln=pysam_alignment,
                intervals_cutread_alignments=consensus.intervals_cutread_alignments,
                core_interval_start=consensus.consensus_padding.padding_size_left,
            )

            # filter svp_cache for all svps that have no supporting reads
            svp_cache = [
                svp
                for svp in svp_cache
                if len(svp.genotypeMeasurement.supporting_reads_start) > 0
            ]

            for idx, svPrimitive in enumerate(svp_cache):
                debug_performance__len_consensus_seuquence = len(
                    consensus.consensus_sequence
                )
                debug_performance__sum_svPrimitive_sizes = sum(
                    (
                        svp.size
                        if svp.sv_type == 0
                        else (50 if svp.sv_type == 1 else 200)
                    )
                    for svp in svp_cache
                )
                debug_performance__num_svPrimitives = len(svp_cache)

                # Start timing
                start_time = time.perf_counter()

                svPrimitive.similar_sequence_intervals_on_consensus = (
                    calculate_similar_sequence_intervals_on_consensus(
                        svPrimitive=svPrimitive,
                        core_padding_left=consensus.consensus_padding.padding_size_left,
                        minimizer_index=minimizer_index,
                    )
                )

                if (
                    not svPrimitive.similar_sequence_intervals_on_consensus
                    or len(svPrimitive.similar_sequence_intervals_on_consensus) == 0
                ):
                    # make the location of the SV primitive on the consensus sequence
                    start = (
                        svPrimitive.read_start
                        - 50
                        - consensus.consensus_padding.padding_size_left
                    )
                    end = (
                        svPrimitive.read_end
                        + 50
                        - consensus.consensus_padding.padding_size_left
                    )
                    svPrimitive.similar_sequence_intervals_on_consensus = [(start, end)]

                # End timing and log performance data
                end_time = time.perf_counter()
                execution_time_ms = (end_time - start_time) * 1000.0

                # Get performance logger
                perf_logger = logging.getLogger("performance")
                perf_logger.info(
                    f"{consensus.ID},{idx},{debug_performance__len_consensus_seuquence},{debug_performance__num_svPrimitives},{debug_performance__sum_svPrimitive_sizes},{svPrimitive.sv_type},{svPrimitive.size},{execution_time_ms:.3f}"
                )

            cache.extend(svp_cache)
        end_querying_time = datetime.now()
        if print_performance_times:
            log.info(
                f"Querying similar sequence intervals took {(end_querying_time - start_querying_time).total_seconds() * 1000:.3f} ms"
            )

    # add adjacencies to the cached SVprimitives
    result.extend(add_adjacencies_to_svPrimitives(cache))
    return result


def generate_SVprimitives_parallel(
    crs_container_results_iter: Iterator[consensus_class.CrsContainerResult],
    dict_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    pysam_cache: dict[int, AlignedSegment],
    intervals_core: dict[str, tuple[int, int]],
    samplename: str,
    n_workers: int = None,
    batch_size: int = 50,
) -> list[SVprimitive]:
    """
    Parallelized version of generate_SVprimitives with on-demand batch generation.
    """
    if n_workers < 2:
        return generate_SVprimitives(
            crs_container_results_iter=crs_container_results_iter,
            dict_alignments=dict_alignments,
            pysam_cache=pysam_cache,
            intervals_core=intervals_core,
            samplename=samplename,
        )
    if n_workers is None:
        n_workers = min(mp.cpu_count() - 1, 8)  # Leave one core free, cap at 8

    # Collect all consensus sequences first
    all_consensus = []
    for crs_container_result in crs_container_results_iter:
        for consensus in crs_container_result.consensus_dicts.values():
            all_consensus.append(consensus)

    if not all_consensus:
        return []

    # Create batches of consensus sequences
    consensus_batches = [
        all_consensus[i : i + batch_size]
        for i in range(0, len(all_consensus), batch_size)
    ]

    log.info(
        f"Processing {len(all_consensus)} consensus sequences in {len(consensus_batches)} batches with {n_workers} workers..."
    )

    # Create a generator function that yields serialized batches on-demand
    def batch_generator():
        for consensus_batch in consensus_batches:
            yield _create_serialized_consensus_batch(
                consensus_batch=consensus_batch,
                dict_alignments=dict_alignments,
                pysam_cache=pysam_cache,
                intervals_core=intervals_core,
                samplename=samplename,
            )

    # Process batches in parallel using the generator
    with mp.Pool(processes=n_workers) as pool:
        batch_results = list(
            tqdm(
                pool.imap(_process_consensus_batch_serialized, batch_generator()),
                total=len(consensus_batches),
                desc="Processing consensus batches",
            )
        )

    # Combine results
    result = []
    for batch_result in batch_results:
        result.extend(batch_result)

    return result
