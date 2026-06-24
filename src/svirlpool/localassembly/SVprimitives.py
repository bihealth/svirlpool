import logging

import attrs
import cattrs
from pysam import AlignedSegment

from ..svcalling import genotyping
from ..util import datatypes
from ..util.signal_loss_logger import get_signal_loss_logger
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

    def _get_description(self) -> str:
        # get description from superclass
        str_super = super()._get_description(
            self.chr
        )  # defines the region, type, size, read_span

        return f"{str_super}, samplename={self.samplename}, consensusID={self.consensusID}, crID={self.consensusID.split('.')[0]}, alignmentID={self.alignmentID}, svID={self.svID}, aln_is_reverse={self.aln_is_reverse}, consensus_aln_interval={self.consensus_aln_interval}, genotypeMeasurement={self.genotypeMeasurement}"

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
        # intervals_cutread_alignments stores core-relative coordinates (0 to core_size).
        # svp.read_start is in padded-consensus coordinates, so we subtract core_interval_start.
        supporting_reads_start = []
        svp_start = svp.read_start - core_interval_start
        for start, end, readname, _forward in intervals_cutread_alignments:
            start, end = (end, start) if start > end else (start, end)
            if start <= svp_start <= end:
                supporting_reads_start.append(readname)
        _svp_crID = (
            svp.consensusID.split(".")[0] if "." in svp.consensusID else svp.consensusID
        )
        log.debug(
            f"TRANSFORMED::add_genotypeMeasurements_to_SVprimitives::GENOTYPE_MEASUREMENT_START\tconsensusID={svp.consensusID}\tcrID={_svp_crID}\tsvID={svp.svID}\tsv_type={svp.sv_type}\t"
            f"svp_start_core={svp_start}\tcore_interval_start={core_interval_start}\t"
            f"n_supporting_reads_start={len(supporting_reads_start)}"
        )
        if svp.sv_type == 1:  # not deletion
            svp.genotypeMeasurement = genotyping.GenotypeMeasurement(
                start_on_consensus=start_on_consensus,
                end_on_consensus=None,  # None for insertions and breakends
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
        log.debug(
            f"TRANSFORMED::add_genotypeMeasurements_to_SVprimitives::GENOTYPE_MEASUREMENT_END\tconsensusID={svp.consensusID}\tcrID={_svp_crID}\tsvID={svp.svID}\tsv_type={svp.sv_type}\t"
            f"svp_end_core={svp_end}\tcore_interval_start={core_interval_start}\t"
            f"n_supporting_reads_end={len(supporting_reads_end)}"
        )
        svp.genotypeMeasurement = genotyping.GenotypeMeasurement(
            start_on_consensus=start_on_consensus,
            end_on_consensus=end_on_consensus,  # None for insertions and breakends
            # estimated_total_depth_start=estimated_total_depth_start,
            # estimated_total_depth_end=estimated_total_depth_end,  # None for insertions and breakends
            supporting_reads_start=supporting_reads_start,
            supporting_reads_end=supporting_reads_end,  # None for insertions and breakends
        )


def generate_SVprimitives(
    samplename: str,
    mergedSVs: list[datatypes.MergedSVSignal],
    consensus: consensus_class.Consensus,
    consensus_alignment: datatypes.Alignment,
    core_interval: tuple[str, int, int],  # chr, ref start, ref end
    alignmentID: int,  # ID following alignment order
) -> list[SVprimitive]:
    """Generates SVprimitives for a given consensus sequence and its alignments."""
    if consensus.consensus_padding is None:
        raise ValueError("Consensus padding is None, cannot generate SVprimitives.")
    svp_cache: list[SVprimitive] = []
    for j, svCandidate in enumerate(mergedSVs):
        consensus_aln_interval = (
            str(consensus_alignment.reference_name),
            core_interval[1],
            core_interval[2],
        )
        svp = SVprimitive.from_merged_sv_signal(
            merged_sv_signal=svCandidate,
            samplename=samplename,
            consensusID=consensus.ID,
            alignmentID=alignmentID,
            svID=j,
            aln_is_reverse=consensus_alignment.is_reverse(),
            consensus_aln_interval=consensus_aln_interval,
        )
        svp_cache.append(svp)
    add_genotypeMeasurements_to_SVprimitives(
        svps=svp_cache,
        pysam_aln=consensus_alignment.to_pysam(),
        intervals_cutread_alignments=consensus.intervals_cutread_alignments,
        core_interval_start=consensus.consensus_padding.padding_size_left,
    )
    # filter svp_cache for all svps that have no supporting reads
    _all_svp_count = len(svp_cache)
    _filtered_no_support = [
        svp
        for svp in svp_cache
        if svp.genotypeMeasurement is None
        or len(svp.genotypeMeasurement.supporting_reads_start) == 0
    ]
    svp_cache = [
        svp
        for svp in svp_cache
        if svp.genotypeMeasurement is not None
        and len(svp.genotypeMeasurement.supporting_reads_start) > 0
    ]
    if _filtered_no_support:
        _loss_logger = get_signal_loss_logger()
        for _fsvp in _filtered_no_support:
            _loss_logger.log_no_support(
                stage="generate_SVprimitives",
                consensusID=consensus.ID,
                reason="no_supporting_reads",
                details={
                    "sv_type": _fsvp.sv_type,
                    "size": _fsvp.size,
                    "ref_start": _fsvp.ref_start,
                    "ref_end": _fsvp.ref_end,
                    "chr": _fsvp.chr,
                    "alignmentID": _fsvp.alignmentID,
                    "svID": _fsvp.svID,
                },
            )
        log.info(
            f"generate_SVprimitives: consensus {consensus.ID}: "
            f"filtered {len(_filtered_no_support)}/{_all_svp_count} SVprimitives with no supporting reads"
        )
    return svp_cache
