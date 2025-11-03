import typing

import attrs
import cattrs
import pysam
from intervaltree import Interval

# is part of a ReadAlignmentFragment object
# types of SVs


@attrs.define
class SVsignal:
    ref_start: int
    ref_end: int
    read_start: int
    read_end: int
    size: int
    sv_type: int  # 0: INS, 1:DEL, 3:BNDL, 4:BNDR # 2 is DEL as well.

    def unstructure(self):
        return cattrs.unstructure(self)

    def __eq__(self, __value: object) -> bool:
        return (
            self.ref_start == __value.ref_start
            and self.ref_end == __value.ref_end
            and self.read_start == __value.read_start
            and self.read_end == __value.read_end
            and self.size == __value.size
            and self.sv_type == __value.sv_type
        )

    def __lt__(self, __value: object) -> bool:
        if self.ref_start > __value.ref_start:
            return False
        if self.ref_start < __value.ref_start:
            return True
        else:
            if self.ref_end > __value.ref_end:
                return False
            if self.ref_end < __value.ref_end:
                return True
        return False


@attrs.define
class ExtendedSVsignal(SVsignal):
    chr: str
    chrID: int
    coverage: int
    readname: str
    sampleID: int
    forward: int
    repeatID: int = -1
    strength: float = 0.0

    def aug_readname(self):
        return f"{self.readname}.{self.sampleID}"

    def unstructure(self):
        return super().unstructure()

    def __eq__(self, __value: object) -> bool:
        return (
            super().__eq__(__value)
            and self.chr == __value.chr
            and self.chrID == __value.chrID
            and self.coverage == __value.coverage
            and self.readname == __value.readname
            and self.sampleID == __value.sampleID
            and self.forward == __value.forward
            and self.repeatID == __value.repeatID
        )


@attrs.define
class ReadAlignmentSignals:
    sampleID: int
    read_name: str
    reference_name: str
    alignment_forward: bool
    SV_signals: typing.List[SVsignal] = []

    def unstructure(self):
        return cattrs.unstructure(self)


# not inherited from ReadAlignmentSignals, because its object of change
@attrs.define
class ReadAlignmentFragment:
    sampleID: int
    read_name: str
    reference_name: str
    referenceID: int
    reference_alignment_start: int
    reference_alignment_end: int
    read_alignment_start: int
    read_alignment_end: int
    mapping_quality: int
    alignment_forward: bool
    alignment_secondary: bool
    alignment_supplementary: bool
    effective_interval: tuple[str, int, int] | None = (
        None  # the span on the reference without loose ends, that are not flaked by unique region
    )
    SV_signals: typing.List[SVsignal] = []

    def unstructure(self):
        return cattrs.unstructure(self)

    def __eq__(self, __value: object) -> bool:
        # check if any SV signal is not equal
        if len(self.SV_signals) != len(__value.SV_signals):
            return False
        if not all([s in __value.SV_signals for s in self.SV_signals]):
            return False
        return (
            self.sampleID == __value.sampleID
            and self.read_name == __value.read_name
            and self.reference_name == __value.reference_name
            and self.referenceID == __value.referenceID
            and self.reference_alignment_start == __value.reference_alignment_start
            and self.reference_alignment_end == __value.reference_alignment_end
            and self.read_alignment_start == __value.read_alignment_start
            and self.read_alignment_end == __value.read_alignment_end
            and self.mapping_quality == __value.mapping_quality
            and self.alignment_forward == __value.alignment_forward
            and self.alignment_secondary == __value.alignment_secondary
            and self.alignment_supplementary == __value.alignment_supplementary
            and self.effective_interval == __value.effective_interval
        )


@attrs.define
class CrSeed:
    chr: str
    chrID: int
    start: int
    end: int
    readIDs: list
    min_extents: tuple[int | None, int | None] = (0, 0)  # on ref coords
    max_extents: tuple[int | None, int | None] = (0, 0)  # on ref coords
    actual_extents: tuple[int | None, int | None] = (0, 0)  # on ref coords
    # 0      1          2        3           4         5       6       7      8         9         10
    # chrID, ref_start, ref_end, read_start, read_end, svtype, svsize, depth, repeatID, sampleID, readID
    signals: list[ExtendedSVsignal] = []
    valid_intervals: list[Interval] = []

    def __hash__(self) -> int:
        return hash((self.chrID, self.start, self.end, tuple(self.readIDs)))

    def __lt__(self, __o: object) -> bool:
        if self.chrID > __o.chrID:
            return False
        else:
            if self.start > __o.start:
                return False
            if self.start < __o.start:
                return True
            else:
                if self.end > __o.end:
                    return False
                if self.end < __o.end:
                    return True
        return False

    def unstructure(self):
        return cattrs.unstructure(self)


# regions can get info about coverage, avg aln scores etc.
# to allow for better filtering later on
@attrs.define
class CandidateRegion:
    crID: int
    chr: str
    referenceID: int
    referenceStart: int
    referenceEnd: int
    sv_signals: list[ExtendedSVsignal]

    def __lt__(self, __o: object) -> bool:
        if self.referenceID > __o.referenceID:
            return False
        if self.referenceID < __o.referenceID:
            return True
        else:
            if self.referenceStart > __o.referenceStart:
                return False
            if self.referenceStart < __o.referenceStart:
                return True
            else:
                if self.referenceEnd > __o.referenceEnd:
                    return False
                if self.referenceEnd < __o.referenceEnd:
                    return True
        return False

    def __hash__(self) -> int:
        return self.crID

    def unstructure(self):
        return cattrs.unstructure(self)

    def region_string(self):
        return f"{self.chr}:{self.referenceStart}-{self.referenceEnd}"


@attrs.define
class SequenceObject:
    name: str
    id: str
    sequence: str
    description: str | None = None
    qualities: list[int] | None = None

    def unstructure(self):
        return cattrs.unstructure(self)


@attrs.define
class Alignment:
    readname: str
    sampleID: int
    reference_name: str
    reference_ID: int
    samdict: dict
    headerdict: dict

    def to_pysam(self):
        header = pysam.AlignmentHeader.from_dict(self.headerdict)
        return pysam.AlignedSegment.from_dict(self.samdict, header)

    def aug_name(self):
        return f"{self.readname}.{self.sampleID}"

    def unstructure(self):
        cattrs.unstructure(self)


# @attrs.define
# class Alignment:
#     readname:str
#     sampleID:int
#     flag:int
#     reference_name:str
#     reference_ID:int
#     reference_start:str
#     mapq:int
#     cigar:str
#     is_reverse:bool
#     is_supplementary:bool
#     is_secondary:bool
#     sequence:str|None
#     qualities:str|None
#     pysamHeader:dict|None
#     def aug_name(self):
#         return f"{self.readname}.{self.sampleID}"
#     def unstructure(self):
#         cattrs.unstructure(self)


@attrs.define
class ConsensusInfo:
    crID: int
    clusterID: int
    subcluster: int
    n_used_reads: int

    def unstructure(self):
        return cattrs.unstructure(self)


@attrs.define
class Consensus:
    ID: str
    crIDs: list[int]
    info: list[ConsensusInfo]
    consensus_sequence: str
    chosen_reads: list[SequenceObject]
    cut_read_alignments: list[Alignment]  # includes cut read sequences
    cut_read_intervals: dict[int]
    initial_read_alignments: list[Alignment]  # excludes read sequences
    # unused_reads:list[SequenceObject]
    consensus_to_initial_reference_alignments: list[Alignment] | None = None

    def unstructure(self):
        return cattrs.unstructure(self)


@attrs.define
class CrsContainerResult:
    consensus_dicts: dict[str, Consensus]  # consenus ID to consensus object
    unused_reads: list[SequenceObject]  # read ID to read object

    def unstructure(self):
        return cattrs.unstructure(self)


@attrs.define
class InDelBndCandidate:
    svType: str  # INS DEL BNDL BNDR
    svSize: int
    consensus_name: str
    refChrID: int
    refChr: str  # reference chromosome
    refStart: int  # start position of SV in reference
    refEnd: int  # end position of SV in reference
    consensus_start: int  # start ppsition of SV in consensus
    consensus_end: int  # end position of SV in consensus
    consensus_length: int  # length of consensus sequence
    depthSupport: dict  # of form {sampleID:[depth_l,depth_r,error_l,error_r,support]}
    repIDs: list  # list of IDs of repeats
    is_reverse: int

    def unstructure(self):
        r = cattrs.unstructure(self)
        r["repIDs"] = list(map(int, r["repIDs"]))
        return r

    def __eq__(self, __o: object) -> bool:
        return (
            self.svType == __o.svType
            and self.svSize == __o.svSize
            and self.consensus_name == __o.consensus_name
            and self.refChrID == __o.refChrID
            and self.refStart == __o.refStart
            and self.refEnd == __o.refEnd
            and self.consensus_start == __o.consensus_start
            and self.consensus_end == __o.consensus_end
            and self.depthSupport == __o.depthSupport
            and self.repIDs == __o.repIDs
            and self.is_reverse == __o.is_reverse
        )


@attrs.define
class IntraSVCandidate:
    svType: str  # INS DEL BND
    svSize: int
    consensus_name: str
    referenceIntervals: list
    referenceStrings: list
    consensusIntervals: list
    consensusStrings: list
    depthSupport: dict  # of form {sampleID:[depth_l,depth_r,error_l,error_r,support]}
    repIDs: typing.List[int]  # list of IDs of repeats

    def unstructure(self):
        r = cattrs.unstructure(self)
        r["repIDs"] = list(map(int, r["repIDs"]))
        return r

    def __eq__(self, __o: object) -> bool:
        return (
            self.svType == __o.svType
            and self.svSize == __o.svSize
            and self.consensus_name == __o.consensus_name
            and self.referenceIntervals == __o.referenceIntervals
            and self.consensusIntervals == __o.consensusIntervals
            and self.depthSupport == __o.depthSupport
            and self.repIDs == __o.repIDs
        )
