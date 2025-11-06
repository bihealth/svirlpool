import attrs
import cattrs
import pysam
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval

# is part of a ReadAlignmentFragment object
# types of SVs

# @attrs.define
# class SVtypes(Enum):
#     INS = 0
#     DEL = 1
#     # 2 is blocked for DELR
#     BNDL = 3
#     BNDR = 4
#     INV = 5
#     DUP = 6
#     def unstructure(self):
#         return cattrs.unstructure(self)

SV_TYPE_DICT = {0: "INS", 1: "DEL", 2: "DEL", 3: "BND", 4: "BND", 5: "INV", 6: "DUP"}


@attrs.define
class SVsignal:
    ref_start: int
    ref_end: int
    read_start: int
    read_end: int
    size: int
    sv_type: int  # can be any of the SVtypes

    # reverse: bool = False
    def unstructure(self):
        return cattrs.unstructure(self)

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

    def __eq__(self, __value: object) -> bool:
        return (
            self.ref_start == __value.ref_start
            and self.ref_end == __value.ref_end
            and self.read_start == __value.read_start
            and self.read_end == __value.read_end
            and self.size == __value.size
            and self.sv_type == __value.sv_type
        )

    def __hash__(self) -> int:
        return hash((
            self.ref_start,
            self.ref_end,
            self.read_start,
            self.read_end,
            self.size,
            self.sv_type,
        ))


@attrs.define
class ExtendedSVsignal(SVsignal):
    chr: str
    chrID: int
    coverage: int
    readname: str
    samplename: str
    forward: int
    repeatID: int = -1
    strength: float = 0.0

    # def aug_readname(self):
    #     return f"{self.readname}.{self.samplename}"
    def unstructure(self):
        return super().unstructure()

    def __eq__(self, __value: object) -> bool:
        return (
            super().__eq__(__value)
            and self.chr == __value.chr
            and self.chrID == __value.chrID
            and self.coverage == __value.coverage
            and self.readname == __value.readname
            and self.samplename == __value.samplename
            and self.forward == __value.forward
            and self.repeatID == __value.repeatID
        )


_converter = cattrs.Converter()


# not inherited from ReadAlignmentSignals, because its object of change
@attrs.define
class ReadAlignmentFragment:
    samplename: str
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
    SV_signals: list[SVsignal] = []

    def unstructure(self):
        return cattrs.unstructure(self)

    @classmethod
    def from_unstructured(cls, data: dict):
        """Create ReadAlignmentFragment from unstructured data."""
        return _converter.structure(data, cls)

    def __eq__(self, __value: object) -> bool:
        # check if any SV signal is not equal
        if len(self.SV_signals) != len(__value.SV_signals):
            return False
        if not all(s in __value.SV_signals for s in self.SV_signals):
            return False
        return (
            self.samplename == __value.samplename
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


@attrs.define
class ReadAlignmentSignals:
    samplename: str
    read_name: str
    reference_name: str
    alignment_forward: bool
    SV_signals: list[SVsignal] = []

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

    @classmethod
    def from_sv_signals(
        cls, crID: int, sv_signals: list[ExtendedSVsignal], buffer: int = 300
    ):
        if len(sv_signals) == 0:
            raise ValueError("Cannot create CandidateRegion from empty sv_signals list")
        chr = sv_signals[0].chr
        referenceID = sv_signals[0].chrID
        referenceStart = max(0, min([s.ref_start for s in sv_signals]) - buffer)
        referenceEnd = max([s.ref_end for s in sv_signals]) + buffer
        return cls(
            crID=crID,
            chr=chr,
            referenceID=referenceID,
            referenceStart=referenceStart,
            referenceEnd=referenceEnd,
            sv_signals=sv_signals,
        )

    def __hash__(self) -> int:
        return self.crID

    def unstructure(self):
        return cattrs.unstructure(self)

    def region_string(self):
        return f"{self.chr}:{self.referenceStart}-{self.referenceEnd}"

    def get_read_names(self) -> set[str]:
        return {s.readname for s in self.sv_signals}

    def get_ReadAlignmentSignals(self) -> list[ReadAlignmentSignals]:
        # construct a ReadAlignmentSignals object from the sv_signals list
        # sort signals by readname and position
        signals = sorted(self.sv_signals, key=lambda x: (x.readname, x.ref_start))
        # iterate and for each subset of signals with the same readname, construct a ReadAlignmentSignals object
        read_alignment_signals = []
        current_readname = signals[0].readname
        current_samplename = signals[0].samplename
        current_reference_name = signals[0].chr
        current_alignment_forward = signals[0].forward
        current_SV_signals = []
        for signal in signals:
            if signal.readname == current_readname:
                current_SV_signals.append(signal)
            else:
                read_alignment_signals.append(
                    ReadAlignmentSignals(
                        read_name=current_readname,
                        samplename=current_samplename,
                        reference_name=current_reference_name,
                        alignment_forward=current_alignment_forward,
                        SV_signals=current_SV_signals,
                    )
                )
                current_readname = signal.readname
                current_samplename = signal.samplename
                current_reference_name = signal.chr
                current_alignment_forward = signal.forward
                current_SV_signals = [signal]
            read_alignment_signals.append(
                ReadAlignmentSignals(
                    samplename=current_samplename,
                    read_name=current_readname,
                    reference_name=current_reference_name,
                    alignment_forward=current_alignment_forward,
                    SV_signals=current_SV_signals,
                )
            )
        return read_alignment_signals


@attrs.define
class SequenceObject:
    name: str
    id: str
    sequence: str
    description: str | None = None
    qualities: list[int] | None = None

    def unstructure(self):
        return cattrs.unstructure(self)

    @classmethod
    def from_biopython(seqobj, seqrecord: SeqRecord):
        name = seqrecord.name
        id = seqrecord.id
        sequence = str(seqrecord.seq)
        description = seqrecord.description
        qualities = seqrecord.letter_annotations.get("phred_quality", None)
        return seqobj(name, id, sequence, description, qualities)


@attrs.define()
class Alignment:
    readname: str
    reference_name: str
    reference_ID: int
    reference_start: int
    reference_end: int
    samdict: dict
    headerdict: dict
    samplename: str | None = None

    @classmethod
    def from_pysam(cls, aln: pysam.AlignedSegment, samplename: str | None = None):
        """Construct an Alignment object from a pysam.AlignedSegment."""
        return cls(
            readname=aln.query_name,
            reference_name=aln.reference_name,
            reference_ID=aln.reference_id,
            reference_start=aln.reference_start,
            reference_end=aln.reference_end,
            samdict=aln.to_dict(),
            headerdict=aln.header.to_dict() if hasattr(aln, "header") else {},
            samplename=samplename,
        )

    def to_pysam(self) -> pysam.AlignedSegment:
        header = pysam.AlignmentHeader.from_dict(self.headerdict)
        return pysam.AlignedSegment.from_dict(self.samdict, header)

    def is_reverse(self) -> bool:
        """Check if the alignment is reverse complemented."""
        return int(self.samdict["flag"]) & 0x10 != 0

    def __eq__(self, value) -> bool:
        # TODO: compare start and cigar strings, not others
        return (
            self.samdict["ref_pos"] == value.samdict["ref_pos"]
            and self.samdict["ref_name"] == value.samdict["ref_name"]
            and self.samdict["flag"] == value.samdict["flag"]
            and self.readname == value.readname
            and self.samdict["cigar"] == value.samdict["cigar"]
        )

    def unstructure(self):
        cattrs.unstructure(self)

    def __hash__(self) -> int:
        return hash((
            self.readname,
            self.reference_name,
            self.reference_ID,
            self.samplename,
            self.samdict["cigar"],
        ))


# given this object and a reference, the original query sequence can be reconstructed
# via consensus_lib.reconstruct_ReconstructibleSequence
# @attrs.define
# class ReconstructibleSequence:
#     alignment:Alignment
#     inserted_subseqs:list[str]
#     reference_name:str
#     description:str=''
#     def unstructure(self):
#         return cattrs.unstructure(self)


# used to merge sv signals in case they are on the same repeat
# all sv signals can be merged to merged_sv_signalsrepeat
# singleton objects: with only one original signal and no repeatID
# merged objects: with multiple original signals and a repeatID
@attrs.define
class MergedSVSignal(SVsignal):
    chr: str
    repeatIDs: list[int]
    original_alt_sequences: list[str]
    original_ref_sequences: list[str]

    def unstructure(self):
        return super().unstructure()

    def __lt__(self, other: object) -> bool:
        return super().__lt__(other)

    def __hash__(self):
        return hash(super().__hash__() + hash(self.chr))

    def get_alt_sequence(self):
        if self.sv_type == 0:
            return "".join(self.original_alt_sequences)[: self.size]
        else:
            return "".join(self.original_alt_sequences)

    def get_ref_sequence(self):
        if self.sv_type == 1:
            return "".join(self.original_ref_sequences)[: self.size]
        else:
            return "".join(self.original_ref_sequences)


# @attrs.define
# class SVcomplex:
#     svprimitives:list[SVprimitive] # list of SVprimitive objects (break ends) that are part of this complex
#     svType:int
#     vcfID:None|str=None
#     alt_sequences:list[str]=[] # list of alternative sequences for the complex
#     ref_sequences:list[str]=[]
#     def unstructure(self):
#         return cattrs.unstructure(self)
#     def get_consensusID(self) -> str:
#         if len(self.svprimitives) == 0:
#             raise ValueError("Cannot get consensus ID from empty SVcomplex")
#         return self.svprimitives[0].consensusID
#     def __lt__(self, other: object) -> bool:
#         # compare reference positions of the minimal SVprimitive in the list
#         if len(self.svprimitives) == 0 or len(other.svprimitives) == 0:
#             raise ValueError("Cannot compare empty SVcomplexes")
#         min_self = min(self.svprimitives,key=lambda x: (x.reference_name, x.ref_start))
#         min_other = min(other.svprimitives,key=lambda x: (x.reference_name, x.ref_start))
#         return min_self.reference_name <= min_other.reference_name
#     def __hash__(self):
#         return hash(sum([hash(p) for p in self.svprimitives]))


# @attrs.define
# class InDelBndCandidate:
#     svType:str # INS DEL BNDL BNDR
#     svSize:int
#     consensus_name:str
#     refChrID:int
#     refChr:str # reference chromosome
#     refStart:int # start position of SV in reference
#     refEnd:int # end position of SV in reference
#     consensus_start:int # start ppsition of SV in consensus
#     consensus_end:int # end position of SV in consensus
#     consensus_length:int # length of consensus sequence
#     depthSupport:dict # of form {sampleID:[depth_l,depth_r,error_l,error_r,support]}
#     repIDs:list # list of IDs of repeats
#     is_reverse:int
#     def unstructure(self):
#         r = cattrs.unstructure(self)
#         r['repIDs'] = list(map(int,r['repIDs']))
#         return r
#     def __eq__(self, __o: object) -> bool:
#         return self.svType == __o.svType \
#             and self.svSize == __o.svSize \
#             and self.consensus_name == __o.consensus_name \
#             and self.refChrID == __o.refChrID \
#             and self.refStart == __o.refStart \
#             and self.refEnd == __o.refEnd \
#             and self.consensus_start == __o.consensus_start \
#             and self.consensus_end == __o.consensus_end \
#             and self.depthSupport == __o.depthSupport \
#             and self.repIDs == __o.repIDs \
#             and self.is_reverse == __o.is_reverse

# @attrs.define
# class IntraSVCandidate:
#     svType:str # INS DEL BND
#     svSize:int
#     consensus_name:str
#     referenceIntervals:list
#     referenceStrings:list
#     consensusIntervals:list
#     consensusStrings:list
#     depthSupport:dict # of form {sampleID:[depth_l,depth_r,error_l,error_r,support]}
#     repIDs:typing.List[int] # list of IDs of repeats
#     def unstructure(self):
#         r = cattrs.unstructure(self)
#         r['repIDs'] = list(map(int,r['repIDs']))
#         return r
#     def __eq__(self, __o: object) -> bool:
#         return self.svType == __o.svType \
#             and self.svSize == __o.svSize \
#             and self.consensus_name == __o.consensus_name \
#             and self.referenceIntervals == __o.referenceIntervals \
#             and self.consensusIntervals == __o.consensusIntervals \
#             and self.depthSupport == __o.depthSupport \
#             and self.repIDs == __o.repIDs
