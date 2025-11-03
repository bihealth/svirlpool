# import cattrs, attrs
# from . import datatypes, genotyping


# @attrs.define
# class Adjacency:
#     chr: str
#     ref_pos: int
#     gapsize: int
#     tp_str: str  # needs to be of form t[p[ or t]p] or [p[t or ]p]t; p is chr:start, t is alt seq
#     mate_id: str

#     def unstructure(self):
#         return cattrs.unstructure(self)


# @attrs.define
# class SVprimitive(datatypes.MergedSVSignal):  # can be ins,del,bndl,bndr
#     samplename: str
#     consensusID: str
#     alignmentID: int  # ID of the alignment to the reference
#     svID: int  # ID of the SV primitive on the consensus alignment
#     reference_name: str
#     aln_is_reverse: bool
#     distortions: list[tuple[int, int, float]]  # list of abstracted SV signals parsed from the cut read to consensus alignments with (svtype, svsize, distance to this call)
#     genotypeMeasurement: genotyping.GenotypeMeasurement | None = None  # samplename: GenotypeMeasurement
#     adjacent_bnd: None | Adjacency = None

#     def unstructure(self):
#         return super().unstructure()

#     def __lt__(self, other: object) -> bool:
#         return super().__lt__(other)

#     def __hash__(self) -> int:
#         return super().__hash__()

#     def get_vcfID(self) -> str:
#         return f"{datatypes.SV_TYPE_DICT[self.sv_type]}.{self.consensusID}__{str(self.alignmentID)}.{str(self.svID)}"
