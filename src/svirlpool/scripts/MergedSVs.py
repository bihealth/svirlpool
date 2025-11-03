# """
# MergedSVs module for unifying SVprimitive and SVcomplex objects.

# This module provides the MergedSVs class that can contain and merge both
# SVprimitive and SVcomplex objects, providing a unified interface for
# handling structural variants of varying complexity.
# """

# import typing
# import cattrs
# from . import datatypes


# class MergedSVs:
#     """
#     Unified container for SVprimitives and SVcomplexes that can merge with other MergedSVs.
#     Does not use attrs to allow for custom serialization/deserialization logic.
#     """

#     def __init__(self,
#                  sv_type: int,
#                  size: int,
#                  chr: str,
#                  ref_start: int,
#                  ref_end: int,
#                  svPrimitives: list = None,
#                  svComplexes: list = None,
#                  merged_from: list = None,
#                  consensusIDs: list = None,
#                  reference_names: set = None,
#                  vcfID: str = None,
#                  total_support: int = 0,
#                  confidence_score: float = 0.0,
#                  complexity_level: int = 0):

#         self.sv_type = sv_type
#         self.size = size
#         self.chr = chr
#         self.ref_start = ref_start
#         self.ref_end = ref_end

#         # Component tracking
#         self.svPrimitives = svPrimitives or []
#         self.svComplexes = svComplexes or []
#         self.merged_from = merged_from or []

#         # Metadata
#         self.consensusIDs = consensusIDs or []
#         self.reference_names = reference_names or set()
#         self.vcfID = vcfID

#         # Quality metrics
#         self.total_support = total_support
#         self.confidence_score = confidence_score
#         self.complexity_level = complexity_level

#         # Auto-update derived properties
#         self._update_derived_properties()

#     def _update_derived_properties(self):
#         """Update derived properties based on current components"""
#         # Update consensusIDs from components
#         consensus_ids = set(self.consensusIDs)

#         for svp in self.svPrimitives:
#             if svp.consensusID:
#                 consensus_ids.add(svp.consensusID)
#             if svp.reference_name:
#                 self.reference_names.add(svp.reference_name)

#         for svc in self.svComplexes:
#             try:
#                 consensus_ids.add(svc.get_consensusID())
#             except ValueError:
#                 pass  # Empty SVcomplex
#             for svp in svc.svprimitives:
#                 if svp.reference_name:
#                     self.reference_names.add(svp.reference_name)

#         self.consensusIDs = list(consensus_ids)

#         # Update support count
#         self.total_support = self._calculate_total_support()

#         # Update complexity level
#         if self.svComplexes:
#             self.complexity_level = max(self.complexity_level, 1)
#         if len(self.svComplexes) > 1 or len(self.svPrimitives) > 2:
#             self.complexity_level = max(self.complexity_level, 2)

#     def _calculate_total_support(self) -> int:
#         """Calculate total support from all components"""
#         total = 0

#         # From SVprimitives
#         for svp in self.svPrimitives:
#             for gm in svp.genotypeMeasurements.values():
#                 total += len(gm.supporting_reads)

#         # From SVcomplexes
#         for svc in self.svComplexes:
#             for svp in svc.svprimitives:
#                 for gm in svp.genotypeMeasurements.values():
#                     total += len(gm.supporting_reads)

#         return total

#     @classmethod
#     def from_primitive(cls, svp: 'SVprimitive_class.SVprimitive') -> 'MergedSVs':
#         """Create MergedSVs from a single SVprimitive"""
#         return cls(
#             sv_type=svp.sv_type,
#             size=svp.size,
#             chr=svp.chr,
#             ref_start=svp.ref_start,
#             ref_end=svp.ref_end,
#             svPrimitives=[svp],
#             consensusIDs=[svp.consensusID] if svp.consensusID else [],
#             reference_names={svp.reference_name} if svp.reference_name else set(),
#             complexity_level=0
#         )

#     @classmethod
#     def from_complex(cls, svc: 'datatypes.SVcomplex') -> 'MergedSVs':
#         """Create MergedSVs from a single SVcomplex"""
#         if not svc.svprimitives:
#             raise ValueError("Cannot create MergedSVs from empty SVcomplex")

#         # Calculate bounds from all primitives
#         ref_starts = [svp.ref_start for svp in svc.svprimitives]
#         ref_ends = [svp.ref_end for svp in svc.svprimitives]

#         # Use the first primitive's chromosome (assuming all are on same chr for complex)
#         chr_name = svc.svprimitives[0].chr

#         # Estimate size based on complex type
#         estimated_size = max(ref_ends) - min(ref_starts)

#         return cls(
#             sv_type=svc.svType,
#             size=estimated_size,
#             chr=chr_name,
#             ref_start=min(ref_starts),
#             ref_end=max(ref_ends),
#             svComplexes=[svc],
#             consensusIDs=[svc.get_consensusID()],
#             reference_names={svp.reference_name for svp in svc.svprimitives if svp.reference_name},
#             complexity_level=1
#         )

#     def get_all_genotype_measurements(self) -> dict[str, list]:
#         """Collect genotype measurements from all components"""
#         measurements = {}

#         # From SVprimitives
#         for svp in self.svPrimitives:
#             for sample, gm in svp.genotypeMeasurements.items():
#                 if sample not in measurements:
#                     measurements[sample] = []
#                 measurements[sample].append(gm)

#         # From SVcomplexes
#         for svc in self.svComplexes:
#             for svp in svc.svprimitives:
#                 for sample, gm in svp.genotypeMeasurements.items():
#                     if sample not in measurements:
#                         measurements[sample] = []
#                     measurements[sample].append(gm)

#         return measurements

#     def merge_with_mergedsv(self, other: 'MergedSVs') -> 'MergedSVs':
#         """Merge with another MergedSVs"""
#         return MergedSVs(
#             sv_type=self._resolve_sv_type_with_mergedsv(other),
#             size=self._estimate_merged_size(other),
#             chr=self.chr,  # Assuming same chromosome
#             ref_start=min(self.ref_start, other.ref_start),
#             ref_end=max(self.ref_end, other.ref_end),
#             svPrimitives=self.svPrimitives + other.svPrimitives,
#             svComplexes=self.svComplexes + other.svComplexes,
#             merged_from=[self, other],
#             consensusIDs=list(set(self.consensusIDs + other.consensusIDs)),
#             reference_names=self.reference_names | other.reference_names,
#             total_support=self.total_support + other.total_support,
#             complexity_level=max(self.complexity_level, other.complexity_level)
#         )

#     def add_svprimitive(self, svp: 'SVprimitive_class.SVprimitive') -> 'MergedSVs':
#         """Add an SVprimitive and return new MergedSVs"""
#         return MergedSVs(
#             sv_type=self._resolve_sv_type_with_primitive(svp),
#             size=self._estimate_size_with_primitive(svp),
#             chr=self.chr,
#             ref_start=min(self.ref_start, svp.ref_start),
#             ref_end=max(self.ref_end, svp.ref_end),
#             svPrimitives=self.svPrimitives + [svp],
#             svComplexes=self.svComplexes.copy(),
#             consensusIDs=list(set(self.consensusIDs + ([svp.consensusID] if svp.consensusID else []))),
#             reference_names=self.reference_names | ({svp.reference_name} if svp.reference_name else set()),
#             complexity_level=max(self.complexity_level, 0)
#         )

#     def add_svcomplex(self, svc: 'datatypes.SVcomplex') -> 'MergedSVs':
#         """Add an SVcomplex and return new MergedSVs"""
#         if not svc.svprimitives:
#             raise ValueError("Cannot add empty SVcomplex")

#         ref_starts = [svp.ref_start for svp in svc.svprimitives]
#         ref_ends = [svp.ref_end for svp in svc.svprimitives]

#         return MergedSVs(
#             sv_type=self._resolve_sv_type_with_complex(svc),
#             size=self._estimate_size_with_complex(svc),
#             chr=self.chr,
#             ref_start=min(self.ref_start, min(ref_starts)),
#             ref_end=max(self.ref_end, max(ref_ends)),
#             svPrimitives=self.svPrimitives.copy(),
#             svComplexes=self.svComplexes + [svc],
#             consensusIDs=list(set(self.consensusIDs + [svc.get_consensusID()])),
#             reference_names=self.reference_names | {svp.reference_name for svp in svc.svprimitives if svp.reference_name},
#             complexity_level=max(self.complexity_level, 1)
#         )

#     def can_merge_with(self, other) -> bool:
#         """Determine if this MergedSVs can merge with another object"""
#         if isinstance(other, MergedSVs):
#             return self._can_merge_with_mergedsv(other)
#         elif isinstance(other, SVprimitive_class.SVprimitive):
#             return self._can_merge_with_primitive(other)
#         elif isinstance(other, datatypes.SVcomplex):
#             return self._can_merge_with_complex(other)
#         return False

#     def get_alt_sequences(self) -> list[str]:
#         """Get all alternative sequences from components"""
#         alt_sequences = []

#         # From SVprimitives
#         for svp in self.svPrimitives:
#             alt_seq = svp.get_alt_sequence()
#             if alt_seq:
#                 alt_sequences.append(alt_seq)

#         # From SVcomplexes
#         for svc in self.svComplexes:
#             alt_sequences.extend(svc.alt_sequences)

#         return alt_sequences

#     def get_ref_sequences(self) -> list[str]:
#         """Get all reference sequences from components"""
#         ref_sequences = []

#         # From SVprimitives
#         for svp in self.svPrimitives:
#             ref_seq = svp.get_ref_sequence()
#             if ref_seq:
#                 ref_sequences.append(ref_seq)

#         # From SVcomplexes
#         for svc in self.svComplexes:
#             ref_sequences.extend(svc.ref_sequences)

#         return ref_sequences

#     # Helper methods for merging logic
#     def _resolve_sv_type_with_mergedsv(self, other: 'MergedSVs') -> int:
#         """Resolve SV type when merging with another MergedSVs"""
#         # Prioritize complex types over simple ones
#         if self.complexity_level > other.complexity_level:
#             return self.sv_type
#         elif other.complexity_level > self.complexity_level:
#             return other.sv_type
#         else:
#             # Same complexity level, use the one with higher support
#             return self.sv_type if self.total_support >= other.total_support else other.sv_type

#     def _resolve_sv_type_with_primitive(self, svp: 'SVprimitive_class.SVprimitive') -> int:
#         """Resolve SV type when adding a primitive"""
#         # If current is complex, keep complex type
#         if self.complexity_level > 0:
#             return self.sv_type
#         # If both are simple, use the one with higher support
#         svp_support = sum(len(gm.supporting_reads) for gm in svp.genotypeMeasurements.values())
#         return self.sv_type if self.total_support >= svp_support else svp.sv_type

#     def _resolve_sv_type_with_complex(self, svc: 'datatypes.SVcomplex') -> int:
#         """Resolve SV type when adding a complex"""
#         # Complex types have priority
#         return svc.svType

#     def _estimate_merged_size(self, other: 'MergedSVs') -> int:
#         """Estimate size after merging"""
#         # Use the span of all components
#         new_span = max(self.ref_end, other.ref_end) - min(self.ref_start, other.ref_start)
#         # Take the maximum of span and individual sizes
#         return max(new_span, self.size, other.size)

#     def _estimate_size_with_primitive(self, svp: 'SVprimitive_class.SVprimitive') -> int:
#         """Estimate size when adding a primitive"""
#         new_span = max(self.ref_end, svp.ref_end) - min(self.ref_start, svp.ref_start)
#         return max(new_span, self.size, svp.size)

#     def _estimate_size_with_complex(self, svc: 'datatypes.SVcomplex') -> int:
#         """Estimate size when adding a complex"""
#         ref_starts = [svp.ref_start for svp in svc.svprimitives]
#         ref_ends = [svp.ref_end for svp in svc.svprimitives]
#         complex_span = max(ref_ends) - min(ref_starts)
#         new_span = max(self.ref_end, max(ref_ends)) - min(self.ref_start, min(ref_starts))
#         return max(new_span, self.size, complex_span)

#     def _can_merge_with_mergedsv(self, other: 'MergedSVs') -> bool:
#         """Check if two MergedSVs can be merged using existing SVprimitive logic"""
#         # Use similar logic as can_merge_svPrimitives_pair from sv_calling.py
#         # For now, simple distance and chromosome check
#         if self.chr != other.chr:
#             return False

#         # Distance check (using similar parameters as in sv_calling.py)
#         distance = abs(self.ref_start - other.ref_start)
#         tolerance = 1000  # Similar to 'close' parameter

#         if distance > tolerance:
#             return False

#         # Size similarity check
#         size_diff = abs(self.size - other.size)
#         min_size = min(self.size, other.size)
#         if min_size > 0 and size_diff / min_size > 0.5:  # 50% size difference threshold
#             return False

#         return True

#     def _can_merge_with_primitive(self, svp: 'SVprimitive_class.SVprimitive') -> bool:
#         """Check if can merge with a primitive"""
#         if self.chr != svp.chr:
#             return False

#         distance = abs(self.ref_start - svp.ref_start)
#         if distance > 1000:
#             return False

#         size_diff = abs(self.size - svp.size)
#         min_size = min(self.size, svp.size)
#         if min_size > 0 and size_diff / min_size > 0.5:
#             return False

#         return True

#     def _can_merge_with_complex(self, svc: 'datatypes.SVcomplex') -> bool:
#         """Check if can merge with a complex"""
#         if not svc.svprimitives:
#             return False

#         # Check chromosome compatibility
#         svc_chr = svc.svprimitives[0].chr
#         if self.chr != svc_chr:
#             return False

#         # Check distance to closest primitive in complex
#         min_distance = float('inf')
#         for svp in svc.svprimitives:
#             distance = abs(self.ref_start - svp.ref_start)
#             min_distance = min(min_distance, distance)

#         return min_distance <= 1000

#     def to_dict(self) -> dict:
#         """Convert MergedSVs to dictionary representation"""
#         return {
#             'sv_type': self.sv_type,
#             'size': self.size,
#             'chr': self.chr,
#             'ref_start': self.ref_start,
#             'ref_end': self.ref_end,
#             'svPrimitives': [svp.unstructure() for svp in self.svPrimitives],
#             'svComplexes': [svc.unstructure() for svc in self.svComplexes],
#             'merged_from': [msv.to_dict() for msv in self.merged_from],
#             'consensusIDs': self.consensusIDs,
#             'reference_names': list(self.reference_names),
#             'vcfID': self.vcfID,
#             'total_support': self.total_support,
#             'confidence_score': self.confidence_score,
#             'complexity_level': self.complexity_level
#         }

#     @classmethod
#     def from_dict(cls, data: dict) -> 'MergedSVs':
#         """Create MergedSVs from dictionary representation"""
#         # Reconstruct SVprimitives
#         svprimitives = []
#         for svp_data in data.get('svPrimitives', []):
#             svprimitives.append(cattrs.structure(svp_data, SVprimitive_class.SVprimitive))

#         # Reconstruct SVcomplexes
#         svComplexes = []
#         for svc_data in data.get('svComplexes', []):
#             svComplexes.append(cattrs.structure(svc_data, datatypes.SVcomplex))

#         # Reconstruct merged_from (recursive)
#         merged_from = []
#         for msv_data in data.get('merged_from', []):
#             merged_from.append(cls.from_dict(msv_data))

#         return cls(
#             sv_type=data['sv_type'],
#             size=data['size'],
#             chr=data['chr'],
#             ref_start=data['ref_start'],
#             ref_end=data['ref_end'],
#             svPrimitives=svprimitives,
#             svComplexes=svComplexes,
#             merged_from=merged_from,
#             consensusIDs=data.get('consensusIDs', []),
#             reference_names=set(data.get('reference_names', [])),
#             vcfID=data.get('vcfID'),
#             total_support=data.get('total_support', 0),
#             confidence_score=data.get('confidence_score', 0.0),
#             complexity_level=data.get('complexity_level', 0)
#         )

#     def to_vcf_representation(self) -> dict:
#         """Convert to VCF-compatible representation"""
#         # Get all genotype measurements
#         genotype_measurements = self.get_all_genotype_measurements()

#         # Determine VCF fields
#         vcf_repr = {
#             'CHROM': self.chr,
#             'POS': self.ref_start + 1,  # VCF is 1-based
#             'ID': self.vcfID or f"MergedSV_{self.chr}_{self.ref_start}",
#             'REF': self.get_ref_sequences()[0] if self.get_ref_sequences() else 'N',
#             'ALT': self.get_alt_sequences()[0] if self.get_alt_sequences() else f'<{datatypes.SV_TYPE_DICT.get(self.sv_type, "UNK")}>',
#             'QUAL': 60,  # Default quality
#             'FILTER': 'PASS',
#             'INFO': {
#                 'SVTYPE': datatypes.SV_TYPE_DICT.get(self.sv_type, 'UNK'),
#                 'SVLEN': self.size if self.sv_type != 1 else -self.size,  # Negative for deletions
#                 'END': self.ref_end,
#                 'SUPPORT': self.total_support,
#                 'COMPLEXITY': self.complexity_level,
#                 'CONSENSUSIDS': ','.join(self.consensusIDs),
#                 'NPRIMITIVES': len(self.svPrimitives),
#                 'NCOMPLEXES': len(self.svComplexes)
#             }
#         }

#         # Add sample-specific information
#         if genotype_measurements:
#             vcf_repr['FORMAT'] = ['GT', 'GQ', 'DR', 'DV']
#             vcf_repr['SAMPLES'] = {}

#             for sample, measurements in genotype_measurements.items():
#                 # Aggregate measurements per sample
#                 total_support = sum(len(gm.supporting_reads) for gm in measurements)
#                 total_depth = sum(gm.estimated_total_depth for gm in measurements)

#                 # Simple genotype calling
#                 if total_support >= 3:  # Threshold for variant call
#                     genotype = '0/1' if total_support < total_depth * 0.8 else '1/1'
#                 else:
#                     genotype = '0/0'

#                 vcf_repr['SAMPLES'][sample] = {
#                     'GT': genotype,
#                     'GQ': min(60, total_support * 10),  # Quality score
#                     'DR': total_depth - total_support,  # Reference depth
#                     'DV': total_support  # Variant depth
#                 }

#         return vcf_repr

#     def __str__(self) -> str:
#         return f"MergedSVs(type={datatypes.SV_TYPE_DICT.get(self.sv_type, self.sv_type)}, size={self.size}, pos={self.chr}:{self.ref_start}-{self.ref_end}, primitives={len(self.svPrimitives)}, complexes={len(self.svComplexes)})"

#     def __repr__(self) -> str:
#         return self.__str__()

#     def __eq__(self, other) -> bool:
#         if not isinstance(other, MergedSVs):
#             return False
#         return (self.sv_type == other.sv_type and
#                 self.chr == other.chr and
#                 self.ref_start == other.ref_start and
#                 self.ref_end == other.ref_end and
#                 self.size == other.size)

#     def __hash__(self) -> int:
#         return hash((self.sv_type, self.chr, self.ref_start, self.ref_end, self.size))

#     def __lt__(self, other) -> bool:
#         if not isinstance(other, MergedSVs):
#             return False
#         if self.chr != other.chr:
#             return self.chr < other.chr
#         return self.ref_start < other.ref_start


# def merge_merged_svs_list(merged_svs_list: list[MergedSVs],
#                          far: int = 10_000,
#                          close: int = 1_000,
#                          very_close: int = 50) -> list[MergedSVs]:
#     """
#     Merge a list of MergedSVs objects using similar logic to merge_svPrimitives.

#     Args:
#         merged_svs_list: List of MergedSVs objects to merge
#         far: Maximum distance for merging with repeat overlap
#         close: Maximum distance for merging without repeat overlap
#         very_close: Maximum distance for very close variants

#     Returns:
#         List of merged MergedSVs objects
#     """
#     if not merged_svs_list:
#         return []

#     # Sort by chromosome and position
#     sorted_svs = sorted(merged_svs_list, key=lambda x: (x.chr, x.ref_start))

#     merged_list = []
#     i = 0

#     while i < len(sorted_svs):
#         current = sorted_svs[i]
#         merged_with_any = False

#         # Try to merge with existing merged SVs
#         for j, existing in enumerate(merged_list):
#             if current.can_merge_with(existing):
#                 merged_list[j] = existing.merge_with_mergedsv(current)
#                 merged_with_any = True
#                 break

#         # If not merged with any existing, add as new
#         if not merged_with_any:
#             merged_list.append(current)

#         i += 1

#     return merged_list


# def create_merged_svs_from_mixed_list(sv_objects: list) -> list[MergedSVs]:
#     """
#     Create MergedSVs objects from a mixed list of SVprimitives and SVcomplexes.

#     Args:
#         sv_objects: List containing SVprimitive and/or SVcomplex objects

#     Returns:
#         List of MergedSVs objects
#     """
#     merged_svs_list = []

#     for obj in sv_objects:
#         if isinstance(obj, SVprimitive_class.SVprimitive):
#             merged_svs_list.append(MergedSVs.from_primitive(obj))
#         elif isinstance(obj, datatypes.SVcomplex):
#             merged_svs_list.append(MergedSVs.from_complex(obj))
#         else:
#             raise ValueError(f"Unsupported object type: {type(obj)}")

#     return merged_svs_list
