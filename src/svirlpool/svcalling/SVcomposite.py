from __future__ import annotations

import json

import attrs  # type: ignore
import numpy as np  # type: ignore
from intervaltree import IntervalTree  # type: ignore
from xxhash import xxh64

from ..localassembly import SVpatterns


@attrs.define
class SVcomposite:
    """Composite structural variant from one or more consensus sequences (across samples), built from one or more SVpattern instances."""

    sv_type: type[SVpatterns.SVpatternType]
    svPatterns: list[SVpatterns.SVpatternType]

    def __attrs_post_init__(self):
        if not self.svPatterns:
            raise ValueError("SVcomposite must contain at least one SVpattern")

    def unstructure(self):
        """Unstructure the SVcomposite for serialization."""
        d = SVpatterns.converter.unstructure(self)
        if "sv_type" in d and isinstance(d["sv_type"], type):
            d["sv_type"] = d["sv_type"].__name__
        return d

    @classmethod
    def from_unstructured(cls, data: dict) -> SVcomposite:
        """Create SVcomposite from unstructured data."""
        if "sv_type" in data and isinstance(data["sv_type"], str):
            data["sv_type"] = getattr(SVpatterns, data["sv_type"])
        return SVpatterns.converter.structure(data, cls)

    def to_json(self) -> str:
        """Convert SVcomposite to JSON string."""
        return json.dumps(self.unstructure(), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> SVcomposite:
        """Create SVcomposite from JSON string."""
        data = json.loads(json_str)
        return cls.from_unstructured(data)

    @classmethod
    def from_SVpattern(cls, svPattern: SVpatterns.SVpatternType) -> SVcomposite:
        """Build a composite from a single SVpattern, inferring the sv_type."""
        return cls(sv_type=type(svPattern), svPatterns=[svPattern])

    @classmethod
    def from_SVpatterns(cls, svPatterns: list[SVpatterns.SVpatternType]) -> SVcomposite:
        """Build a complex composite from multiple SVpatterns. All SVpatterns must share the same samplename and consensusID. They are parsed as CPX (complex) SVs."""
        if not svPatterns or len(svPatterns) == 0:
            raise ValueError("Must provide at least one SVpattern")
        # if the sv types of the SVpatterns are INS and DEL, compute the total number of insertions and deletions, where its abs(ins - del)
        if len(svPatterns) == 1:
            return cls(sv_type=type(svPatterns[0]), svPatterns=svPatterns)
        sv_type = type(svPatterns[0])
        all_sv_types = {type(sv) for sv in svPatterns}
        # re-calculate the SV size of horizontally merged INS+DEL runs
        if all_sv_types == {
            SVpatterns.SVpatternInsertion,
            SVpatterns.SVpatternDeletion,
        }:  # compositions of INS and DEL can occur, e.g. in tandem repeats
            ins_count = sum(
                sv.get_size() for sv in svPatterns if isinstance(sv, SVpatterns.SVpatternInsertion)
            )
            del_count = sum(
                sv.get_size() for sv in svPatterns if isinstance(sv, SVpatterns.SVpatternDeletion)
            )
            inserted_bases = ins_count - del_count
            if inserted_bases > 0:
                sv_type = SVpatterns.SVpatternInsertion
            else:
                sv_type = SVpatterns.SVpatternDeletion

        return cls(sv_type=sv_type, svPatterns=svPatterns)

    @property
    def ref_start(self) -> tuple[str, int]:
        return min((p.chr, p.ref_start) for p in self.svPatterns)

    @property
    def ref_end(self) -> tuple[str, int]:
        return max((p.chr, p.ref_end) for p in self.svPatterns)

    @property
    def repeatIDs(self) -> set[int]:
        """Get the repeatIDs from all SVpatterns in the composite."""
        repeatIDs = set()
        for svp in self.svPatterns:
            repeatIDs.update(svp.repeatIDs)
        return repeatIDs

    def __lt__(self, other: SVcomposite) -> bool:
        return self.ref_start < other.ref_start

    def calculate_size(self) -> int:
        """Get the size of the SVcomposite, which is the weighted mean size of the SVpatterns it contains."""
        sizes = [sv.get_size() for sv in self.svPatterns]
        weights = [len(sv.get_supporting_reads()) for sv in self.svPatterns]
        weighted_mean_size = sum(
            size * weight for size, weight in zip(sizes, weights, strict=True)
        ) / sum(weights)
        return int(round(weighted_mean_size))

    def get_size_populations(self) -> list[int]:
        """Returns a list of sizes of distortion signals. Size is neg. for del and pos. for ins."""
        return [
            int(size)
            for svp in self.svPatterns
            if svp.size_distortions is not None
            for size in svp.size_distortions.values()
        ]

    # def get_most_supported_svPattern(self) -> SVpatterns.SVpatternType:
    #     """Get the SVpattern with the most supporting reads."""
    #     return max(self.svPatterns, key=lambda svp: len(svp.get_supporting_reads()))

    def get_regions(self, tolerance_radius: int = 25) -> list[tuple[str, int, int]]:
        if tolerance_radius <= 0:
            raise ValueError("Tolerance radius must be greater than 0")
        return [
            region
            for svp in self.svPatterns
            for region in svp.get_regions(tolerance_radius=tolerance_radius)
        ]

    def get_alt_readnames_per_sample(self) -> dict[str, set[str]]:
        """Get all supporting reads per sample (from each SVpattern, from each SVprimitive) returns a dict {samplename: set(readnames)}"""
        alt_readnames: dict[str, set[str]] = {}
        for svPattern in self.svPatterns:
            alt_readnames.setdefault(svPattern.samplename, set()).update(
                svPattern.get_supporting_reads()
            )
        return alt_readnames

    def get_alt_readnamehashes_per_sample(self) -> dict[str, set[int]]:
        """Get all supporting reads per sample (from each SVpattern, from each SVprimitive) returns a dict {samplename: set(hashed readnames)}"""
        alt_readnames: dict[str, set[int]] = {}
        for svPattern in self.svPatterns:
            alt_readnames.setdefault(svPattern.samplename, set()).update(
                int(xxh64(s).intdigest()) for s in svPattern.get_supporting_reads()
            )
        return alt_readnames

    def get_size(self) -> int:
        # Group SVpatterns by (samplename, consensusID) - same logic as get_alt_sequence
        groups = {}
        for svPattern in self.svPatterns:
            key = (svPattern.samplename, svPattern.consensusID)
            if key not in groups:
                groups[key] = []
            groups[key].append(svPattern)

        if not groups:
            return 0

        # Find the group with maximum supporting unique reads - same logic as get_alt_sequence
        best_group = None
        max_supporting_reads = 0

        for group_patterns in groups.values():
            # Get all unique supporting reads from this group
            all_supporting_reads = set()
            for svPattern in group_patterns:
                all_supporting_reads.update(svPattern.get_supporting_reads())

            if len(all_supporting_reads) > max_supporting_reads:
                max_supporting_reads = len(all_supporting_reads)
                best_group = group_patterns

        if best_group is None:
            return 0

        # Calculate summed size from the best group only
        summed_bp = 0
        for svPattern in best_group:
            value = svPattern.get_size()
            if isinstance(svPattern, SVpatterns.SVpatternDeletion):
                summed_bp -= value
            else:
                summed_bp += value

        return abs(summed_bp)  # Return absolute value for size


    def get_alt_sequence(self) -> str:
        # Group SVpatterns by (samplename, consensusID)
        groups = {}
        for svPattern in self.svPatterns:
            key = (svPattern.samplename, svPattern.consensusID)
            if key not in groups:
                groups[key] = []
            groups[key].append(svPattern)

        if not groups:
            return ""

        # Find the group with maximum supporting unique reads
        best_group = None
        max_supporting_reads = 0

        for group_patterns in groups.values():
            # Get all unique supporting reads from this group
            all_supporting_reads = set()
            for svPattern in group_patterns:
                all_supporting_reads.update(svPattern.get_supporting_reads())

            if len(all_supporting_reads) > max_supporting_reads:
                max_supporting_reads = len(all_supporting_reads)
                best_group = group_patterns

        if best_group is None:
            return ""

        # Concatenate alt sequences from the best group
        concatenated_sequence = ""
        for svPattern in best_group:
            if isinstance(svPattern, SVpatterns.SVpatternInsertion) or isinstance(
                svPattern, SVpatterns.SVpatternInversion
            ):
                seq = svPattern.get_sequence()
                if seq is None:
                    raise ValueError(
                        f"SVpattern {svPattern} in SVcomposite has no alt sequence set."
                    )
                concatenated_sequence += seq
            # Add other pattern types as needed
        size = self.get_size()
        return concatenated_sequence[:size] if size > 0 else ""

    def get_ref_sequence(self) -> str:
        # Group SVpatterns by (samplename, consensusID)
        groups = {}
        for svPattern in self.svPatterns:
            key = (svPattern.samplename, svPattern.consensusID)
            if key not in groups:
                groups[key] = []
            groups[key].append(svPattern)

        if not groups:
            return ""

        # Find the group with maximum supporting unique reads
        best_group = None
        max_supporting_reads = 0

        for group_patterns in groups.values():
            # Get all unique supporting reads from this group
            all_supporting_reads = set()
            for svPattern in group_patterns:
                all_supporting_reads.update(svPattern.get_supporting_reads())

            if len(all_supporting_reads) > max_supporting_reads:
                max_supporting_reads = len(all_supporting_reads)
                best_group = group_patterns

        if best_group is None:
            return ""

        # Concatenate reference sequences from the best group
        concatenated_sequence = ""
        for svPattern in best_group:
            if type(svPattern) is SVpatterns.SVpatternDeletion:
                seq = svPattern.get_sequence()
                if seq is None:
                    raise ValueError(
                        f"SVpattern {svPattern} in SVcomposite has no ref sequence set."
                    )
                concatenated_sequence += seq
            # Add other pattern types as needed

        size = self.get_size()
        return concatenated_sequence[:size] if size > 0 else "" 

    def overlaps_any(self, other: SVcomposite, tolerance_radius: int = 50) -> bool:
        if self.repeatIDs.intersection(other.repeatIDs):
            return True
        regions_a = self.get_regions(tolerance_radius=tolerance_radius)
        regions_b = other.get_regions(tolerance_radius=tolerance_radius)
        intervaltrees: dict[str, IntervalTree] = {}
        for region in regions_a:
            if region[0] not in intervaltrees:
                intervaltrees[region[0]] = IntervalTree()
            intervaltrees[region[0]].addi(
                region[1], region[2], {0}
            )  # add 0 as it is the identifier for the first SVcomposite
        for region in regions_b:
            if region[0] not in intervaltrees:
                intervaltrees[region[0]] = IntervalTree()
            intervaltrees[region[0]].addi(
                region[1], region[2], {1}
            )  # add 1 as it is the identifier for the second SVcomposite
        for _chr, tree in intervaltrees.items():
            tree.merge_overlaps(data_reducer=lambda x, y: x.union(y))
            # then iterate all merged intervals in the tree. If data is {0,1}, then there is an overlap
            for interval in tree:
                if len(interval.data) == 2:
                    return True
        return False

    def get_inserted_sequences(self) -> list[str]:
        """Get the inserted sequences of all SVpatterns in the composite."""
        alt_sequences = []
        for svPattern in self.svPatterns:
            seq = ""
            if isinstance(svPattern, SVpatterns.SVpatternInsertion) or isinstance(
                svPattern, SVpatterns.SVpatternInversion
            ):
                seq = svPattern.get_sequence()
            else:
                pass
            if seq is None:
                raise ValueError(
                    f"SVpattern {svPattern} in SVcomposite has no alt sequence set."
                )
            alt_sequences.append(seq)
        return alt_sequences

    def get_inserted_complexity_tracks(self) -> list[np.ndarray]:
        """Get the complexity tracks of all SVpatterns in the composite. Analogous to get_inserted_sequences."""
        complexity_tracks = []
        for svPattern in self.svPatterns:
            if isinstance(svPattern, SVpatterns.SVpatternInsertion) or isinstance(
                svPattern, SVpatterns.SVpatternInversion
            ):
                comp_track = svPattern.get_sequence_complexity()
                complexity_tracks.append(comp_track)
        return complexity_tracks

    def get_reference_sequences(self) -> list[str]:
        """Get the deleted sequences of all SVpatterns in the composite."""
        ref_sequences = []
        for svPattern in self.svPatterns:
            seq = ""
            if isinstance(svPattern, SVpatterns.SVpatternDeletion):
                seq = svPattern.get_sequence()
                ref_sequences.append(seq)
            else:
                pass
        return ref_sequences

    def get_reference_complexity_tracks(self) -> list[np.ndarray]:
        """Get the complexity tracks of all SVpatterns in the composite. Analogous to get_reference_sequences."""
        complexity_tracks = []
        for svPattern in self.svPatterns:
            if isinstance(svPattern, SVpatterns.SVpatternDeletion):
                comp_track = svPattern.get_sequence_complexity()
                complexity_tracks.append(comp_track)
        return complexity_tracks
