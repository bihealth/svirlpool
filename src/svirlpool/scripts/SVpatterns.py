from __future__ import annotations

import logging
import pickle
import sqlite3
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path

import attrs  # type: ignore
import cattrs
import numpy as np
from Bio.Seq import Seq  # type: ignore
from intervaltree import Interval  # type: ignore

from .consensus_class import Consensus
from .SVprimitives import SVprimitive
from .util import complexity_local_track, exponential_weight

# Initialize logger
log = logging.getLogger(__name__)


# region SVpattern
@attrs.define
class SVpattern(ABC):
    """Base class for all SV-specific patterns."""

    SVprimitives: list[SVprimitive]
    size_distortions: dict[str, float] | None = (
        None  # is set after filtering by similar regions on the consensus sequence
    )

    def __attrs_post_init__(self):
        if not self.SVprimitives or len(self.SVprimitives) == 0:
            raise ValueError("SVpattern must contain at least one SVprimitive")

    @abstractmethod
    def get_sv_type(self) -> str:
        pass

    @abstractmethod
    def get_size(self) -> int:
        pass

    @property
    def consensusID(self) -> str:
        return self.SVprimitives[0].consensusID

    @property
    def repeatIDs(self) -> set[int]:
        return {
            repID
            for svp in self.SVprimitives
            for repID in svp.repeatIDs
            if repID is not None
        }

    @property
    def samplename(self) -> str:
        return self.SVprimitives[0].samplename

    @property
    def chr(self) -> str:
        return self.SVprimitives[0].chr

    @property
    def ref_start(self) -> int:
        return min(p.ref_start for p in self.SVprimitives)

    @property
    def ref_end(self) -> int:
        return max(p.ref_end for p in self.SVprimitives)

    @property
    def read_start(self) -> int:
        return min(p.read_start for p in self.SVprimitives)

    @property
    def read_end(self) -> int:
        return max(p.read_end for p in self.SVprimitives)

    @property
    def samplenamed_consensusID(self) -> str:
        return self.samplename + ":" + self.consensusID

    def get_consensus_aln_intervals_from_svPrimitives(
        self,
    ) -> list[tuple[str, int, int]]:
        """Gets the consensus alignment intervals from the SVprimitives."""
        return [
            svp.consensus_aln_interval
            for svp in self.SVprimitives
            if svp.consensus_aln_interval is not None
        ]

    def get_similar_sequence_intervals_from_svPrimitives(self) -> list[list[int, int]]:
        """gets all similar_sequence_intervals_on_consensus from the SVprimitives and merges them."""
        intervals = sorted([
            interval
            for svp in self.SVprimitives
            for interval in svp.similar_sequence_intervals_on_consensus
        ])
        # since the intervals have different lengths, the merging window can only advance once the next interval starts after the previous one ends
        if len(intervals) == 0:
            return [[svp.read_start, svp.read_end] for svp in self.SVprimitives]
        current_start = intervals[0][0]
        merged_intervals = []
        current_end = intervals[0][1]
        for start, end in intervals[1:]:
            if start > current_end:
                merged_intervals.append([current_start, current_end])
                current_start = start
                current_end = end
            else:
                current_end = max(current_end, end)
        merged_intervals.append([current_start, current_end])
        return merged_intervals

    def get_supporting_reads(self) -> list[str]:
        start = [
            readname
            for svp in self.SVprimitives
            for readname in svp.genotypeMeasurement.supporting_reads_start
        ]
        try:
            end = []
            for svp in self.SVprimitives:
                if svp.genotypeMeasurement.supporting_reads_end is not None:
                    end.extend(svp.genotypeMeasurement.supporting_reads_end)
            # end = [readname for svp in self.SVprimitives for readname in svp.genotypeMeasurement.supporting_reads_end if svp.genotypeMeasurement.supporting_reads_end is not None]
        except AttributeError:
            raise AttributeError(
                "Some SVprimitives do not have a genotypeMeasurement with supporting_reads_end. self= "
                + str(self)
            )
        if not end:
            return list(set(start))
        return list(set(start + end))

    # def get_total_coverage(self) -> int:
    #     return max([svp.get_total_coverage() for svp in self.SVprimitives])

    def to_json(self) -> str:
        """Convert SVpattern to JSON string."""
        import json

        return json.dumps(converter.unstructure(self), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> SVpatternType:
        """Create SVpattern from JSON string."""
        import json

        data = json.loads(json_str)
        return converter.structure(data, SVpatternType)

    def get_regions(self, tolerance_radius: int = 25) -> list[tuple[str, int, int]]:
        # gets all regions from all SVprimitives
        regions: list[tuple[str, int, int]] = []
        for svp in self.SVprimitives:
            if svp.sv_type == 1 or svp.sv_type == 2:
                regions.append((svp.chr, svp.ref_start, svp.ref_end))
            else:
                regions.append((
                    svp.chr,
                    max(0, svp.ref_start - tolerance_radius),
                    svp.ref_end + tolerance_radius,
                ))
        return regions


@attrs.define
class SVpatternSingleBreakend(SVpattern):
    def get_sv_type(self) -> str:
        return "BND"


@attrs.define
class SVpatternDeletion(SVpattern):
    deleted_sequence: bytes | None = None
    sequence_complexity: bytes | None = None  # Pickled numpy array of float32

    def get_sv_type(self) -> str:
        return "DEL"

    def is_inter_alignment(self) -> bool:
        """Check if two adjacent BNDs are present by checking if the svPrimitives are sv_type 3 or 4"""
        if all(svp.sv_type == 3 or svp.sv_type == 4 for svp in self.SVprimitives):
            return True
        return False

    def set_sequence(
        self, sequence: str, sequence_complexity_max_length: int = 300
    ) -> None:
        """Set the deleted sequence as a pickled byte string."""
        self.deleted_sequence = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length:
            dna_iter = iter(sequence)
            self.sequence_complexity = pickle.dumps(
                complexity_local_track(
                    dna_iter=dna_iter, w=11, K=[2, 3, 4, 5], padding=True
                )
            )
        else:
            # Create dummy placeholder with 1.0 values for long sequences
            dummy_complexity = np.ones(len(sequence), dtype=np.float16)
            self.sequence_complexity = pickle.dumps(dummy_complexity)

    def get_sequence(self) -> str | None:
        """Retrieve the deleted sequence by unpickling."""
        if self.deleted_sequence is None:
            return None
        try:
            return pickle.loads(self.deleted_sequence)
        except Exception:
            return None

    def get_sequence_complexity(self) -> np.ndarray | None:
        """Retrieve the sequence complexity scores."""
        if self.sequence_complexity is None:
            return None
        try:
            return pickle.loads(self.sequence_complexity)
        except Exception:
            return None

    def get_mean_complexity(self) -> float | None:
        """Get the mean complexity score across the sequence."""
        complexity = self.get_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.mean(complexity))
        return None

    def get_min_complexity(self) -> float | None:
        """Get the minimum complexity score in the sequence."""
        complexity = self.get_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.min(complexity))
        return None

    def get_reference_region(self) -> tuple[str, int, int]:
        chr = self.SVprimitives[0].chr
        start = min(svp.ref_start for svp in self.SVprimitives)
        end = max(svp.ref_end for svp in self.SVprimitives)
        if start > end:
            raise ValueError(
                f"SVpatternDeletion has inconsistent reference coordinates! This one does: {self}"
            )
        return (chr, start, end)

    def get_size(self) -> int:
        return abs(self.SVprimitives[-1].ref_end - self.SVprimitives[0].ref_start)

    def get_combined_reference_sequences(self) -> list[str]:
        if all(svp.sv_type == 3 or svp.sv_type == 4 for svp in self.SVprimitives):
            if self.deleted_sequence is None:
                raise ValueError("Deleted sequence is not set for SVpatternDeletion")
            return [pickle.loads(self.deleted_sequence)]
        return [
            seq
            for svp in self.SVprimitives
            for seq in svp.original_ref_sequences
            if svp.original_ref_sequences is not None and svp.sv_type == 1
        ]


@attrs.define
class SVpatternInsertion(SVpattern):
    inserted_sequence: bytes | None = None
    sequence_complexity: bytes | None = None  # Pickled numpy array of float32

    def get_sv_type(self) -> str:
        return "INS"

    def set_sequence(
        self, sequence: str, sequence_complexity_max_length: int = 300
    ) -> None:
        """Set the inserted sequence and compute its complexity."""
        self.inserted_sequence = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length:
            dna_iter = iter(sequence)
            self.sequence_complexity = pickle.dumps(
                complexity_local_track(
                    dna_iter=dna_iter, w=11, K=[2, 3, 4, 5], padding=True
                )
            )
        else:
            # Create dummy placeholder with 1.0 values for long sequences
            dummy_complexity = np.ones(len(sequence), dtype=np.float16)
            self.sequence_complexity = pickle.dumps(dummy_complexity)

    def get_sequence(self) -> str | None:
        """Retrieve the inserted sequence by unpickling."""
        if self.inserted_sequence is None:
            return None
        try:
            return pickle.loads(self.inserted_sequence)
        except Exception:
            return None

    def get_sequence_complexity(self) -> np.ndarray | None:
        """Retrieve the sequence complexity scores."""
        if self.sequence_complexity is None:
            return None
        try:
            return pickle.loads(self.sequence_complexity)
        except Exception:
            return None

    def get_mean_complexity(self) -> float | None:
        """Get the mean complexity score across the sequence."""
        complexity = self.get_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.mean(complexity))
        return None

    def get_min_complexity(self) -> float | None:
        """Get the minimum complexity score in the sequence."""
        complexity = self.get_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.min(complexity))
        return None

    def is_inter_alignment(self) -> bool:
        if all(svp.sv_type == 3 or svp.sv_type == 4 for svp in self.SVprimitives):
            return True
        return False

    def get_sequence_from_consensus(self, consensus: Consensus) -> str:
        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        s = self.SVprimitives[0].read_start - core_sequence_start
        e = (
            self.SVprimitives[-1].read_end - core_sequence_start
        )  # this applies to bot a single insertion or a bnd->insertion
        seq = consensus.consensus_sequence[s:e]
        if self.SVprimitives[0].aln_is_reverse:
            seq = str(Seq(seq).reverse_complement())
        return seq

    def get_size(self) -> int:
        return abs(self.SVprimitives[-1].read_end - self.SVprimitives[0].read_start)

    def get_reference_region(self) -> tuple[str, int, int]:
        # get the svp that has min(chr, ref_start)
        first = min(self.SVprimitives, key=lambda svp: (svp.chr, svp.ref_start))
        last = max(self.SVprimitives, key=lambda svp: (svp.chr, svp.ref_end))
        if first.chr != last.chr:
            raise ValueError(
                f"SVpatternInsertion cannot span over different chromosomes! This one does: {self}"
            )
        if first.ref_start > last.ref_end:
            raise ValueError(
                f"SVpatternInsertion has inconsistent reference coordinates! This one does: {self}"
            )
        return first.chr, first.ref_start, last.ref_end

    def get_combined_alt_sequences(self) -> list[str]:
        # check if the svPrimitives are of type break ends. If so, return inserted sequence
        if all(svp.sv_type == 3 or svp.sv_type == 4 for svp in self.SVprimitives):
            if self.inserted_sequence is None:
                raise ValueError("Inserted sequence is not set for SVpatternInsertion")
            return [pickle.loads(self.inserted_sequence)]
        return [
            seq
            for svp in self.SVprimitives
            for seq in svp.original_alt_sequences
            if svp.original_alt_sequences is not None and svp.sv_type == 0
        ]


@attrs.define
class SVpatternInvertedTranslocation(SVpattern):
    def get_sv_type(self) -> str:
        return "INVTRANS"

    def get_size(self) -> int:
        """Returns the size of the inverted translocation."""
        return abs(self.SVprimitives[-1].read_end - self.SVprimitives[0].read_start)


@attrs.define
class SVpatternTranslocation(SVpattern):
    def get_sv_type(self) -> str:
        return "TRANS"

    def get_size(self) -> int:
        return abs(self.SVprimitives[-1].read_start - self.SVprimitives[0].read_end)


@attrs.define
class SVpatternComplex(SVpattern):
    def get_sv_type(self) -> str:
        return "CPX"

    def get_size(self) -> int:
        """Returns the size of the complex SV re-arrangement."""
        return abs(self.SVprimitives[-1].read_end - self.SVprimitives[0].read_start)


@attrs.define
class SVpatternInversion(SVpattern):
    inserted_sequence: bytes | None = None
    deleted_sequence: bytes | None = None
    sequence_complexity_inserted_sequence: bytes | None = (
        None  # Pickled numpy array of float32
    )
    sequence_complexity_deleted_sequence: bytes | None = (
        None  # Pickled numpy array of float32
    )

    # needs to have an interface to retrieve the inverted sequence
    # can be a method that receives the compressed core sequence
    # and infers if it is rev-comp or needs to rev-comp it given the svprimitives.
    # cases:
    #   1) 4-relation:
    #       - middle aligned fragment (aln_is_reverse) -> correct orientation of the consensus sequence
    #       - middle aligned fragment (aln_is_reverse == False) -> needs to rev-comp the consensus sequence
    def get_sv_type(self) -> str:
        return "INV"

    def set_inserted_sequence(
        self, sequence: str, sequence_complexity_max_length: int = 300
    ) -> None:
        self.inserted_sequence = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length:
            dna_iter = iter(sequence)
            self.sequence_complexity_inserted_sequence = pickle.dumps(
                complexity_local_track(
                    dna_iter=dna_iter, w=11, K=[2, 3, 4, 5], padding=True
                )
            )
        else:
            # Create dummy placeholder with 1.0 values for long sequences
            dummy_complexity = np.ones(len(sequence), dtype=np.float16)
            self.sequence_complexity_inserted_sequence = pickle.dumps(dummy_complexity)

    def get_inserted_sequence(self) -> str | None:
        """Retrieve the inserted sequence by unpickling."""
        if self.inserted_sequence is None:
            return None
        try:
            return pickle.loads(self.inserted_sequence)
        except Exception:
            return None

    def set_deleted_sequence(
        self, sequence: str, sequence_complexity_max_length: int = 300
    ) -> None:
        self.deleted_sequence = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length:
            dna_iter = iter(sequence)
            self.sequence_complexity_deleted_sequence = pickle.dumps(
                complexity_local_track(
                    dna_iter=dna_iter, w=11, K=[2, 3, 4, 5], padding=True
                )
            )
        else:
            # Create dummy placeholder with 1.0 values for long sequences
            dummy_complexity = np.ones(len(sequence), dtype=np.float16)
            self.sequence_complexity_deleted_sequence = pickle.dumps(dummy_complexity)

    def get_deleted_sequence(self) -> str | None:
        """Retrieve the deleted sequence by unpickling."""
        if self.deleted_sequence is None:
            return None
        try:
            return pickle.loads(self.deleted_sequence)
        except Exception:
            return None

    def get_inserted_sequence_complexity(self) -> np.ndarray | None:
        """Retrieve the inserted sequence complexity scores."""
        if self.sequence_complexity_inserted_sequence is None:
            return None
        try:
            return pickle.loads(self.sequence_complexity_inserted_sequence)
        except Exception:
            return None

    def get_deleted_sequence_complexity(self) -> np.ndarray | None:
        """Retrieve the deleted sequence complexity scores."""
        if self.sequence_complexity_deleted_sequence is None:
            return None
        try:
            return pickle.loads(self.sequence_complexity_deleted_sequence)
        except Exception:
            return None

    def get_mean_inserted_complexity(self) -> float | None:
        """Get the mean complexity score for inserted sequence."""
        complexity = self.get_inserted_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.mean(complexity))
        return None

    def get_mean_deleted_complexity(self) -> float | None:
        """Get the mean complexity score for deleted sequence."""
        complexity = self.get_deleted_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.mean(complexity))
        return None

    def get_min_inserted_complexity(self) -> float | None:
        """Get the minimum complexity score for inserted sequence."""
        complexity = self.get_inserted_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.min(complexity))
        return None

    def get_min_deleted_complexity(self) -> float | None:
        """Get the minimum complexity score for deleted sequence."""
        complexity = self.get_deleted_sequence_complexity()
        if complexity is not None and len(complexity) > 0:
            return float(np.min(complexity))
        return None

    def get_sequence_from_consensus(self, consensus: Consensus) -> str:
        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        s = self.SVprimitives[0].read_start - core_sequence_start
        e = (
            self.SVprimitives[-1].read_end - core_sequence_start
        )  # this applies to bot a single insertion or a bnd->insertion
        return consensus.consensus_sequence[s:e]

    def get_reference_region(self) -> tuple[str, int, int]:
        chr = self.SVprimitives[0].chr
        start = self.SVprimitives[1].ref_start
        end = self.SVprimitives[2].ref_start
        start, end = (end, start) if start > end else (start, end)
        return (chr, start, end)

    def get_inv_sequence(self, core_sequence: bytes, core_sequence_start: int) -> str:
        """
        Returns the inverted (relative to the reference strand) sequence based on the core sequence and the SVprimitives.
        All SVprimitives must be of the same consensusID and sorted by their ref_start.
        """
        if not all(
            svpr.consensusID == self.SVprimitives[0].consensusID
            for svpr in self.SVprimitives
        ):
            raise ValueError("All SVprimitives must have the same consensusID")
        if not all(
            self.SVprimitives[i].ref_start <= self.SVprimitives[i + 1].ref_start
            for i in range(len(self.SVprimitives) - 1)
        ):
            raise ValueError("SVprimitives must be sorted by ref_start")

        seq: str = pickle.loads(core_sequence)
        if len(self.SVprimitives) == 4:
            seq = seq[
                self.SVprimitives[1].ref_start
                - core_sequence_start : self.SVprimitives[2].ref_start
                - core_sequence_start
            ]
            if self.SVprimitives[1].aln_is_reverse:
                return seq
            else:
                return str(Seq(seq).reverse_complement())
        else:
            raise ValueError("Inversion attributes must have 4 primitives")

    def get_size(self, inner: bool = True) -> int:
        """Returns the size of the inversion."""
        if len(self.SVprimitives) != 4:
            raise ValueError("Inversion attributes must have 4 primitives")
        if inner:
            return abs(self.SVprimitives[2].ref_start - self.SVprimitives[1].ref_start)
        else:
            return abs(self.SVprimitives[3].ref_end - self.SVprimitives[0].ref_start)


# Use Union for better type hints
SVpatternType = (
    SVpatternInsertion
    | SVpatternDeletion
    | SVpatternInvertedTranslocation
    | SVpatternTranslocation
    | SVpatternSingleBreakend
    | SVpatternComplex
    | SVpatternInversion
)

# endregion

# region 2-,4-relations


class TWORELATIONS(Enum):
    READCLOSED = 0  # read closed BNDs are directly adjacent without read sequence betwen them, e.g. deletions, translocations
    REFCLOSED = 1  # ref closed BNDs are directly adjacent without a gap in the reference sequence, e.g. inversions, insertions
    INVERSION = 2  # inversions are two BNDs that have different read orientations
    TRANSLOCATION = 3  # two BNDs don't share the same chr
    OVERLAP = 4  # the alignments of the two BNDs overlap on the reference
    REFGAP = (
        5  # two BNDs are separated by a gap in the reference sequence, e.g. deletions
    )
    READGAP = 6  # two BNDs are separated by a gap in the read sequence, e.g. insertions

    def __lt__(self, other):
        return self.value < other.value

    def __eq__(self, other):
        return self.value == other.value

    def __hash__(self):
        return hash(self.value)


class FOURRELATIONS(Enum):
    HOP = 0  # breakends 2,3 are outside of the reference interval 1,4 given 1,4 share the same chr.
    INVERSION = 1  # reverse aln is equal between breakends 1,4 and 2,3 bot not between 1,2 and 3,4


# two-relations allow the detection of translocations, insertions and deletions
def two_relations_of_group(
    group: list[SVprimitive], max_del_size: int, distance_tolerance: int = 30
) -> dict[tuple[int, int], set[TWORELATIONS]]:
    if len(group) < 2:
        return {}
    assert all(x.consensusID == group[0].consensusID for x in group), (
        "Not all SVprimitives have the same consensusID"
    )
    assert all(x.sv_type == 3 or x.sv_type == 4 for x in group), (
        "Not all SVprimitives are of type BND"
    )
    assert all(
        group[i].read_start <= group[i + 1].read_start for i in range(len(group) - 1)
    ), "SVprimitives are not sorted by read_start"
    assert all(x.consensus_aln_interval is not None for x in group), (
        "Not all SVprimitives have an alignment_to_ref"
    )

    two_relations: dict[tuple[int, int], set[TWORELATIONS]] = {}
    for i in range(len(group) - 1):
        a: SVprimitive = group[i]
        b: SVprimitive = group[i + 1]
        # both BNDs need to be on different alignments
        if a.alignmentID == b.alignmentID:
            continue
        tags: set[TWORELATIONS] = set()
        # READCLOSED - two adjacent bnds are no farther apart than 'distance_tolerance' bp on the read
        if abs(a.read_start - b.read_start) <= distance_tolerance:
            tags.add(TWORELATIONS.READCLOSED)
        # REFCLOSED - two adjacent bnds are no farther apart than 'distance_tolerance' bp on the reference
        if abs(a.ref_start - b.ref_start) <= distance_tolerance:
            tags.add(TWORELATIONS.REFCLOSED)
        # INVERSION - two adjacent bnds have different aln_is_reverse
        if a.aln_is_reverse != b.aln_is_reverse:
            tags.add(TWORELATIONS.INVERSION)
        # TRANSLOCATION - two adjacent bnds are not on the same chromosome
        if a.chr != b.chr or (
            a.chr == b.chr and abs(a.ref_start - b.ref_start) > max_del_size
        ):
            tags.add(TWORELATIONS.TRANSLOCATION)

        if a.chr == b.chr:
            a_it = Interval(a.consensus_aln_interval[1], a.consensus_aln_interval[2])
            b_it = Interval(b.consensus_aln_interval[1], b.consensus_aln_interval[2])
            # check if a_it really is an interval of two elements of type int
            assert isinstance(a_it.begin, int) and isinstance(a_it.end, int), (
                f"Invalid interval for {str(a.consensus_aln_interval)}"
            )
            assert isinstance(b_it.begin, int) and isinstance(b_it.end, int), (
                f"Invalid interval for {str(b.consensus_aln_interval)}"
            )
            # OVERLAP - two adjacent bnds overlap on the reference
            if a_it.overlap_size(b_it) >= distance_tolerance:
                tags.add(TWORELATIONS.OVERLAP)
            # REFGAP - two adjacent bnds are separated by a gap in the reference sequence
            elif (
                not a_it.overlaps(b_it)
                and abs(a.ref_start - b.ref_start) > distance_tolerance
            ):
                tags.add(TWORELATIONS.REFGAP)
        # READGAP - two adjacent bnds are separated by a gap in the read sequence
        if abs(a.read_start - b.read_start) > distance_tolerance:
            tags.add(TWORELATIONS.READGAP)

        two_relations[(i, i + 1)] = tags
    return two_relations


# four-relations allow the detection of complex SVs like inversions
# some insertions, e.g. of Mobile Elements, can be aligned to a different locus
# the generalization is a "hop" where a middle aigned segment between two adjacent segments is translocated.
def four_relations_of_group(
    group: list[SVprimitive], distance_tolerance: int = 5
) -> dict[tuple[int, int, int, int], set[FOURRELATIONS]]:
    if len(group) < 4:
        return {}
    assert all(x.consensusID == group[0].consensusID for x in group), (
        "Not all SVprimitives have the same consensusID"
    )
    assert all(x.sv_type == 3 or x.sv_type == 4 for x in group), (
        "Not all SVprimitives are of type BND"
    )
    assert all(
        group[i].read_start <= group[i + 1].read_start for i in range(len(group) - 1)
    ), "SVprimitives are not sorted by read_start"

    four_relations: dict[tuple[int, int, int, int], set[FOURRELATIONS]] = {}
    for i in range(len(group) - 3):
        a: SVprimitive = group[i]
        b: SVprimitive = group[i + 1]
        c: SVprimitive = group[i + 2]
        d: SVprimitive = group[i + 3]
        # breakend pairs a,b and c,d need to be on different alignments
        if a.alignmentID == b.alignmentID or c.alignmentID == d.alignmentID:
            continue
        tags: set[FOURRELATIONS] = set()
        # HOP - breakends b,c are outside of the reference interval a,d given a,d share the same chr.
        if (
            a.alignmentID != b.alignmentID
            and b.alignmentID == c.alignmentID
            and c.alignmentID != d.alignmentID
        ):
            interval_ad = Interval(a.ref_start, d.ref_start)
            if a.chr != b.chr:
                tags.add(FOURRELATIONS.HOP)
            elif not interval_ad.overlaps(
                b.ref_start - distance_tolerance, b.ref_start + distance_tolerance
            ) and not interval_ad.overlaps(
                c.ref_start - distance_tolerance, c.ref_start + distance_tolerance
            ):
                tags.add(FOURRELATIONS.HOP)
        # INVERSION - reverse aln is equal between breakends a,d and b,c but not between a,b and c,d
        if (
            a.aln_is_reverse == d.aln_is_reverse
            and b.aln_is_reverse == c.aln_is_reverse
            and a.aln_is_reverse != b.aln_is_reverse
        ):
            tags.add(FOURRELATIONS.INVERSION)

        four_relations[(i, i + 1, i + 2, i + 3)] = tags
    return four_relations


def possible_inversions_from_BNDs(
    svps: list[SVprimitive],
    fourrelations: dict[tuple[int, int, int, int], set[FOURRELATIONS]],
) -> list[tuple[int, int, int, int]]:
    """Calls inversions from a group of SVprimitives and their four-relations. The returned list of 4-tuples contains the indices of the SVprimitives that can form an inversion."""
    if len(svps) < 4:
        raise ValueError(
            "At least 4 SVprimitives are required to call inversions"
        )  # or just return an empty list?
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError(
            "All SVprimitives must be of type BND to call inter-alignment inversions"
        )
    # an inversion is present if there exists a 4-hop and a 4-inv
    results: list[tuple[int, int, int, int]] = []
    for (a, b, c, d), x in fourrelations.items():
        if x == set([FOURRELATIONS.HOP, FOURRELATIONS.INVERSION]):
            # create a SVcomplex from the group
            results.append((a, b, c, d))
    return results


def possible_insertions_from_BNDs(
    svps: list[SVprimitive], tworelations: dict[tuple[int, int], set[TWORELATIONS]]
) -> list[tuple[int, int]]:
    """Calls insertions from a group of SVprimitives and their two-relations. The returned list of 2-tuples contains the indices of the SVprimitives that can form an insertion."""
    # an insertion is present if there exists a 2-refclosed and a 2-readgap
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError(
            "All SVprimitives must be of type BND to call inter-alignment insertions"
        )
    if len(svps) < 2:
        return []
    results: list[tuple[int, int]] = []
    for (a, b), x in tworelations.items():
        if x == set([TWORELATIONS.REFCLOSED, TWORELATIONS.READGAP]):
            results.append((a, b))
    return results


def possible_deletions_from_BNDs(
    svps: list[SVprimitive], tworelations: dict[tuple[int, int], set[TWORELATIONS]]
) -> list[tuple[int, int]]:
    """Calls deletions from a group of SVprimitives and their two-relations. The returned list of 2-tuples contains the indices of the SVprimitives that can form a deletion."""
    # a deletion is present if there exists a 2-readclosed and a 2-REFGAP
    # check if all sv types are BND
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError(
            "All SVprimitives must be of type BND to call inter-alignment deletions"
        )
    if len(svps) < 2:
        return []
    results: list[tuple[int, int]] = []
    for (a, b), x in tworelations.items():
        if x == set([TWORELATIONS.READCLOSED, TWORELATIONS.REFGAP]):
            results.append((a, b))
    return results


def possible_inverted_translocations_from_BNDs(
    svps: list[SVprimitive], tworelations: dict[tuple[int, int], set[TWORELATIONS]]
) -> list[tuple[int, int]]:
    """calls translocations from a group of SVprimitives. The returned list of 2-tuples contains the indices of the SVprimitives that can form a translocation."""
    # a translocation is present if there exists a 2-readclosed or a 2-refgap
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError("All SVprimitives must be of type BND to call translocations")
    if len(svps) < 1:
        return []
    results: list[tuple[int, int]] = []
    for (
        a,
        b,
    ), x in tworelations.items():  # ignores readgap or refgap (more complex SVs)
        if TWORELATIONS.TRANSLOCATION in x and TWORELATIONS.INVERSION in x:
            results.append((a, b))
    return results


def possible_translocations_from_BNDs(
    svps: list[SVprimitive], tworelations: dict[tuple[int, int], set[TWORELATIONS]]
) -> list[tuple[int, int]]:
    """calls translocations from a group of SVprimitives. The returned list of 2-tuples contains the indices of the SVprimitives that can form a translocation."""
    # a translocation is present if there exists a 2-readclosed or a 2-refgap
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError("All SVprimitives must be of type BND to call translocations")
    if len(svps) < 1:
        return []
    results: list[tuple[int, int]] = []
    for (
        a,
        b,
    ), x in tworelations.items():  # ignores readgap or refgap (more complex SVs)
        if TWORELATIONS.TRANSLOCATION in x:
            results.append((a, b))
    return results


def possible_single_ended_breakends_from_BNDs(
    svps: list[SVprimitive],
) -> list[tuple[int]]:
    """calls single-ended breakends from a group of SVprimitives. The returned list of 1-tuples contains the indices of the SVprimitives that can form a single-ended breakend."""
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError(
            "All SVprimitives must be of type BND to call single-ended breakends"
        )
    if len(svps) < 1:
        return []
    results: list[tuple[int]] = []
    for i, svp in enumerate(svps):
        # check if the SVprimitive is a single-ended breakend (no adjacency)
        if svp.adjacent_bnd is None:
            results.append((i,))
    return results


# endregion


def parse_SVprimitives_to_SVpatterns(
    SVprimitives: list[SVprimitive], max_del_size: int = 100_000
) -> list[SVpatternType]:
    """
    Parses Sv patterns from the SVprimitives of one consensus. All SVprimitives must have the same consensusID.
    """
    if not all(svp.consensusID == SVprimitives[0].consensusID for svp in SVprimitives):
        raise ValueError(
            "All SVprimitives must have the same consensusID to parse SVpatterns"
        )

    # parse primitives to simple SV types (INS, DEL, BND)
    result: list[SVpatternType] = []

    # at first, separate INDELS from breakends
    indels = [svp for svp in SVprimitives if svp.sv_type <= 2]
    breakends = [svp for svp in SVprimitives if svp.sv_type > 2]

    used_indices: set[int] = set()

    for svp in indels:
        if svp.sv_type == 0:
            result.append(SVpatternInsertion(SVprimitives=[svp]))
        else:
            result.append(SVpatternDeletion(SVprimitives=[svp]))

    if len(breakends) == 0:
        return result

    if len(breakends) == 4:
        fourrelations = four_relations_of_group(group=breakends)
        indices_inversions: list[tuple[int, int, int, int]] = (
            possible_inversions_from_BNDs(svps=breakends, fourrelations=fourrelations)
        )
        for a, b, c, d in indices_inversions:
            if (
                a in used_indices
                or b in used_indices
                or c in used_indices
                or d in used_indices
            ):
                continue
            result.append(
                SVpatternInversion(
                    SVprimitives=[
                        breakends[a],
                        breakends[b],
                        breakends[c],
                        breakends[d],
                    ]
                )
            )
            used_indices.update([a, b, c, d])

    if len(breakends) >= 2 and len(used_indices) < len(breakends):
        tworelations = two_relations_of_group(
            group=breakends, max_del_size=max_del_size
        )

        idx_tuples_insertions: list[tuple[int, int]] = possible_insertions_from_BNDs(
            svps=breakends, tworelations=tworelations
        )
        for a, b in idx_tuples_insertions:
            if a in used_indices or b in used_indices:
                continue
            used_indices.update([a, b])
            result.append(SVpatternInsertion(SVprimitives=[breakends[a], breakends[b]]))

        idx_tuples_deletions: list[tuple[int, int]] = possible_deletions_from_BNDs(
            svps=breakends, tworelations=tworelations
        )
        for a, b in idx_tuples_deletions:
            if a in used_indices or b in used_indices:
                continue
            used_indices.update([a, b])
            result.append(SVpatternDeletion(SVprimitives=[breakends[a], breakends[b]]))

        idx_tuples_inverted_translocations: list[tuple[int, int]] = (
            possible_inverted_translocations_from_BNDs(
                svps=breakends, tworelations=tworelations
            )
        )
        for a, b in idx_tuples_inverted_translocations:
            if a in used_indices or b in used_indices:
                continue
            used_indices.update([a, b])
            result.append(
                SVpatternInvertedTranslocation(
                    SVprimitives=[breakends[a], breakends[b]]
                )
            )

        idx_tuples_translocations: list[tuple[int, int]] = (
            possible_translocations_from_BNDs(svps=breakends, tworelations=tworelations)
        )
        for a, b in idx_tuples_translocations:
            if a in used_indices or b in used_indices:
                continue
            used_indices.update([a, b])
            result.append(
                SVpatternTranslocation(SVprimitives=[breakends[a], breakends[b]])
            )

    # if there is only one breakend left, it is a single-ended breakend
    if len(used_indices) - len(breakends) == 1:
        unused_index = (set(range(len(breakends))) - used_indices).pop()
        result.append(SVpatternSingleBreakend(SVprimitives=[breakends[unused_index]]))
        used_indices.add(unused_index)

    # if there are still unused breakends, create a complex SVpattern
    unused_indices = set(range(len(breakends))) - used_indices
    if unused_indices:
        remaining_breakends = [breakends[i] for i in unused_indices]
        result.append(SVpatternComplex(SVprimitives=remaining_breakends))

    return result


def distortions_by_svPattern(
    svPattern: SVpatternType,
    consensus: Consensus,
    distance_scale: float,
    falloff: float,
) -> dict[str, float]:
    r"""Return distortion magnitudes weighted by their distance to the SV pattern.

    The weighting uses exponential decay and returns weighted means per read name."""

    supporting_reads = svPattern.get_supporting_reads()
    result = dict.fromkeys(supporting_reads, 0.0)
    distortions = [
        d
        for d in consensus.get_consensus_distortions()
        if d.readname in supporting_reads
    ]
    if not distortions:
        log.debug(
            f"No distortions found for consensus {consensus.ID}. Returning trivial result for {len(supporting_reads)} reads."
        )
    if len(distortions) == 0:
        return result

    # Calculate SV pattern boundaries
    sv_start = svPattern.SVprimitives[0].ref_start
    sv_end = svPattern.SVprimitives[-1].ref_end

    # Group distortions by read name
    distortions_by_read: dict[str, list[tuple[float, float]]] = {}
    for distortion in distortions:
        distance = min(
            abs(distortion.position - sv_start), abs(distortion.position - sv_end)
        )
        weight = exponential_weight(
            distance=distance, scale=distance_scale, falloff=falloff
        )

        if distortion.readname not in distortions_by_read:
            distortions_by_read[distortion.readname] = []
        distortions_by_read[distortion.readname].append((distortion.size, weight))

    # Calculate weighted mean for each read
    weighted_means: dict[str, float] = dict.fromkeys(
        svPattern.get_supporting_reads(), 0.0
    )
    for readname, size_weight_pairs in distortions_by_read.items():
        total_weighted_size = sum(size * weight for size, weight in size_weight_pairs)
        total_weight = sum(weight for _, weight in size_weight_pairs)

        if total_weight > 0:
            weighted_means[readname] = total_weighted_size / total_weight
        else:
            weighted_means[readname] = 0.0

    return weighted_means


# def build_size_population_by_svPattern(
#         base_size: int,
#         svPattern: SVpatternType) -> list[int]:
#     if not svPattern.size_distortions:
#         raise ValueError("SVpattern has no size distortions defined.")
#     return [base_size + value for value in svPattern.size_distortions.values()]


# ======== CATTRS CONFIGURATION FOR SERIALIZATION ======== #

# Configure cattrs to handle the SVpattern union
converter = cattrs.Converter()


# Custom hooks to handle pickled bytes data for JSON serialization
def unstructure_bytes_field(field_value):
    """Convert pickled bytes to base64 string for JSON serialization."""
    if field_value is None:
        return None
    try:
        # Try to unpickle - if it's a numpy array, convert to list
        unpickled = pickle.loads(field_value)
        if isinstance(unpickled, np.ndarray):
            return unpickled.tolist()
        return unpickled
    except Exception as e:
        # If unpickling fails, encode as base64
        import base64

        log.error(f"Failed to unpickle bytes field: {e}")
        return base64.b64encode(field_value).decode("utf-8")


def structure_bytes_field(field_value, _):
    """Convert data back to pickled bytes for object reconstruction."""
    if field_value is None:
        return None
    try:
        # If it's a list (from numpy array), convert back to numpy array and pickle
        if isinstance(field_value, list):
            array = np.array(
                field_value, dtype=np.float16
            )  # Changed from float32 to float16
            return pickle.dumps(array)
        # If it's a string, pickle it directly (it's the unpickled sequence)
        elif isinstance(field_value, str):
            return pickle.dumps(field_value)
        else:
            return pickle.dumps(field_value)
    except Exception:
        return None


# Register hooks for bytes fields in SVpattern classes
converter.register_unstructure_hook(bytes, unstructure_bytes_field)

converter.register_structure_hook(bytes, structure_bytes_field)


# Register unstructure hooks for each SVpattern type
def _unstructure_svpattern_insertion(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternInsertion", "data": base_converter.unstructure(obj)}


def _unstructure_svpattern_deletion(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternDeletion", "data": base_converter.unstructure(obj)}


def _unstructure_svpattern_inversion(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternInversion", "data": base_converter.unstructure(obj)}


def _unstructure_svpattern_translocation(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternTranslocation", "data": base_converter.unstructure(obj)}


def _unstructure_svpattern_inverted_translocation(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {
        "type": "SVpatternInvertedTranslocation",
        "data": base_converter.unstructure(obj),
    }


def _unstructure_svpattern_single_breakend(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternSingleBreakend", "data": base_converter.unstructure(obj)}


def _unstructure_svpattern_complex(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternComplex", "data": base_converter.unstructure(obj)}


converter.register_unstructure_hook(
    SVpatternInsertion, _unstructure_svpattern_insertion
)
converter.register_unstructure_hook(SVpatternDeletion, _unstructure_svpattern_deletion)
converter.register_unstructure_hook(
    SVpatternInversion, _unstructure_svpattern_inversion
)
converter.register_unstructure_hook(
    SVpatternTranslocation, _unstructure_svpattern_translocation
)
converter.register_unstructure_hook(
    SVpatternInvertedTranslocation, _unstructure_svpattern_inverted_translocation
)
converter.register_unstructure_hook(
    SVpatternSingleBreakend, _unstructure_svpattern_single_breakend
)
converter.register_unstructure_hook(SVpatternComplex, _unstructure_svpattern_complex)


# Register structure hooks
def structure_svpattern(data, _):
    pattern_type = data["type"]
    pattern_data = data["data"]

    # Create a base converter for structuring with bytes hook
    base_converter = cattrs.Converter()
    base_converter.register_structure_hook(bytes, structure_bytes_field)

    if pattern_type == "SVpatternInsertion":
        return base_converter.structure(pattern_data, SVpatternInsertion)
    elif pattern_type == "SVpatternDeletion":
        return base_converter.structure(pattern_data, SVpatternDeletion)
    elif pattern_type == "SVpatternInversion":
        return base_converter.structure(pattern_data, SVpatternInversion)
    elif pattern_type == "SVpatternTranslocation":
        return base_converter.structure(pattern_data, SVpatternTranslocation)
    elif pattern_type == "SVpatternInvertedTranslocation":
        return base_converter.structure(pattern_data, SVpatternInvertedTranslocation)
    elif pattern_type == "SVpatternSingleBreakend":
        return base_converter.structure(pattern_data, SVpatternSingleBreakend)
    elif pattern_type == "SVpatternComplex":
        return base_converter.structure(pattern_data, SVpatternComplex)
    else:
        raise ValueError(f"Unknown SVpattern type: {pattern_type}")


converter.register_structure_hook(SVpatternType, structure_svpattern)


# ======== DATABASE FUNCTIONS ======== #


def create_svPatterns_db(database: Path):
    """Create database for storing SVpatterns."""
    with sqlite3.connect("file:" + str(database) + "?mode=rwc", uri=True) as conn:
        conn.execute("DROP TABLE IF EXISTS svPatterns")
        conn.execute(
            """CREATE TABLE svPatterns (
                        svPatternID VARCHAR(50) PRIMARY KEY,
                        crID INTEGER,
                        consensusID VARCHAR(30),
                        svPattern BLOB)"""
        )
        conn.execute(
            "CREATE INDEX idx_svPatterns_consensusID ON svPatterns (consensusID)"
        )
        conn.execute("CREATE INDEX idx_svPatterns_crID ON svPatterns (crID)")
        conn.commit()


def write_svPatterns_to_db(
    database: Path,
    data: list[SVpatternType],
    timeout: float = 10.0,
    batch_size: int = 1000,
):
    """Write SVpatterns to database in batches to avoid memory issues.

    Args:
        database: Path to the SQLite database
        data: List of SVpatterns to write
        timeout: SQLite connection timeout
        batch_size: Number of records to write per batch (default: 1000)
    """
    query = "INSERT INTO svPatterns (svPatternID, consensusID, crID, svPattern) VALUES (?, ?, ?, ?)"

    with sqlite3.connect(
        "file:" + str(database) + "?mode=rwc", uri=True, timeout=timeout
    ) as conn:
        c = conn.cursor()
        cache = []

        for i, svPattern in enumerate(data):
            # Generate a unique ID for the SVpattern
            if svPattern.SVprimitives:
                consensusID = svPattern.SVprimitives[0].consensusID
                crID = int(consensusID.split(".")[0])
                svPatternID = f"{svPattern.get_sv_type()}.{consensusID}.{i}"
            else:
                raise ValueError(f"SVpattern has no primitives: {svPattern}")

            cache.append((
                svPatternID,
                consensusID,
                crID,
                pickle.dumps(converter.unstructure(svPattern)),
            ))

            # Write batch when cache reaches batch_size
            if len(cache) >= batch_size:
                c.executemany(query, cache)
                conn.commit()
                cache.clear()

        # Write remaining records
        if cache:
            c.executemany(query, cache)
            conn.commit()

        c.close()


def read_svPatterns_from_db(
    database: Path, consensusIDs: list[str] | None = None, crIDs: set[int] | None = None
) -> list[SVpatternType]:
    """Read SVpatterns from database."""
    with sqlite3.connect(str(database), timeout=30.0) as conn:
        svPatterns: list[SVpatternType] = []
        c = conn.cursor()
        try:
            if crIDs:
                placeholders = ",".join(["?" for _ in crIDs])
                query = (
                    f"SELECT svPattern FROM svPatterns WHERE crID IN ({placeholders})"
                )
                log.info(
                    f"reading svPatterns from {str(database)} for {len(crIDs)} crIDs"
                )
                c.execute(
                    query,
                    [
                        *crIDs,
                    ],
                )
            elif consensusIDs:
                placeholders = ",".join(["?" for _ in consensusIDs])
                query = f"SELECT svPattern FROM svPatterns WHERE consensusID IN ({placeholders})"
                log.info(
                    f"reading svPatterns from {str(database)} for {len(consensusIDs)} consensusIDs"
                )
                c.execute(
                    query,
                    [
                        *consensusIDs,
                    ],
                )
            else:
                query = "SELECT svPattern FROM svPatterns"
                log.info(f"reading all svPatterns from {str(database)}")
                c.execute(query)

            for row in c.fetchall():
                # Deserialize using the custom converter
                pattern_data = pickle.loads(row[0])
                svPatterns.append(converter.structure(pattern_data, SVpatternType))
        except sqlite3.Error as e:
            log.error(f"Error reading svPatterns from {str(database)}: {e}")
            # print all tables found in the database
            c.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = c.fetchall()
            log.error(f"Tables in database: {tables}")
            c.close()
            raise e
        c.close()

    return svPatterns
