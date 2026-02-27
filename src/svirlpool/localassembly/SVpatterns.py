from __future__ import annotations

import json
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

from ..util.util import complexity_local_track, exponential_weight
from .consensus_class import Consensus
from .SVprimitives import SVprimitive

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

    @classmethod
    @abstractmethod
    def get_sv_type(cls) -> str:
        pass

    @abstractmethod
    def get_size(self) -> int:
        pass

    @abstractmethod
    def get_reference_region(self) -> tuple[str, int, int]:
        pass

    @abstractmethod
    def _log_id(self) -> str:
        """Returns a string identifier for logging purposes, e.g. samplename:consensusIDs:type:size:regions"""
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
        return json.dumps(converter.unstructure(self), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> SVpatternType:
        """Create SVpattern from JSON string."""
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


def _get_first_svp_on_reference_start_pos(svp: SVpattern) -> tuple[str, int, int]:
    # pick first of SVpatterns after sorting them by chr, start
    _svps = sorted(svp.SVprimitives, key=lambda svp: (svp.chr, svp.ref_start))
    chr = _svps[0].chr
    start = _svps[0].ref_start
    end = _svps[-1].ref_end
    return (chr, min(start, end), max(start, end))


@attrs.define
class SVpatternSingleBreakend(SVpattern):
    CONTEXT_SIZE: int = 200
    sequence_context: bytes | None = (
        None  # the sequence on the consensus +/- 200 pb around the break end
    )
    clipped_tail: bytes | None = (
        None  # the clipped tail sequence on the consenus sequence
    )
    sequence_complexity: bytes | None = None  # Pickled numpy array of float32

    @classmethod
    def get_sv_type(cls) -> str:
        return "BND"

    def get_size(self) -> int:
        return 0

    def get_reference_region(self) -> tuple[str, int, int]:
        chr = self.SVprimitives[0].chr
        start = self.SVprimitives[0].ref_start
        end = self.SVprimitives[0].ref_end
        return (chr, min(start, end), max(start, end))

    def get_sequence(self) -> str | None:
        """returns the context by default. Use get_sequence_clipped() to get the clipped tail."""
        return self.get_sequence_context()

    def get_sequence_clipped(self) -> str | None:
        """Retrieve the clipped sequence by unpickling."""
        if self.clipped_tail is None:
            return None
        else:
            try:
                return pickle.loads(self.clipped_tail)
            except Exception:
                return None

    def get_sequence_context(self) -> str | None:
        """Retrieve the clipped sequence by unpickling."""
        if self.sequence_context is None:
            return None
        else:
            try:
                return pickle.loads(self.sequence_context)
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

    def get_context_sequence_from_consensus(self, consensus: Consensus) -> str:
        """Extract context region around the breakpoint (±CONTEXT_SIZE bp)."""
        if consensus.consensus_padding is None:
            raise ValueError("Consensus sequence is not set in the Consensus object")

        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        svp = self.SVprimitives[0]  # Single breakend has exactly one primitive
        consensus_length = len(consensus.consensus_sequence)

        # Convert read coordinates to consensus coordinates
        read_start_on_consensus = svp.read_start - core_sequence_start
        read_end_on_consensus = svp.read_end - core_sequence_start

        # Extract context region around the breakpoint
        s = max(0, read_start_on_consensus - self.CONTEXT_SIZE)
        e = min(consensus_length, read_end_on_consensus + self.CONTEXT_SIZE)

        seq = consensus.consensus_sequence[s:e]

        # Apply reverse complement if alignment is reversed
        if svp.aln_is_reverse:
            seq = str(Seq(seq).reverse_complement())

        return seq

    def set_context_sequence(
        self, sequence: str, sequence_complexity_max_length: int = CONTEXT_SIZE
    ) -> None:
        """Set the inserted sequence and compute its complexity."""
        self.sequence_context = pickle.dumps(sequence)
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

    def get_clipped_sequence_from_consensus(self, consensus: Consensus) -> str:
        """Extract the clipped tail sequence based on sv_type and aln_is_reverse."""
        if consensus.consensus_padding is None:
            raise ValueError("Consensus sequence is not set in the Consensus object")

        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        svp = self.SVprimitives[0]  # Single breakend has exactly one primitive
        consensus_length = len(consensus.consensus_sequence)

        # Convert read coordinates to consensus coordinates
        read_start_on_consensus = svp.read_start - core_sequence_start
        read_end_on_consensus = svp.read_end - core_sequence_start

        # Extract the clipped tail based on sv_type and aln_is_reverse
        # Confusion matrix:
        # sv_type=3 (left), aln_is_reverse=False  → left side clipped  → [0 : read_start]
        # sv_type=3 (left), aln_is_reverse=True   → right side clipped → [read_end : end]
        # sv_type=4 (right), aln_is_reverse=False → right side clipped → [read_end : end]
        # sv_type=4 (right), aln_is_reverse=True  → left side clipped  → [0 : read_start]

        if svp.sv_type == 3:  # Left breakend on reference
            if not svp.aln_is_reverse:
                # Forward alignment: left on ref = left on consensus
                s = 0
                e = read_start_on_consensus
            else:
                # Reverse alignment: left on ref = right on consensus
                s = read_end_on_consensus
                e = consensus_length
        elif svp.sv_type == 4:  # Right breakend on reference
            if not svp.aln_is_reverse:
                # Forward alignment: right on ref = right on consensus
                s = read_end_on_consensus
                e = consensus_length
            else:
                # Reverse alignment: right on ref = left on consensus
                s = 0
                e = read_start_on_consensus
        else:
            raise ValueError(
                f"SVpatternSingleBreakend must have sv_type 3 or 4, got {svp.sv_type}"
            )

        # Ensure valid coordinates
        s = max(0, s)
        e = min(consensus_length, e)

        seq = consensus.consensus_sequence[s:e]

        # Apply reverse complement if alignment is reversed
        if svp.aln_is_reverse:
            seq = str(Seq(seq).reverse_complement())

        return seq

    def set_clipped_sequence(self, sequence: str) -> None:
        """Set the clipped sequence."""
        self.clipped_tail = pickle.dumps(sequence)

    def _log_id(self) -> str:
        chr, start, end = self.get_reference_region()
        return f"sample={self.samplename}|consensusID={self.consensusID}|crID={self.consensusID.split('.')[0]}|type={self.get_sv_type()}|size={self.get_size()}|region={chr}:{start}-{end}"


@attrs.define
class SVpatternDeletion(SVpattern):
    deleted_sequence: bytes | None = None
    sequence_complexity: bytes | None = None  # Pickled numpy array of float32

    @classmethod
    def get_sv_type(cls) -> str:
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

    def _log_id(self) -> str:
        chr, start, end = self.get_reference_region()
        return f"sample={self.samplename}|consensusID={self.consensusID}|crID={self.consensusID.split('.')[0]}|type={self.get_sv_type()}|size={self.get_size()}|region={chr}:{start}-{end}"


@attrs.define
class SVpatternInsertion(SVpattern):
    inserted_sequence: bytes | None = None
    sequence_complexity: bytes | None = None  # Pickled numpy array of float32

    @classmethod
    def get_sv_type(cls) -> str:
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
        if consensus.consensus_padding is None:
            raise ValueError("Consensus sequence is not set in the Consensus object")
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

    def _log_id(self) -> str:
        chr, start, end = self.get_reference_region()
        return f"sample={self.samplename}|consensusID={self.consensusID}|crID={self.consensusID.split('.')[0]}|type={self.get_sv_type()}|size={self.get_size()}|region={chr}:{start}-{end}"


@attrs.define
class SVpatternInversion(SVpattern):
    inserted_sequence: bytes | None = None
    deleted_sequence: bytes | None = None
    sequence_complexity: bytes | None = None  # Pickled numpy array of float32

    @classmethod
    def get_sv_type(cls) -> str:
        return "INV"

    def set_sequence(
        self,
        sequence: str,
        sequence_complexity_max_length: int = 300,
        write_complexity: bool = True,
    ) -> None:
        """Set the inserted sequence and compute its complexity."""
        self.inserted_sequence = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length and write_complexity:
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

    def set_deleted_sequence(
        self,
        sequence: str,
        sequence_complexity_max_length: int = 300,
        write_complexity: bool = True,
    ) -> None:
        """Set the deleted sequence and compute its complexity."""
        self.deleted_sequence = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length and write_complexity:
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

    def get_deleted_sequence(self) -> str | None:
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

    def get_sequence_from_consensus(self, consensus: Consensus) -> str:
        if consensus.consensus_padding is None:
            raise ValueError("Consensus sequence is not set in the Consensus object")
        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        s = self.SVprimitives[0].read_start - core_sequence_start
        e = (
            self.SVprimitives[-1].read_end - core_sequence_start
        )  # this applies to bot a single insertion or a bnd->insertion
        return consensus.consensus_sequence[s:e]

    # alway return the inner region of the inversion
    def get_reference_region(self, inner: bool = False) -> tuple[str, int, int]:
        idx = (1, 2) if inner else (0, 3)
        chr = self.SVprimitives[idx[0]].chr
        start = self.SVprimitives[idx[0]].ref_start
        end = self.SVprimitives[idx[1]].ref_start
        start, end = (end, start) if start > end else (start, end)
        return (chr, start, end)

    def get_size(self, inner: bool = True) -> int:
        """Returns the size of the inversion."""
        if len(self.SVprimitives) != 4:
            raise ValueError("Inversion attributes must have 4 primitives")
        if inner:
            return abs(self.SVprimitives[2].ref_start - self.SVprimitives[1].ref_start)
        else:
            return abs(self.SVprimitives[3].ref_end - self.SVprimitives[0].ref_start)

    def _log_id(self) -> str:
        chr, start, end = self.get_reference_region()
        return f"sample={self.samplename}|consensusID={self.consensusID}|crID={self.consensusID.split('.')[0]}|type={self.get_sv_type()}|size={self.get_size()}|region={chr}:{start}-{end}"


@attrs.define
class SVpatternInversionDeletion(SVpatternInversion):
    @classmethod
    def get_sv_type(cls) -> str:
        return "INV-DEL"

    def get_size(self, inner: bool = False) -> int:
        """Returns the size of the inverted deletion."""
        return super().get_size(inner=inner)

    def get_reference_region(self, inner: bool = False) -> tuple[str, int, int]:
        return super().get_reference_region(inner=inner)

    def _log_id(self) -> str:
        chr, start, end = self.get_reference_region()
        return f"sample={self.samplename}|consensusID={self.consensusID}|crID={self.consensusID.split('.')[0]}|type={self.get_sv_type()}|size={self.get_size()}|region={chr}:{start}-{end}"


@attrs.define
class SVpatternInversionDuplication(SVpatternInversion):
    @classmethod
    def get_sv_type(cls) -> str:
        return "INV-DUP"

    def get_size(self, inner: bool = True) -> int:
        """Returns the size of the inverted deletion."""
        return super().get_size(inner=inner)

    def get_reference_region(self, inner: bool = True) -> tuple[str, int, int]:
        return super().get_reference_region(inner=inner)

    def _log_id(self) -> str:
        chr, start, end = self.get_reference_region()
        return f"sample={self.samplename}|consensusID={self.consensusID}|crID={self.consensusID.split('.')[0]}|type={self.get_sv_type()}|size={self.get_size()}|region={chr}:{start}-{end}"


@attrs.define
class SVpatternAdjacency(SVpattern):
    """A pattern that defines an adjacency defined by two connected break ends."""

    CONTEXT_SIZE: int = 200
    sequence_contexts: dict[int, bytes] | None = None
    """A dictionary mapping SVprimitive indices to their sequence contexts (pickled byte strings)."""
    sequence_contexts_complexities: dict[int, bytes] | None = None
    """A dictionary mapping SVprimitive indices to their sequence complexities (pickled numpy arrays)."""
    inserted_sequence: bytes | None = None
    """The inserted sequence at the adjacency (pickled byte string)."""
    sequence_complexity: bytes | None = None
    """Pickled numpy array of float32 representing the complexity of the inserted sequence."""

    def __attrs_post_init__(self):
        super().__attrs_post_init__()
        if len(self.SVprimitives) != 2:
            raise ValueError(
                f"SVpatternAdjacency must contain exactly two SVprimitives, got {len(self.SVprimitives)}"
            )

    @classmethod
    def get_sv_type(cls) -> str:
        return "BND"

    def get_size(self) -> int:
        """Returns the size of the adjacency on the reference."""
        return abs(self.SVprimitives[-1].ref_end - self.SVprimitives[0].ref_start)

    def get_reference_region(self) -> tuple[str, int, int]:
        return _get_first_svp_on_reference_start_pos(self)

    def get_reference_regions(self) -> list[tuple[str, int, int]]:
        """Gets all reference regions from all SVprimitives."""
        if len(self.SVprimitives) != 2:
            raise ValueError("Adjacency must have exactly two SVprimitives.")
        return [
            (
                self.SVprimitives[0].chr,
                self.SVprimitives[0].ref_start,
                self.SVprimitives[0].ref_end,
            ),
            (
                self.SVprimitives[1].chr,
                self.SVprimitives[1].ref_start,
                self.SVprimitives[1].ref_end,
            ),
        ]

    def set_sequence_context(
        self,
        svp_index: int,
        sequence: str,
        sequence_complexity_max_length: int = 300,
    ) -> None:
        """Set the sequence context and compute its complexity for a specific SVprimitive.

        Args:
            svp_index: Index of the SVprimitive (must be 0 or 1)
            sequence: The sequence context to set
            sequence_complexity_max_length: Maximum length for complexity computation
        """
        if svp_index not in (0, 1):
            raise ValueError(f"svp_index must be 0 or 1, got {svp_index}")

        if self.sequence_contexts is None:
            self.sequence_contexts = {}
        if self.sequence_contexts_complexities is None:
            self.sequence_contexts_complexities = {}

        self.sequence_contexts[svp_index] = pickle.dumps(sequence)
        if len(sequence) <= sequence_complexity_max_length:
            dna_iter = iter(sequence)
            self.sequence_contexts_complexities[svp_index] = pickle.dumps(
                complexity_local_track(
                    dna_iter=dna_iter, w=11, K=[2, 3, 4, 5], padding=True
                )
            )
        else:
            # Create dummy placeholder with 1.0 values for long sequences
            dummy_complexity = np.ones(len(sequence), dtype=np.float16)
            self.sequence_contexts_complexities[svp_index] = pickle.dumps(
                dummy_complexity
            )

    def set_all_sequence_contexts(
        self,
        contexts: dict[int, str],
        sequence_complexity_max_length: int = CONTEXT_SIZE,
    ) -> None:
        """Set the sequence contexts and compute their complexities for both SVprimitives.

        Args:
            contexts: Dictionary mapping indices (0 and 1) to sequence contexts
            sequence_complexity_max_length: Maximum length for complexity computation
        """
        if not contexts:
            return

        # Validate that we have contexts for indices 0 and 1 only
        invalid_indices = set(contexts.keys()) - {0, 1}
        if invalid_indices:
            raise ValueError(
                f"contexts dict contains invalid indices: {invalid_indices}. Only 0 and 1 are allowed."
            )

        for svp_index in (0, 1):
            if svp_index in contexts:
                self.set_sequence_context(
                    svp_index,
                    contexts[svp_index],
                    sequence_complexity_max_length=sequence_complexity_max_length,
                )

    def get_sequence_context(self, svp_index: int) -> str | None:
        """Retrieve the sequence context for a specific SVprimitive by unpickling.

        Args:
            svp_index: Index of the SVprimitive (must be 0 or 1)

        Returns:
            The sequence context string, or None if not set
        """
        if svp_index not in (0, 1):
            raise ValueError(f"svp_index must be 0 or 1, got {svp_index}")

        if self.sequence_contexts is None or svp_index not in self.sequence_contexts:
            return None
        try:
            return pickle.loads(self.sequence_contexts[svp_index])
        except Exception:
            return None

    def get_all_sequence_contexts(self) -> dict[int, str]:
        """Retrieve both sequence contexts by unpickling.

        Returns:
            Dictionary mapping indices (0, 1) to sequence context strings.
            Missing or unpickleable contexts are set to None.
        """
        if self.sequence_contexts is None:
            return {}

        contexts = {}
        # Process both indices explicitly (0 and 1)
        for svp_index in (0, 1):
            if svp_index in self.sequence_contexts:
                try:
                    contexts[svp_index] = pickle.loads(
                        self.sequence_contexts[svp_index]
                    )
                except Exception:
                    contexts[svp_index] = None
        return contexts

    def get_sequence_contexts_from_consensus(
        self, consensus: Consensus
    ) -> dict[int, str]:
        """Extract context regions around both breakpoints (±CONTEXT_SIZE bp).

        Returns:
            Dictionary with keys 0 and 1 mapping to sequence contexts for each SVprimitive
        """
        if consensus.consensus_padding is None:
            raise ValueError("Consensus sequence is not set in the Consensus object")

        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        consensus_length = len(consensus.consensus_sequence)
        contexts = {}

        # Process both SVprimitives (indices 0 and 1)
        svp_0, svp_1 = self.SVprimitives[0], self.SVprimitives[1]

        for idx, svp in [(0, svp_0), (1, svp_1)]:
            # Convert read coordinates to consensus coordinates
            read_start_on_consensus = svp.read_start - core_sequence_start
            read_end_on_consensus = svp.read_end - core_sequence_start

            # Extract context region around the breakpoint
            s = max(0, read_start_on_consensus - self.CONTEXT_SIZE)
            e = min(consensus_length, read_end_on_consensus + self.CONTEXT_SIZE)

            seq = consensus.consensus_sequence[s:e]

            # Apply reverse complement if alignment is reversed
            if svp.aln_is_reverse:
                seq = str(Seq(seq).reverse_complement())

            contexts[idx] = seq

        return contexts

    def get_sequence_context_complexity(self, svp_index: int) -> np.ndarray | None:
        """Retrieve the sequence context complexity for a specific SVprimitive by unpickling.

        Args:
            svp_index: Index of the SVprimitive (must be 0 or 1)

        Returns:
            Numpy array of complexity scores, or None if not set
        """
        if svp_index not in (0, 1):
            raise ValueError(f"svp_index must be 0 or 1, got {svp_index}")

        if (
            self.sequence_contexts_complexities is None
            or svp_index not in self.sequence_contexts_complexities
        ):
            return None
        try:
            return pickle.loads(self.sequence_contexts_complexities[svp_index])
        except Exception:
            return None

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

    def get_sequence_from_consensus(self, consensus: Consensus) -> str:
        """Extract the sequence between the two breakends from the consensus."""
        if consensus.consensus_padding is None:
            raise ValueError("Consensus sequence is not set in the Consensus object")
        core_sequence_start: int = consensus.consensus_padding.padding_size_left
        s = self.SVprimitives[0].read_start - core_sequence_start
        e = self.SVprimitives[1].read_end - core_sequence_start
        seq = consensus.consensus_sequence[s:e]
        if self.SVprimitives[0].aln_is_reverse:
            seq = str(Seq(seq).reverse_complement())
        return seq

    def _log_id(self) -> str:
        regions = self.get_reference_regions()
        region_strs = [f"{chr}:{start}-{end}" for chr, start, end in regions]
        return f"sample={self.samplename}|consensusID={self.consensusID}|type={self.get_sv_type()}|size={self.get_size()}|regions={','.join(region_strs)}"


# Use Union for better type hints
SVpatternType = (
    SVpatternInsertion
    | SVpatternDeletion
    | SVpatternSingleBreakend
    | SVpatternInversion  # Balanced Inversion
    | SVpatternInversionDeletion  # Inverted Deletion
    | SVpatternInversionDuplication  # Inverted Duplication
    | SVpatternAdjacency  # Adjacent SVs (2-relations)
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
    # assert all(
    #     group[i].read_start <= group[i + 1].read_start for i in range(len(group) - 1)
    # ), "SVprimitives are not sorted by read_start"
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
    group: list[SVprimitive], distance_tolerance: int = 5, max_gap_size: int = 500_000
) -> dict[tuple[int, int, int, int], set[FOURRELATIONS]]:
    if len(group) < 4:
        return {}
    assert all(x.consensusID == group[0].consensusID for x in group), (
        "Not all SVprimitives have the same consensusID"
    )
    assert all(x.sv_type == 3 or x.sv_type == 4 for x in group), (
        "Not all SVprimitives are of type BND"
    )
    # assert all(
    #     group[i].read_start <= group[i + 1].read_start for i in range(len(group) - 1)
    # ), "SVprimitives are not sorted by read_start"

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
            interval_ad = Interval(
                min(a.ref_start, d.ref_start), max(a.ref_end, d.ref_end)
            )
            if len({a.chr, b.chr, c.chr, d.chr}) > 1:
                tags.add(FOURRELATIONS.HOP)
                log.debug(
                    f"HOP detected on different chromosomes: a.chr={a.chr}, b.chr={b.chr}, c.chr={c.chr}, d.chr={d.chr}"
                )
            elif not interval_ad.overlaps(
                min(b.ref_start, c.ref_end), max(b.ref_end, c.ref_start)
            ):
                tags.add(FOURRELATIONS.HOP)
                log.debug(
                    f"HOP detected on intervals: interval_ad={interval_ad}, min(b.ref_start, c.ref_end)={min(b.ref_start, c.ref_end)}, max(b.ref_end, c.ref_start)={max(b.ref_end, c.ref_start)}"
                )
            # HOP is also appiled if a,b,c,d are on the same chr, but max(a,b,c,d) - min(a,b,c,d) > max_gap_size (max SV size check)
            elif (
                max(a.ref_start, b.ref_start, c.ref_start, d.ref_start)
                - min(a.ref_start, b.ref_start, c.ref_start, d.ref_start)
                > max_gap_size
            ):
                tags.add(FOURRELATIONS.HOP)
                log.debug(
                    f"HOP detected on same chromosome but large distance (>{max_gap_size}bp): a.ref_start={a.ref_start}, b.ref_start={b.ref_start}, c.ref_start={c.ref_start}, d.ref_start={d.ref_start}"
                )
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
    # an inversion is present if there exists a 4-inv
    results: list[tuple[int, int, int, int]] = []
    for (a, b, c, d), x in fourrelations.items():
        if FOURRELATIONS.INVERSION in x and FOURRELATIONS.HOP not in x:
            # create a SVcomplex from the group
            results.append((a, b, c, d))
    return results


def possible_inversions_deletions_from_BNDs(
    svps: list[SVprimitive],
    fourrelations: dict[tuple[int, int, int, int], set[FOURRELATIONS]],
    margin: int = 5,
) -> list[tuple[int, int, int, int]]:
    """Calls inverted deletions from a group of SVprimitives. The returned list of 4-tuples contains the indices of the SVprimitives that can form an inverted deletion.

    An inverted deletion is an inversion where the inner interval (b,c) is shorter than the outer interval (a,d) by more than the margin.
    """
    if len(svps) < 4:
        return []
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError(
            "All SVprimitives must be of type BND to call inverted deletions"
        )

    results: list[tuple[int, int, int, int]] = []
    for (a, b, c, d), relations in fourrelations.items():
        if FOURRELATIONS.INVERSION not in relations or FOURRELATIONS.HOP in relations:
            continue

        svp_a, svp_b, svp_c, svp_d = svps[a], svps[b], svps[c], svps[d]

        # Calculate intervals on the reference
        interval_ad_size = abs(
            max(svp_a.ref_end, svp_d.ref_end) - min(svp_a.ref_start, svp_d.ref_start)
        )
        interval_bc_size = abs(
            max(svp_b.ref_end, svp_c.ref_end) - min(svp_b.ref_start, svp_c.ref_start)
        )

        # Check if interval (b,c) is shorter than (a,d) by more than margin
        if interval_bc_size < interval_ad_size - margin:
            results.append((a, b, c, d))

    return results


def possible_inversions_duplications_from_BNDs(
    svps: list[SVprimitive],
    fourrelations: dict[tuple[int, int, int, int], set[FOURRELATIONS]],
    margin: int = 5,
) -> list[tuple[int, int, int, int]]:
    """Calls inverted duplications from a group of SVprimitives. The returned list of 4-tuples contains the indices of the SVprimitives that can form an inverted duplication.

    An inverted duplication is an inversion where the inner interval (b,c) is longer than the outer interval (a,d) by more than the margin AND there is no HOP.
    """
    if len(svps) < 4:
        return []
    if not all(sv.sv_type == 3 or sv.sv_type == 4 for sv in svps):
        raise ValueError(
            "All SVprimitives must be of type BND to call inverted duplications"
        )

    results: list[tuple[int, int, int, int]] = []
    for (a, b, c, d), relations in fourrelations.items():
        if FOURRELATIONS.INVERSION not in relations or FOURRELATIONS.HOP in relations:
            continue

        svp_a, svp_b, svp_c, svp_d = svps[a], svps[b], svps[c], svps[d]

        # Calculate intervals on the reference
        interval_ad_size = abs(
            max(svp_a.ref_end, svp_d.ref_end) - min(svp_a.ref_start, svp_d.ref_start)
        )
        interval_bc_size = abs(
            max(svp_b.ref_end, svp_c.ref_end) - min(svp_b.ref_start, svp_c.ref_start)
        )

        # Check if interval (b,c) is longer than (a,d) by more than margin
        if interval_bc_size > interval_ad_size + margin:
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
        if x == {TWORELATIONS.REFCLOSED, TWORELATIONS.READGAP}:
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
        if x == {TWORELATIONS.READCLOSED, TWORELATIONS.REFGAP}:
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
    SVprimitives: list[SVprimitive],
    max_del_size: int = 500_000,
    log_level_override: int | None = None,
    max_fourrelations_gap_size: int = 500_000,
) -> list[SVpatternType]:
    """
    Parses Sv patterns from the SVprimitives of one consensus. All SVprimitives must have the same consensusID.
    IMPORTANT: SVprimitives must not be sorted! They are supposed to be in read order as they were parsed in the order of alignments.
    """
    if not all(svp.consensusID == SVprimitives[0].consensusID for svp in SVprimitives):
        raise ValueError(
            "All SVprimitives must have the same consensusID to parse SVpatterns"
        )
    if log_level_override is not None:
        log.setLevel(log_level_override)
    log.debug("consensusIDs = %s", {svp.consensusID for svp in SVprimitives})

    result: list[SVpatternType] = []

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

    # 4-relations based parsing for inversions and complex SVs
    if len(breakends) == 4:
        fourrelations = four_relations_of_group(
            group=breakends, max_gap_size=max_fourrelations_gap_size
        )
        log.debug(f"Four-relations: {fourrelations}")

        checks = [
            (possible_inversions_deletions_from_BNDs, SVpatternInversionDeletion),
            (possible_inversions_duplications_from_BNDs, SVpatternInversionDuplication),
            (possible_inversions_from_BNDs, SVpatternInversion),
        ]

        for check, pattern_class in checks:
            # early stop if not more than 4 unsused breakends are left
            if len(used_indices) + 4 > len(breakends):
                break
            indices = check(svps=breakends, fourrelations=fourrelations)
            for a, b, c, d in indices:
                if not (
                    a in used_indices
                    or b in used_indices
                    or c in used_indices
                    or d in used_indices
                ):
                    bnds = [breakends[i] for i in (a, b, c, d)]
                    result.append(pattern_class(SVprimitives=bnds))
                    used_indices.update({a, b, c, d})

    # 2-realtions based parsing for insertions, deletions, translocations
    if len(breakends) >= 2 and len(used_indices) < len(breakends):
        tworelations = two_relations_of_group(
            group=breakends, max_del_size=max_del_size
        )
        log.debug(f"Two-relations: {tworelations}")

        idx_tuples_insertions: list[tuple[int, int]] = possible_insertions_from_BNDs(
            svps=breakends, tworelations=tworelations
        )
        for a, b in idx_tuples_insertions:
            if a in used_indices or b in used_indices:
                continue
            used_indices.update([a, b])
            result.append(SVpatternInsertion(SVprimitives=[breakends[a], breakends[b]]))
            log.debug(f"Parsed insertion from BNDs: {result[-1]}")

        idx_tuples_deletions: list[tuple[int, int]] = possible_deletions_from_BNDs(
            svps=breakends, tworelations=tworelations
        )
        for a, b in idx_tuples_deletions:
            if a in used_indices or b in used_indices:
                continue
            used_indices.update([a, b])
            result.append(SVpatternDeletion(SVprimitives=[breakends[a], breakends[b]]))
            log.debug(f"Parsed deletion from BNDs: {result[-1]}")

    # if there is only one breakend left, it is a single-ended breakend
    if len(used_indices) - len(breakends) == 1:
        unused_index = (set(range(len(breakends))) - used_indices).pop()
        result.append(SVpatternSingleBreakend(SVprimitives=[breakends[unused_index]]))
        used_indices.add(unused_index)
        log.debug(f"Parsed single-ended breakend from BND: {result[-1]}")

    # parse all leftover breakends to single ended BNDs or Adjacency BNDs
    unused_indices = set(range(len(breakends))) - used_indices
    if unused_indices:
        # _loss_logger = get_signal_loss_logger()
        _consensusID = SVprimitives[0].consensusID if SVprimitives else ""
        # log every consensusID and every SVprimitives that is part of the unused indices for breakends
        # log.debug(
        #     f"TRANFORMED::parse_SVprimitives_to_SVpatterns:(unused break ends to adjacencies/singletons), unused_indices={sorted(unused_indices)}, svprimitives:{[breakends[i]._get_description() for i in sorted(unused_indices)]}"
        # )
        # _loss_logger.log_skipped(
        #     stage="parse_SVprimitives_to_SVpatterns",
        #     consensusID=_consensusID,
        #     reason="breakends_not_parsed_into_complex_SVpattern",
        #     count=len(unused_indices),
        #     details={
        #         "unused_breakend_indices": sorted(unused_indices),
        #         "total_breakends": len(breakends),
        #         "breakend_details": [
        #             {
        #                 "sv_type": breakends[i].sv_type,
        #                 "chr": breakends[i].chr,
        #                 "ref_start": breakends[i].ref_start,
        #                 "ref_end": breakends[i].ref_end,
        #                 "alignmentID": breakends[i].alignmentID,
        #             }
        #             for i in sorted(unused_indices)
        #         ],
        #     },
        # )
        log.info(
            f"parse_SVprimitives_to_SVpatterns: consensus {_consensusID}: "
            f"{len(unused_indices)} out of {len(breakends)} breakends could not be parsed "
            f"into complex SVpatterns and were broken up into adjacencies/singletons"
        )
        result.extend(
            break_up_CPX(SVprimitives=breakends, unused_indices=unused_indices)
        )

    return result


def break_up_CPX(
    SVprimitives: list[SVprimitive], unused_indices: set[int]
) -> list[SVpatternType]:
    """Break up a complex SVpattern into multiple SVpatternAdjacency patterns that have each two
    connected break end sv primitives or into single ended breaks.
    The set of unused indices tells if two primitives are adjacent."""
    # the already SVprimitives are sorted by location on the read. Their order must not be changed.
    # if two sv primitives are of type 3 or 4 and have different alignmentIDs and are adjacent (indices differ by 1),
    # then they are adjacencies and are grouped into one complex pattern
    # if they share the same alignmentID, the are on the same alignment fragment and not a novel adjacency.
    # if they cannot be paried with other breakends, they are single-ended breakends
    # the first and the last breakend can be single break ends
    if not all(svp.sv_type == 3 or svp.sv_type == 4 for svp in SVprimitives):
        raise ValueError(
            "All SVprimitives must be of type BND to break up complex SVpatterns"
        )

    if not unused_indices:
        return []

    result: list[SVpatternType] = []
    sorted_indices = sorted(unused_indices)
    paired_indices: set[int] = set()

    # Iterate through consecutive pairs of unused indices
    for i in range(len(sorted_indices) - 1):
        idx_a = sorted_indices[i]
        idx_b = sorted_indices[i + 1]

        # Skip if already paired
        if idx_a in paired_indices or idx_b in paired_indices:
            continue

        svp_a = SVprimitives[idx_a]
        svp_b = SVprimitives[idx_b]

        # Check if they form a novel adjacency:
        # 1. Both must be breakends (type 3 or 4)
        # 2. Must be consecutive indices (differ by 1)
        # 3. Must have different alignmentIDs
        if (abs(idx_b - idx_a) == 1) and svp_a.alignmentID != svp_b.alignmentID:
            # Create a 2-breakend adjacency pattern
            result.append(SVpatternAdjacency(SVprimitives=[svp_a, svp_b]))
            paired_indices.update([idx_a, idx_b])
            log.debug(
                f"TRANSFORMED::break_up_CPX:(broken into 2 adjacencies), first={svp_a._get_description()}, second={svp_b._get_description()}"
            )

    # Handle unpaired breakends - create single breakends
    unpaired_indices = set(sorted_indices) - paired_indices
    for idx in sorted(unpaired_indices):
        result.append(SVpatternSingleBreakend(SVprimitives=[SVprimitives[idx]]))
        log.debug(
            f"TRANSFORMED::break_up_CPX:(broken into 1 single BND), first={SVprimitives[idx]._get_description()}"
        )

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


def _unstructure_svpattern_single_breakend(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternSingleBreakend", "data": base_converter.unstructure(obj)}


def _unstructure_svpattern_inversion_deletion(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {
        "type": "SVpatternInversionDeletion",
        "data": base_converter.unstructure(obj),
    }


def _unstructure_svpattern_inversion_duplication(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {
        "type": "SVpatternInversionDuplication",
        "data": base_converter.unstructure(obj),
    }


def _unstructure_svpattern_adjacency(obj):
    base_converter = cattrs.Converter()
    base_converter.register_unstructure_hook(bytes, unstructure_bytes_field)
    return {"type": "SVpatternAdjacency", "data": base_converter.unstructure(obj)}


converter.register_unstructure_hook(
    SVpatternInsertion, _unstructure_svpattern_insertion
)
converter.register_unstructure_hook(SVpatternDeletion, _unstructure_svpattern_deletion)
converter.register_unstructure_hook(
    SVpatternInversion, _unstructure_svpattern_inversion
)
converter.register_unstructure_hook(
    SVpatternAdjacency, _unstructure_svpattern_adjacency
)
converter.register_unstructure_hook(
    SVpatternSingleBreakend, _unstructure_svpattern_single_breakend
)
converter.register_unstructure_hook(
    SVpatternInversionDeletion, _unstructure_svpattern_inversion_deletion
)
converter.register_unstructure_hook(
    SVpatternInversionDuplication, _unstructure_svpattern_inversion_duplication
)


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
    elif pattern_type == "SVpatternSingleBreakend":
        return base_converter.structure(pattern_data, SVpatternSingleBreakend)
    elif pattern_type == "SVpatternAdjacency":
        return base_converter.structure(pattern_data, SVpatternAdjacency)
    elif pattern_type == "SVpatternInversionDeletion":
        return base_converter.structure(pattern_data, SVpatternInversionDeletion)
    elif pattern_type == "SVpatternInversionDuplication":
        return base_converter.structure(pattern_data, SVpatternInversionDuplication)
    else:
        raise ValueError(f"Unknown SVpattern type: {pattern_type}")


converter.register_structure_hook(SVpatternType, structure_svpattern)


# Helper functions for SVpattern class name mapping (used by SVcomposite)
def svpattern_class_to_name(cls_obj):
    """Convert a SVpattern class to its name string."""
    if cls_obj is None:
        return None
    return cls_obj.__name__


def svpattern_name_to_class(class_name):
    """Convert a class name string back to the SVpattern class."""
    if class_name is None:
        return None

    # Map class names to actual class objects
    class_map = {
        "SVpatternInsertion": SVpatternInsertion,
        "SVpatternDeletion": SVpatternDeletion,
        "SVpatternInversion": SVpatternInversion,
        "SVpatternSingleBreakend": SVpatternSingleBreakend,
        "SVpatternAdjacency": SVpatternAdjacency,
        "SVpatternInversionDeletion": SVpatternInversionDeletion,
        "SVpatternInversionDuplication": SVpatternInversionDuplication,
    }

    if class_name not in class_map:
        raise ValueError(f"Unknown SVpattern class name: {class_name}")

    return class_map[class_name]


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
    svPatterns_file: Path | str,
    timeout: float = 10.0,
):
    """Write SVpatterns to database from file in batches to avoid memory issues.

    Args:
        database: Path to the SQLite database
        svPatterns_file: Path to a file with serialized SVpatterns (one JSON per line)
        timeout: SQLite connection timeout
        batch_size: Number of records to write per batch (deprecated, size-based batching is used instead)
    """

    svPatterns_file = Path(svPatterns_file)
    if not svPatterns_file.exists():
        raise FileNotFoundError(f"SVpatterns file not found: {svPatterns_file}")

    # SQLite max insert size is 1 GB, use 500 MB threshold for safety
    MAX_CACHE_SIZE_BYTES = 500_000_000  # 500 MB
    MAX_SINGLE_PATTERN_SIZE_BYTES = MAX_CACHE_SIZE_BYTES  # 500 MB

    query = "INSERT OR IGNORE INTO svPatterns (svPatternID, consensusID, crID, svPattern) VALUES (?, ?, ?, ?)"

    with sqlite3.connect(
        "file:" + str(database) + "?mode=rwc", uri=True, timeout=timeout
    ) as conn:
        c = conn.cursor()
        cache = []
        pattern_counter = 0
        total_cache_size = 0  # Track cache size in bytes

        with open(svPatterns_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip("\n")
                if not line:  # Skip empty lines
                    continue

                try:
                    # Deserialize SVpattern from JSON
                    svPattern_dict = json.loads(line)
                    # Use the converter to structure the SVpattern from the dict
                    svPattern = converter.structure(svPattern_dict, SVpatternType)  # type: ignore

                    # Generate a unique ID for the SVpattern
                    if svPattern.SVprimitives:
                        consensusID = svPattern.SVprimitives[0].consensusID
                        crID = int(consensusID.split(".")[0])
                        svPatternID = (
                            f"{pattern_counter}-{svPattern.get_sv_type()}-{consensusID}"
                        )
                    else:
                        log.warning(
                            f"SVpattern at line {line_num} has no primitives, skipping"
                        )
                        continue

                    # Serialize and check individual pattern size
                    pickled_pattern = pickle.dumps(converter.unstructure(svPattern))
                    svPattern_size = len(pickled_pattern)

                    if svPattern_size > MAX_SINGLE_PATTERN_SIZE_BYTES:
                        log.warning(
                            f"SVpattern {svPatternID} is very large ({float(svPattern_size) / 1_000_000.0:.2f} MB) and therefore needs to be skipped."
                        )
                        continue  # Skip overly large patterns

                    cache.append((
                        svPatternID,
                        consensusID,
                        crID,
                        pickled_pattern,
                    ))

                    # Update running cache size
                    total_cache_size += svPattern_size
                    pattern_counter += 1

                    # Write batch when cache exceeds size threshold
                    if total_cache_size >= MAX_CACHE_SIZE_BYTES:
                        log.info(
                            f"Writing batch of {len(cache)} SVpatterns to database (total size: {float(total_cache_size) / 1_000_000.0:.2f} MB)"
                        )

                        try:
                            c.executemany(query, cache)
                            conn.commit()
                            log.info(
                                f"Wrote batch of {len(cache)} SVpatterns to database with size {total_cache_size / 1_000_000.0:.2f} MB"
                            )
                        except sqlite3.Error as e:
                            log.error(
                                f"Failed to write batch ending at line {line_num}: {e}"
                            )
                            log.error(
                                f"Total batch size: {total_cache_size / 1_000_000.0:.2f} MB"
                            )
                            raise
                        finally:
                            log.info("Clearing cache after batch write")
                            cache.clear()
                            total_cache_size = 0  # Reset cache size counter

                except json.JSONDecodeError as e:
                    raise ValueError(f"Failed to parse JSON at line {line_num}: {e}")
                except Exception as e:
                    raise ValueError(
                        f"Failed to process SVpattern with ID {svPatternID} at line {line_num}: {e}"
                    )

        # Write remaining records
        if len(cache) > 0:
            log.info(
                f"Writing final batch of {len(cache)} SVpatterns (total size: {float(total_cache_size) / 1_000_000.0:.2f} MB)"
            )
            c.executemany(query, cache)
            conn.commit()
            log.info(f"Wrote final batch of {len(cache)} SVpatterns to database")
            cache.clear()
            total_cache_size = 0

        c.close()


def read_svPatterns_from_db(
    database: Path, consensusIDs: set[str] | None = None, crIDs: set[int] | None = None
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
