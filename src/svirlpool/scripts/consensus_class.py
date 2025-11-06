"""
Consensus classes and data structures for structural variant detection.

This module contains the class definitions for Consensus and ConsensusPadding,
which were moved here from the main datatypes module for better organization.
"""

import logging
import pickle
import sqlite3
from pathlib import Path
from typing import Union

import attrs
import cattrs

log = logging.getLogger(__name__)

# Import the base datatypes that these classes depend on
from . import datatypes


# Structure hook for handling str | int | float union types
def structure_string_int_float_union(obj, _):
    """Structure hook for str | int | float union types."""
    if isinstance(obj, (str, int, float)):
        return obj
    # If it's a string representation of a number, try to convert it
    if isinstance(obj, str):
        try:
            # Try to convert to int first
            if "." not in obj:
                return int(obj)
            else:
                return float(obj)
        except ValueError:
            return obj
    return obj


# Register the structure hook for the union type
cattrs.register_structure_hook(Union[str, int, float], structure_string_int_float_union)


@attrs.define
class ConsensusAlignment:
    alignment: datatypes.Alignment
    uid: int
    trf_intervals: list[tuple[int, int, int]] = []
    reference_name: str = ""
    proto_svs: list[datatypes.MergedSVSignal] = []
    core_reference_interval: tuple[int, int] = (
        0,
        0,
    )  # start,end on the reference of the core interval

    def unstructure(self):
        return cattrs.unstructure(self)


@attrs.define
class ConsensusPadding:
    sequence: str
    readname_left: str
    readname_right: str
    padding_size_left: int  # core interval start
    padding_size_right: int
    consensus_interval_on_sequence_with_padding: tuple[int, int]


@attrs.define
class ConsensusDistortion:
    position: int
    size: int
    type: int  # 0 for insertion, 1 for deletion
    readname: str
    forward: bool  # True if the read is in the forward direction


@attrs.define
class Consensus:
    ID: str
    crIDs: list[int]
    original_regions: list[tuple[str, int, int]]  # chr,start,end
    consensus_sequence: str
    consensus_padding: ConsensusPadding | None = None
    intervals_cutread_alignments: list[
        tuple[int, int, str, bool]
    ] = []  # start,end,readname,forward
    cut_read_alignment_signals: list[
        datatypes.ReadAlignmentSignals
    ] = []  # ReadAlignmentSignals objects of the aligned cut reads to their consensus
    clustering_meta_data: dict[
        str, str | int | float
    ] = {}  # metadata from the clustering step

    def unstructure(self):
        return cattrs.unstructure(self)

    def get_used_readnames(self) -> set[str]:
        """Returns a set of read names that were used to generate this consensus."""
        return {
            readname
            for start, end, readname, forward in self.intervals_cutread_alignments
        }

    def get_consensus_distortions(self) -> list[ConsensusDistortion]:
        if len(self.cut_read_alignment_signals) == 0:
            return []
        return [
            ConsensusDistortion(
                position=signal.ref_start,
                size=signal.size,
                type=signal.sv_type,
                readname=ras.read_name,
                forward=ras.alignment_forward,
            )
            for ras in self.cut_read_alignment_signals
            for signal in ras.SV_signals
            if signal.sv_type in (0, 1)
        ]

    def get_max_cutread_coverage(self) -> int:
        # compute the maximum depth of overlaps of cut reads on the consensus
        events = []
        for start, end, readname, forward in self.intervals_cutread_alignments:
            events.append((start, 1))
            events.append((end, -1))
        events = sorted(events)
        max_depth = 0
        current_depth = 0
        for position, event_type in events:
            current_depth += event_type
            max_depth = max(max_depth, current_depth)
        return max_depth


@attrs.define
class CrsContainerResult:
    consensus_dicts: dict[str, Consensus]  # consenus ID to consensus object
    unused_reads: dict[int, list[datatypes.SequenceObject]]  # read ID to read object

    def unstructure(self):
        return cattrs.unstructure(self)


def create_consensus_db(database: Path):
    """Create database for storing Consensus objects."""
    with sqlite3.connect("file:" + str(database) + "?mode=rwc", uri=True) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS consensus (
                consensusID VARCHAR(50) PRIMARY KEY,
                crIDs TEXT NOT NULL,
                consensus_data BLOB NOT NULL
            )
        """
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_consensus_id ON consensus (consensusID)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_consensus_crids ON consensus (crIDs)"
        )
        conn.commit()
    log.debug(f"Created consensus database at {database}")


def write_consensus_to_db(
    database: Path, data: list[Consensus], timeout: float = 10.0, batch_size: int = 1000
):
    """Write Consensus objects to database in batches.

    Args:
        database: Path to the SQLite database
        data: List of Consensus objects to write
        timeout: SQLite connection timeout
        batch_size: Number of records to write per batch (default: 1000)
    """
    query = "INSERT OR REPLACE INTO consensus (consensusID, crIDs, consensus_data) VALUES (?, ?, ?)"

    with sqlite3.connect(
        "file:" + str(database) + "?mode=rwc", uri=True, timeout=timeout
    ) as conn:
        cursor = conn.cursor()

        for i in range(0, len(data), batch_size):
            batch = data[i : i + batch_size]
            batch_data = []

            for consensus in batch:
                # Serialize crIDs as comma-separated string for indexing
                crids_str = ",".join(str(crid) for crid in consensus.crIDs)
                # Serialize the entire consensus object using cattrs
                consensus_blob = cattrs.unstructure(consensus)
                consensus_pickled = pickle.dumps(consensus_blob)

                batch_data.append((consensus.ID, crids_str, consensus_pickled))

            cursor.executemany(query, batch_data)
            conn.commit()

            log.debug(
                f"Wrote batch {i // batch_size + 1}/{(len(data) - 1) // batch_size + 1} ({len(batch)} consensus objects)"
            )

    log.info(f"Successfully wrote {len(data)} consensus objects to database")


def read_consensus_from_db(
    database: Path,
    consensusIDs: list[str] | None = None,
    crIDs: set[int] | None = None,
    timeout: float = 30.0,
) -> list[Consensus]:
    """Read Consensus objects from database.

    Args:
        database: Path to the SQLite database
        consensusIDs: Optional list of consensus IDs to filter by
        crIDs: Optional set of crIDs to filter by (returns consensus containing any of these crIDs)
        timeout: SQLite connection timeout

    Returns:
        List of Consensus objects matching the filters, or all if no filters provided
    """
    with sqlite3.connect(
        "file:" + str(database) + "?mode=r", uri=True, timeout=timeout
    ) as conn:
        cursor = conn.cursor()

        # Build query based on filters
        if consensusIDs is not None and len(consensusIDs) > 0:
            placeholders = ",".join("?" * len(consensusIDs))
            query = f"SELECT consensus_data FROM consensus WHERE consensusID IN ({placeholders})"
            cursor.execute(query, consensusIDs)
        elif crIDs is not None and len(crIDs) > 0:
            # For crIDs, we need to check if any of the crIDs are in the comma-separated string
            # Using LIKE with OR conditions
            conditions = []
            params = []
            for crid in crIDs:
                # Match crID at start, middle, or end of comma-separated list
                conditions.append(
                    "(crIDs LIKE ? OR crIDs LIKE ? OR crIDs LIKE ? OR crIDs = ?)"
                )
                crid_str = str(crid)
                params.extend([
                    f"{crid_str},%",
                    f"%,{crid_str},%",
                    f"%,{crid_str}",
                    crid_str,
                ])

            query = (
                f"SELECT consensus_data FROM consensus WHERE {' OR '.join(conditions)}"
            )
            cursor.execute(query, params)
        else:
            # Return all consensus objects
            query = "SELECT consensus_data FROM consensus"
            cursor.execute(query)

        # Deserialize all results
        consensus_objects = []
        for row in cursor.fetchall():
            consensus_pickled = row[0]
            consensus_dict = pickle.loads(consensus_pickled)
            consensus_obj = cattrs.structure(consensus_dict, Consensus)
            consensus_objects.append(consensus_obj)

        log.debug(f"Retrieved {len(consensus_objects)} consensus objects from database")

    return consensus_objects
