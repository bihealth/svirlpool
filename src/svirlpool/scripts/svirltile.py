# this module provides functionality to create a database for a sample. The database should contain:
# 1) svpatterns (taken from another database)
# 2) coverage data (taken from another database)
# 3) consensus sequences (taken from an uncompressed fasta file)
# 4) metadata (currently just the sample name, but it should be a python dict)

import argparse
import logging
import pickle
import sqlite3
from pathlib import Path

from Bio import SeqIO

from . import combine_databases

log = logging.getLogger(__name__)


def create_sample_database(
    svpatterns_db_path: str | Path,
    coverage_db_path: str | Path,
    consensus_fasta_path: str | Path,
    samplename: str,
    output_db_path: str | Path,
    copynumber_db_path: str | Path | None = None,
):
    """
    Create a sample database containing svpatterns, coverage data, consensus sequences, copy number tracks, and metadata.

    Args:
        svpatterns_db_path (str): Path to the svpatterns database.
        coverage_db_path (str): Path to the coverage database.
        consensus_fasta_path (str): Path to the consensus FASTA file.
        samplename (str): Name of the sample.
        output_db_path (str): Path to save the output database.
        copynumber_db_path (str, optional): Path to the copy number tracks database.

    Returns:
        None
    """
    output_path = Path(output_db_path)

    # Step 1: Combine svpatterns, coverage, and optionally copynumber databases
    db_paths_to_combine = [Path(svpatterns_db_path), Path(coverage_db_path)]

    if copynumber_db_path is not None:
        copynumber_path = Path(copynumber_db_path)
        if copynumber_path.exists():
            db_paths_to_combine.append(copynumber_path)
            log.debug(f"Including copy number database: {copynumber_db_path}")
        else:
            log.warning(
                f"Copy number database not found: {copynumber_db_path}. Continuing without it."
            )

    log.debug(f"Combining databases: {db_paths_to_combine}")
    combine_databases.combine_databases(
        db_paths=db_paths_to_combine, output_db_path=output_path, overwrite=True
    )

    # Step 2: Add consensus sequences from FASTA file
    log.debug(f"Adding consensus sequences from {consensus_fasta_path}")
    _add_consensus_sequences_to_db(
        db_path=output_path, fasta_path=Path(consensus_fasta_path)
    )

    # Step 3: Add metadata
    log.debug(f"Adding metadata for sample: {samplename}")
    _add_metadata_to_db(db_path=output_path, samplename=samplename)

    log.debug(f"Sample database created successfully: {output_db_path}")


def _add_consensus_sequences_to_db(db_path: Path, fasta_path: Path):
    with sqlite3.connect(str(db_path)) as conn:
        # Create consensus_sequences table
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS consensus_sequences (
                consensusID VARCHAR(50) PRIMARY KEY,
                sequence BLOB NOT NULL,
                description TEXT
            )
        """
        )

        # Read FASTA file and insert sequences
        sequences_data = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Compress sequence with pickle before storing
            compressed_sequence = pickle.dumps(str(record.seq))
            sequences_data.append((
                record.id,
                compressed_sequence,
                record.description if record.description else "",
            ))

        # Insert all sequences
        conn.executemany(
            """
            INSERT OR REPLACE INTO consensus_sequences
            (consensusID, sequence, description)
            VALUES (?, ?, ?)
        """,
            sequences_data,
        )

        # Create index for faster lookups
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_consensus_sequences_id ON consensus_sequences (consensusID)"
        )
        conn.commit()

        log.debug(f"Added {len(sequences_data)} consensus sequences to database")


def _add_metadata_to_db(db_path: Path, samplename: str):
    """
    Add sample metadata to the database.

    Args:
        db_path (Path): Path to the database.
        samplename (str): Name of the sample.
    """
    metadata = {"samplename": samplename}

    with sqlite3.connect(str(db_path)) as conn:
        # Create metadata table
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS metadata (
                key VARCHAR(50) PRIMARY KEY,
                value TEXT NOT NULL
            )
        """
        )

        # Insert metadata
        for key, value in metadata.items():
            conn.execute(
                """
                INSERT OR REPLACE INTO metadata (key, value)
                VALUES (?, ?)
            """,
                (key, str(value)),
            )

        conn.commit()

        log.debug(f"Added metadata for sample: {samplename}")


def get_consensus_sequences(db_path: Path, consensus_ids: list[str] | None = None):
    """
    Generator that yields consensus sequences one at a time.

    Args:
        db_path: Path to the database
        consensus_ids: Optional list of consensus IDs to filter. If None, yields all sequences.

    Yields:
        tuple[str, dict]: Tuple of (consensusID, dict) where dict contains
                         'id', 'description', 'sequence' (compressed).
    """
    with sqlite3.connect(str(db_path)) as conn:
        c = conn.cursor()
        try:
            if consensus_ids is not None and len(consensus_ids) > 0:
                placeholders = ",".join("?" * len(consensus_ids))
                c.execute(
                    f"SELECT consensusID, sequence, description FROM consensus_sequences WHERE consensusID IN ({placeholders})",
                    consensus_ids,
                )
            else:
                c.execute(
                    "SELECT consensusID, sequence, description FROM consensus_sequences"
                )

            count = 0
            for row in c.fetchall():
                consensus_id, compressed_sequence, description = row
                yield (
                    consensus_id,
                    {
                        "id": consensus_id,
                        "sequence": compressed_sequence,  # Keep compressed for efficiency
                        "description": description,
                    },
                )
                count += 1
            log.debug(f"Yielded {count} consensus sequences from database")
        except sqlite3.Error as e:
            log.error(f"Error reading consensus sequences from database: {e}")
            raise e
        finally:
            c.close()


def get_metadata(db_path: Path) -> dict[str, str]:
    """Retrieve metadata from the database ( dict[str, str]: Dictionary with metadata key-value pairs )"""
    result = {}

    with sqlite3.connect(str(db_path)) as conn:
        c = conn.cursor()
        try:
            c.execute("SELECT key, value FROM metadata")
            for row in c.fetchall():
                key, value = row
                result[key] = value
        except sqlite3.Error as e:
            log.error(f"Error reading metadata from database: {e}")
            raise e
        finally:
            c.close()

    log.debug(f"Retrieved {len(result)} metadata entries from database")
    return result


def decompress_sequence(compressed_sequence: bytes) -> str:
    return pickle.loads(compressed_sequence)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Align consensus sequences to reference and extract SV primitives"
    )
    parser.add_argument("--samplename", type=str, required=True, help="Sample name.")
    parser.add_argument(
        "--consensus", type=Path, required=True, help="Path to consensus FASTA file."
    )
    parser.add_argument(
        "--svpatterns-db", type=Path, required=True, help="Path to svpatterns database."
    )
    parser.add_argument(
        "--coverage-db", type=Path, required=True, help="Path to coverage database."
    )
    parser.add_argument(
        "--copynumber-db",
        type=Path,
        default=None,
        help="Path to copy number tracks database (optional).",
    )
    parser.add_argument(
        "--output-db", type=Path, required=True, help="Path to output database."
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Path to log file. If not provided, logs will be printed to stdout.",
    )
    return parser


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if args.log_file:
        # Add file handler if log file is specified
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setLevel(getattr(logging, args.log_level.upper()))
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

    create_sample_database(
        svpatterns_db_path=args.svpatterns_db,
        coverage_db_path=args.coverage_db,
        consensus_fasta_path=args.consensus,
        samplename=args.samplename,
        output_db_path=args.output_db,
        copynumber_db_path=args.copynumber_db,
    )
    return


if __name__ == "__main__":
    main()
