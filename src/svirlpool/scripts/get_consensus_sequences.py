import argparse
import logging
import shutil
import sqlite3
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from ..localassembly import svirltile


def write_consensus_sequences(
    db_path: Path,
    output_fasta_path: Path,
    batch_size: int = 1000,
    consensus_ids: list[str] | None = None,
):
    """Write consensus sequences from the database to a FASTA file.

    Args:
        db_path: Path to the database
        output_fasta_path: Path to output FASTA file
        batch_size: Number of sequences to batch before writing
        consensus_ids: Optional list of consensus IDs to extract. If None, extracts all sequences.
    """
    is_gzipped = output_fasta_path.suffix == ".gz"

    # Create a temporary file for writing sequences
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_file:
        tmp_path = Path(tmp_file.name)

        try:
            logging.info("counting the sequences to be extracted..")
            # Get total count for progress bar
            with sqlite3.connect(str(db_path)) as conn:
                cursor = conn.cursor()
                if consensus_ids is not None and len(consensus_ids) > 0:
                    placeholders = ",".join("?" * len(consensus_ids))
                    cursor.execute(
                        f"SELECT COUNT(*) FROM consensus_sequences WHERE consensusID IN ({placeholders})",
                        consensus_ids,
                    )
                else:
                    cursor.execute("SELECT COUNT(*) FROM consensus_sequences")
                total = cursor.fetchone()[0]

            # Write sequences using generator with batching
            records_written = 0
            batch = []

            for _consensus_id, consensus_data in tqdm(
                svirltile.get_consensus_sequences(db_path, consensus_ids=consensus_ids),
                total=total,
                desc="Writing consensus sequences",
            ):
                # Decompress sequence
                sequence = svirltile.decompress_sequence(consensus_data["sequence"])

                # Create SeqRecord
                record = SeqRecord(
                    seq=Seq(sequence),
                    id=consensus_data["id"],
                    description=consensus_data["description"],
                )

                batch.append(record)

                # Write batch when it reaches batch_size
                if len(batch) >= batch_size:
                    SeqIO.write(batch, tmp_file, "fasta")
                    records_written += len(batch)
                    batch = []

            # Write remaining records
            if batch:
                SeqIO.write(batch, tmp_file, "fasta")
                records_written += len(batch)

            logging.info(
                f"Wrote {records_written} consensus sequences to temporary file"
            )

        except Exception as e:
            # Clean up temp file on error
            tmp_path.unlink(missing_ok=True)
            raise e

    # Process the temporary file based on output format
    if is_gzipped:
        logging.info("Compressing FASTA file with bgzip...")
        try:
            # Bgzip the temporary file to the output path
            subprocess.run(
                ["bgzip", "-c", str(tmp_path)],
                stdout=open(output_fasta_path, "wb"),
                check=True,
                stderr=subprocess.PIPE,
                text=False,
            )
            logging.info(f"Successfully compressed to {output_fasta_path}")

            # Index with samtools faidx
            logging.info("Indexing FASTA file with samtools faidx...")
            subprocess.run(
                ["samtools", "faidx", str(output_fasta_path)],
                check=True,
                capture_output=True,
                text=True,
            )
            logging.info(f"Successfully created index: {output_fasta_path}.fai")
        except subprocess.CalledProcessError as e:
            logging.error(
                f"Failed to compress or index FASTA file: {e.stderr if hasattr(e, 'stderr') else str(e)}"
            )
            raise
        except FileNotFoundError as e:
            logging.error(
                f"Required tool not found in PATH (bgzip or samtools): {str(e)}"
            )
            raise
        finally:
            # Clean up temp file
            tmp_path.unlink(missing_ok=True)
    else:
        # Just copy the temp file to the output path
        logging.info(f"Copying FASTA file to {output_fasta_path}")
        shutil.move(str(tmp_path), str(output_fasta_path))
        logging.info(f"Successfully wrote {output_fasta_path}")


def run(args):
    """Run function compatible with argparse from __main__.py"""
    logging.basicConfig(level=logging.INFO)
    consensus_ids = (
        args.consensus_ids
        if hasattr(args, "consensus_ids") and args.consensus_ids
        else None
    )
    write_consensus_sequences(
        Path(args.input),
        Path(args.output),
        batch_size=args.batch_size,
        consensus_ids=consensus_ids,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Extract consensus sequences from a SVIRLPOOL database."
    )
    parser.add_argument(
        "database", type=Path, help="Path to the SVIRLPOOL database file."
    )
    parser.add_argument(
        "output_fasta",
        type=Path,
        help="Path to the output FASTA file for consensus sequences.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Number of sequences to batch before writing to file (default: 1000).",
    )
    parser.add_argument(
        "--consensus-ids",
        type=str,
        nargs="+",
        default=None,
        help="Optional list of consensus IDs to extract. If not provided, all sequences will be extracted.",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    write_consensus_sequences(
        args.database,
        args.output_fasta,
        batch_size=args.batch_size,
        consensus_ids=args.consensus_ids,
    )


if __name__ == "__main__":
    main()
