# combine databases
# this script combines the tables from all given databases into a single database

import argparse
import sqlite3
from pathlib import Path

from logzero import logfile
from logzero import logger as log

# %%

# path_crs = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.crs.db")
# path_consensus = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.consensuses.db")
# path_alignments = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.consensus_alignments.db")

# path_output_db = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.aligned.db")


def combine_databases(
    db_paths: list[Path], output_db_path: Path, overwrite: bool = False
):
    """
    Combine all tables from multiple SQLite databases into one output database.

    Args:
        db_paths (list): List of paths to input databases.
        output_db_path (str): Path to the output database.
    """
    if overwrite:
        if output_db_path.exists():
            output_db_path.unlink()
    # Connect to the output database
    conn_output = sqlite3.connect(output_db_path)
    cur_output = conn_output.cursor()

    # Iterate over each database path
    for i, db_path in enumerate(db_paths):
        log.info(f"Processing {db_path}...")

        # Attach the current database to the output database connection
        cur_output.execute(f"ATTACH DATABASE '{db_path}' AS db_{i}")

        # Get the list of tables in the current database
        cur_output.execute(f"SELECT name FROM db_{i}.sqlite_master WHERE type='table'")
        tables = cur_output.fetchall()

        # Iterate over each table and copy it to the output database
        for table in tables:
            table_name = table[0]
            log.info(f"Copying table {table_name} from {db_path}...")

            # Check if the table already exists in the output database (it shouldn't, as table names are unique)
            cur_output.execute(
                f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}'"
            )
            if not cur_output.fetchone():
                # Copy the table structure and data
                cur_output.execute(
                    f"""
                CREATE TABLE {table_name} AS SELECT * FROM db_{i}.{table_name}
                """
                )
            else:
                log.info(
                    f"Table {table_name} already exists in the output database. Skipping."
                )

        # Detach the current database
        cur_output.execute(f"DETACH DATABASE db_{i}")

    # Commit changes and close the connection to the output database
    conn_output.commit()
    cur_output.close()
    conn_output.close()

    log.info(f"All tables have been combined into {output_db_path} successfully.")


# # Example usage:
# db_paths = [path_consensus, path_crs, path_alignments]
# path_output_db
# combine_databases(db_paths, path_output_db)


def run(args, **kwargs):
    combine_databases(
        db_paths=args.input, output_db_path=args.output, overwrite=args.overwrite
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Combines all tables from all databases and writes to the output DB."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        nargs="+",
        help="Path to the input database with the consensus objects.",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to the output database."
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        help="Overwrite the output database if it already exists.",
    )
    parser.add_argument("-l", "--logfile", required=False, help="Path to the logfile.")
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    if args.logfile:
        logfile(str(args.logfile))
    run(args)
    return


if __name__ == "__main__":
    main()
