# %%
import argparse
import logging as log
import pickle
import sqlite3
import tempfile
from pathlib import Path, PosixPath

from numpy import uint64
from xxhash import xxh64

from ..util import datatypes, util

log.basicConfig(level=log.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# %%


def parse_lines(byte_lines):
    parsed_data = []
    for i, line in enumerate(byte_lines):
        try:
            decoded_str = line.decode().strip()  # Decode bytes and strip newline
            parts = decoded_str.split("\t")  # Split by tab

            if len(parts) != 4:
                raise ValueError(f"Line {i} does not have exactly 4 columns: {parts}")

            # Convert values with error handling
            chr = parts[0]  # Keep as string
            start = int(parts[1])  # Convert to int
            end = int(parts[2])  # Convert to int
            readname = uint64(
                parts[3]
            )  # parts[3] is already the hash value as string, just convert to uint64
            parsed_data.append([chr, start, end, readname])
        except (UnicodeDecodeError, ValueError) as e:
            print(f"Skipping line {i} due to error: {e}")
    return parsed_data


def load_cov_from_db(db_path: Path, itemname: str = "effective") -> list:
    """
    Reads the chunked coverage data from the database and parses the lines.
    """
    conn = sqlite3.connect("file:" + str(db_path) + "?mode=ro", uri=True)
    cursor = conn.cursor()

    # Get number of chunks
    cursor.execute(
        "SELECT data FROM coverage WHERE id = ?", (f"{itemname}_num_chunks",)
    )
    row = cursor.fetchone()
    if row is None:
        conn.close()
        return []

    num_chunks = pickle.loads(row[0])

    # Load all chunks
    all_lines = []
    for chunk_idx in range(num_chunks):
        chunk_id = f"{itemname}_chunk_{chunk_idx}"
        cursor.execute("SELECT data FROM coverage WHERE id = ?", (chunk_id,))
        row = cursor.fetchone()
        if row:
            all_lines.extend(pickle.loads(row[0]))

    conn.close()
    return parse_lines(all_lines)


def create_database(output: bytes) -> None:
    # create a new database file if it does not exist
    # with table "coverage" with column "id" that is always 0 and "data" that is a BLOB
    with sqlite3.connect("file:" + str(output) + "?mode=rwc", uri=True) as conn:
        conn.execute(
            """CREATE TABLE IF NOT EXISTS coverage
            (id VARCHAR(40) PRIMARY KEY,
            data BLOB)"""
        )
        conn.commit()


# reading with util.yield_cov_from_db
def write_cov_to_db(
    db_path: Path, tmp_bg_file: str | Path, itemname: str, chunk_size: int = 1_000_000
) -> None:
    """
    Write the tmp_bg_file to the database as chunked BLOBs.

    :param db_path: Path to the SQLite database
    :param tmp_bg_file: Path to the temporary coverage file
    :param itemname: Identifier for this coverage data
    :param chunk_size: Number of lines per chunk (default: 1,000,000)
    """
    conn = sqlite3.connect("file:" + str(db_path) + "?mode=rwc", uri=True)

    with open(tmp_bg_file, "rb") as file:
        chunk_idx = 0
        while True:
            # Read chunk_size lines
            lines = []
            for _ in range(chunk_size):
                line = file.readline()
                if not line:
                    break
                lines.append(line)

            if not lines:
                break  # No more data

            # Compress and store chunk
            data = pickle.dumps(lines)
            chunk_id = f"{itemname}_chunk_{chunk_idx}"
            conn.execute(
                "INSERT OR REPLACE INTO coverage (id, data) VALUES (?, ?)",
                (chunk_id, data),
            )
            chunk_idx += 1

        # Store metadata about number of chunks
        conn.execute(
            "INSERT OR REPLACE INTO coverage (id, data) VALUES (?, ?)",
            (f"{itemname}_num_chunks", pickle.dumps(chunk_idx)),
        )

    conn.commit()
    conn.close()


# can add another item with total raf coverage, not just the effective intervals' coverage. Test first, then augment.
def rafs_to_coverage(input: Path, output: Path) -> None:
    """Iterates all rafs and writes all intervals to a tmp file. The intervals are then written to a database"""
    # tmp_intervals_total = tempfile.NamedTemporaryFile(delete=True, suffix=".bed")
    tmp_intervals_effective = tempfile.NamedTemporaryFile(delete=True, suffix=".bed.gz")
    # with gzip.open(tmp_intervals_total.name, "wt") as f:
    with open(tmp_intervals_effective.name, "w") as g:
        for raf in util.yield_from_raf(input):
            x: datatypes.ReadAlignmentFragment = raf
            # print(x.reference_name, x.reference_alignment_start, x.reference_alignment_end, x.read_name, sep="\t", file=f)
            if x.effective_interval:
                print(
                    *x.effective_interval,
                    uint64(xxh64(x.read_name).intdigest()),
                    sep="\t",
                    file=g,
                )

    # create DB
    create_database(output=output)
    # write both (total, ) to DB
    # write_cov_to_db(db_path=output, tmp_bg_file=tmp_intervals_total.name, itemname="total")
    write_cov_to_db(
        db_path=output, tmp_bg_file=tmp_intervals_effective.name, itemname="effective"
    )
    # remove tmp files
    # tmp_intervals_total.close()
    tmp_intervals_effective.close()


# %%


def run(args, **kwargs):
    rafs_to_coverage(input=args.input, output=args.output)
    return


def get_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        type=PosixPath,
        required=True,
        help="Path to read alignments file (rafs), optimally after filtering.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=PosixPath,
        required=True,
        help="Path to output database file. Is created if it does not already exist. Adds a table 'effective' by default.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
