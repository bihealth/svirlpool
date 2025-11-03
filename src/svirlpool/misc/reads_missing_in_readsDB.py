# %%
# this script finds all readnames that are in the intervals DB but not in the reads DB
import argparse
import sqlite3
from pathlib import Path

from logzero import logger as log
from tqdm import tqdm


# %%
def find_readnames_not_in_reads_db(path_intervals_db: Path, path_reads_db: Path):

    conn_reads = sqlite3.connect(path_reads_db)
    conn_intervals = sqlite3.connect(path_intervals_db)

    # get all readnames from the intervals DB
    log.info("Getting readnames from intervals DB")
    readnames_intervals = set()
    readnames_sampleID = {}
    cursor = conn_intervals.cursor()
    cursor.execute("SELECT readname FROM intervals")
    for row in tqdm(cursor):
        readnames_intervals.add(row[0])
        readnames_sampleID[row[0]] = row[1]
    cursor.close()

    # get all readnames from the reads DB
    log.info("Getting readnames from reads DB")
    readnames_reads = set()
    cursor = conn_reads.cursor()
    cursor.execute("SELECT readname FROM reads")
    for row in tqdm(cursor):
        readnames_reads.add(row[0])
    cursor.close()

    log.info("Comparing readnames")
    a_not_b = readnames_intervals - readnames_reads
    a_and_b = readnames_intervals & readnames_reads
    b_not_a = readnames_reads - readnames_intervals

    # print results
    print(f"Readnames in intervals DB but not in reads DB: {len(a_not_b)}")
    print(f"Readnames in both DBs: {len(a_and_b)}")
    print(f"Readnames in reads DB but not in intervals DB: {len(b_not_a)}")
    log.info("Done")


def run(args, **kwargs):
    find_readnames_not_in_reads_db(
        path_intervals_db=args.intervals_db, path_reads_db=args.reads_db
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Find number of readnames that are in the intervals DB but not in the reads DB or vice versa."
    )
    parser.add_argument(
        "-i",
        "--intervals_db",
        type=Path,
        required=True,
        help="path to the intervals DB",
    )
    parser.add_argument(
        "-r", "--reads_db", type=Path, required=True, help="path to the reads DB"
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
