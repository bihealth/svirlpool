# traceback of trandem repeats from reference to read interval
# this script receives a .crs (candidate regions) file and a .bed file (tandem repeats)
# and an alignments file (in .bam format)
# and writes intervals of tandem repeats on reads to a sqlite3 database
# %%
import argparse
import sqlite3
import subprocess
import sys
import tempfile
import typing
from multiprocessing import Pool
from pathlib import Path
from shlex import split

import pandas as pd
import pysam
from logzero import logger as log
from scripts.consensus_to_IDBCs import Direction, get_ref_position_on_read
from tqdm import tqdm

from . import crs_to_dbs, datatypes, util

# %%


def traceback_tandem_repeats(
    bam_file: Path,
    region: typing.Tuple[str, int, int],
    rows: pd.DataFrame,
    dict_representatives: dict,
) -> typing.List[typing.Tuple[int, int, str, int, int, int, int, int]]:
    log.info(f"tracing back repeats on reads in region {region}..")
    # dict_tr_readnames is of the form {trID:{readnames}}
    #
    # fetch alignments in region. fill a buffer with the alignments.
    # loop over rows. each row is a tandem repeat. rows has the columns: chr,start,end,trID,crID.
    # The alignments buffer is a queue of alignments, that are sorted by their reference end position.
    read_repeat_intervals = []
    alignments_buffer = []
    # create an iterator on the fetched alignment in region
    alignments = list(pysam.AlignmentFile(bam_file, "rb").fetch(*region))
    if len(alignments) > 10_000:
        log.warning(
            f"{len(alignments)} in region {region} consume about {sys.sizeof(alignments)/1_000_000} MB"
        )
    iter_alignments = iter(alignments)
    try:
        alignment = next(iter_alignments)
    except StopIteration:
        return []
    for k, (chr, start, end, trID, crID) in rows.iterrows():
        # pop alignments from the buffer, that end before or on the start of the current repeat
        while (
            len(alignments_buffer) > 0 and alignments_buffer[0].reference_end <= start
        ):
            alignments_buffer.pop(0)
        # fill the buffer with new alignments
        alignments_to_add = []
        while True:
            # if alignment.query_name in reads_to_filter:
            if alignment.reference_start < end and alignment.reference_end > start:
                alignments_to_add.append(alignment)
            try:
                alignment = next(iter_alignments)
                continue
            except StopIteration:
                break
        # merge the two lists and sort the resulting list by reference end. This is faster then bisect.insort because insertion costs O(n).
        alignments_buffer = sorted(
            alignments_buffer + alignments_to_add, key=lambda x: x.reference_end
        )
        # now the buffer contains all alignments that overlap the current repeat
        # for each alignment, compute the interval on the read that overlaps the repeat
        for alignment in alignments_buffer:
            pos_left = get_ref_position_on_read(
                alignment=alignment, position=start, direction=Direction.LEFT
            )
            pos_right = get_ref_position_on_read(
                alignment=alignment, position=end, direction=Direction.RIGHT
            )
            read_repeat_intervals.append(
                (
                    chr,
                    start,
                    end,
                    dict_representatives[crID],
                    crID,
                    trID,
                    alignment.query_name,
                    pos_left,
                    pos_right,
                )
            )
            # chr,start,end,mcrID,crID,trID,readname,pos_left,pos_right
    return read_repeat_intervals


# implement a much faster alternative: use a set of readnames as a frist filter, then use a dictionary of the form: readname:[(chr,start,end,trID,mcrID,crID),..]
def intervals_of_tandem_repeats_on_reads_to_db(
    bam_file: Path,
    path_crs: Path,
    path_connections: Path,
    tandem_repeats_bed: Path,
    threshold_connections: int,
    threads: int,
    path_db: Path,
) -> None:
    # find all relevant tandem repeats and the reads that overlap them
    # the aim is to find the set of read-repeat pairs within a certain region.
    # this region is constructed from agglomerating adjacent repeats
    # (repeats within 50kb of each other are considered adjacent).
    # load crs
    crs: typing.List[datatypes.CandidateRegion] = util.load_crs(path_crs=path_crs)
    # get representatives
    _, UF, __ = crs_to_dbs.read_connections(
        path_connections=path_connections, threshold=threshold_connections
    )
    dict_representatives = crs_to_dbs.representative_crs(UF, [cr.crID for cr in crs])
    # write crs to bed (to intersect with repeats)
    tmp_crs_bed = tempfile.NamedTemporaryFile(suffix=".bed")
    log.info(f"writing crs to {tmp_crs_bed.name}")
    with open(tmp_crs_bed.name, "w") as f:
        for cr in crs:
            print(cr.chr, cr.referenceStart, cr.referenceEnd, cr.crID, sep="\t", file=f)
    # create tmp bed file of enumerated tandem repeats
    tmp_tandem_repeats = tempfile.NamedTemporaryFile(suffix=".bed")
    log.info(f"writing enumerated tandem repeats to {tmp_tandem_repeats.name}")
    with open(tandem_repeats_bed, "r") as f:
        with open(tmp_tandem_repeats.name, "w") as g:
            for i, line in enumerate(f.readlines()):
                # if line is empty, continue
                l = line.strip().split("\t")
                if len(l) != 3:
                    log.warning(f"line {line} does not have 3 columns")
                    continue
                print(*l, str(i), sep="\t", file=g)
    # intersect crs with tandem repeats
    cmd_intersect = split(
        f"bedtools intersect -wa -wb -a {tmp_tandem_repeats.name} -b {tmp_crs_bed.name}"
    )
    tmp_selected_tandem_repeats = tempfile.NamedTemporaryFile(suffix=".bed")
    log.info("intersecting tandem repeats with crs")
    # header: chr,start,end,trID,chr,cr_start,cr_end,crID
    # execute and write to tmp file, then load to pandas dataframe but only columns 0,1,2,3,7
    with open(tmp_selected_tandem_repeats.name, "w") as f:
        subprocess.check_call(cmd_intersect, stdout=f)
    df_selected_repeats = pd.read_csv(
        tmp_selected_tandem_repeats.name, sep="\t", header=None, usecols=[0, 1, 2, 3, 7]
    )
    df_selected_repeats.columns = ["chr", "start", "end", "trID", "crID"]
    # agglomerate adjacent repeats. Two repeats are adjacent if they are within 50kb of each other.
    # store start and end index of each agglomeration to a list
    log.info("condensing regions..")
    agglomerations = []
    current_start = 0
    for i in range(len(df_selected_repeats) - 1):
        # continue
        if (
            df_selected_repeats.iloc[i].chr == df_selected_repeats.iloc[i + 1].chr
            and abs(
                df_selected_repeats.iloc[i].end - df_selected_repeats.iloc[i + 1].start
            )
            < 50_000
        ):
            continue
        else:
            agglomerations.append((current_start, i))
            current_start = i + 1
    # add last agglomeration if the last item in agglomarations is not (current_start,len(df_selected_repeats)-1)
    last_item = (current_start, len(df_selected_repeats) - 1)
    if last_item != agglomerations[-1]:
        agglomerations.append(last_item)
    # for each agglomeration construct a set of rows (chr,start,end,trID,crID) and a region (chr,start,end)
    # the region is the min start and max end of all rows in the agglomeration
    process_regions = []
    for istart, iend in agglomerations:
        rows = df_selected_repeats.iloc[istart : iend + 1]
        region = (rows.chr.iloc[0], rows.start.min(), rows.end.max())
        process_regions.append((region, rows))
    log.info(f"condensed {len(process_regions)} regions")

    # each process_region will be given to a job and processed.
    # Each job fetches the read alignments in the region and then iterates over the rows (repeats) to find the reads that overlap them
    # if a read overlaps a repeat, the interval of the repeat on the read is written to a list.
    # this list is returned finally.
    # the arguments are passed to a process_func that will be executed in parallel or in a single thread if threads=1
    jobs = [
        {
            "bam_file": bam_file,
            "region": region,
            "rows": rows,
            "dict_representatives": dict_representatives,
        }
        for i, (region, rows) in enumerate(process_regions)
    ]
    log.info(f"processing {len(jobs)} regions")
    if threads == 1:
        results = [process_func(kwargs) for kwargs in tqdm(jobs)]
    else:
        with Pool(threads) as pool:
            results = pool.map(process_func, jobs)
    # create database and write all results to it
    log.info(f"writing results to {path_db}")
    create_db(path_db=path_db)
    for arr in results:
        write_reads_to_dbs(intervals=arr, path_db=path_db)
    print(results)


def process_func(kwargs):
    return traceback_tandem_repeats(**kwargs)


# create database with columns: chr,start,end,mcrID,crID,trID,readname,pos_left,pos_right
def create_db(path_db: Path):
    # first empty the database if it already exists
    if path_db.exists():
        path_db.unlink()
    conn = sqlite3.connect(path_db)
    conn.execute(
        """CREATE TABLE IF NOT EXISTS tr_intervals
                (chr VARCHAR(50),
                start INT,
                end INT,
                mcrID INT,
                crID INT,
                trID INT,
                readname VARCHAR(60),
                pos_left INT,
                pos_right INT)"""
    )
    conn.commit()
    conn.close()


def write_reads_to_dbs(
    intervals: typing.List[typing.Tuple[str, int, int, int, int, int, str, int, int]],
    path_db: Path,
):
    # check if intervals is empty. If so, return
    if not intervals or len(intervals) == 0:
        return
    conn = sqlite3.connect(path_db)
    cur = conn.cursor()
    cur.executemany("INSERT INTO tr_intervals VALUES (?,?,?,?,?,?,?,?,?)", intervals)
    conn.commit()
    cur.close()
    conn.close()


def run(args, **kwargs):
    intervals_of_tandem_repeats_on_reads_to_db(
        bam_file=args.path_bam,
        path_crs=args.path_crs,
        path_connections=args.path_connections,
        tandem_repeats_bed=args.tandem_repeats_bed,
        threshold_connections=args.threshold_connections,
        threads=args.threads,
        path_db=args.path_db,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Creates a database with all reads and their intervals in the tandem repeats that overlap their candidate regions."
    )
    parser.add_argument(
        "-i", "--path_crs", type=Path, required=True, help="Path to the .crs file."
    )
    parser.add_argument(
        "-c",
        "--path_connections",
        type=Path,
        required=True,
        help="Path to the .connections file.",
    )
    parser.add_argument(
        "-r",
        "--tandem_repeats_bed",
        type=Path,
        required=True,
        help="Path to the .bed file with tandem repeats.",
    )
    parser.add_argument(
        "-b", "--path_bam", type=Path, required=True, help="Path to the .bam file."
    )
    parser.add_argument(
        "-o", "--path_db", type=Path, required=True, help="Path to the .db file."
    )
    parser.add_argument(
        "-w",
        "--threshold_connections",
        type=int,
        required=False,
        default=2,
        help="Threshold for connections.",
    )
    parser.add_argument(
        "-t", "--threads", type=int, required=True, help="Number of threads to use."
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
