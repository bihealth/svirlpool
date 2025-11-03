# annotate Consensus Objects
# - align all consensus sequences to the target reference
# - anootate the Consensus objects with the target sequence
# - and the alignments
# - write all annotated Consensus Objects to a database

# %%
import argparse
import json
import pickle
import sqlite3
from pathlib import Path, PosixPath

import cattrs
from logzero import logfile
from logzero import logger as log
from tqdm import tqdm

from . import consensus_class, datatypes

# %%


def parse_crs_container_results(path: PosixPath):
    # check if file is gzipped
    if path.suffix == ".gz":
        import gzip

        open_file = gzip.open
    else:
        open_file = open
    with open_file(path, "rt") as f:
        for line in f:
            crs_container_result = json.loads(line)
            yield cattrs.structure(
                crs_container_result, consensus_class.CrsContainerResult
            )

    # with open(path, "r") as f:
    #     for line in f:
    #         crs_container_result = json.loads(line)
    #         yield cattrs.structure(crs_container_result, consensus_class.CrsContainerResult)


def construct_final_database(path_database: Path) -> None:
    # use sqlite3 to create a database with tables:
    # - consensusess [ID,consensus object] (from crsContainerResult)
    # - unused_reads [readname,crIDs,SequenceObject] (from crsContainerResult) # crIDs is a list of INTs
    conn = sqlite3.connect("file:" + str(path_database) + "?mode=rwc", uri=True)
    conn.execute(
        """CREATE TABLE IF NOT EXISTS consensuses 
        (id VARCHAR(200) PRIMARY KEY,
        consensus BLOB)"""
    )
    conn.commit()
    conn.execute(
        """CREATE TABLE IF NOT EXISTS unused_reads
        (readname_crID VARCHAR(200) PRIMARY KEY,
        crID INTEGER,
        sequenceObject BLOB)"""
    )
    conn.commit()
    conn.close()


def write_consensus_crsContainerResult_to_database(
    path_database: Path, crsContainerResults: list[consensus_class.CrsContainerResult]
):
    # iterate over crsContainerResults
    # and write all consensus objects to the database
    # and write all unused reads to the database. append to the aug_name of each read its crID
    # log.info(f"writing {len(crsContainerResults)} consensus objects to the database")
    conn = sqlite3.connect("file:" + str(path_database) + "?mode=rwc", uri=True)
    for crsContainerResult in tqdm(crsContainerResults):
        unused_reads_data = []
        crIDs_unused_reads = []
        for crID, list_sequenceObjects in crsContainerResult.unused_reads.items():
            if len(list_sequenceObjects) > 0:
                crIDs_unused_reads.append(crID)
            for sequenceObject in list_sequenceObjects:
                unused_reads_data.append(
                    [
                        f"{sequenceObject.name}.{crID}",
                        int(crID),
                        pickle.dumps(sequenceObject.unstructure()),
                    ]
                )
        c = conn.cursor()
        c.executemany(
            "INSERT OR REPLACE INTO unused_reads (readname_crID,crID,sequenceObject) VALUES (?,?,?)",
            unused_reads_data,
        )
        # [[f"{sequenceObject.name}.{sequenceObject.description}",int(sequenceObject.description),pickle.dumps(sequenceObject.unstructure())]
        #     for sequenceObject in crsContainerResult.unused_reads])
        conn.commit()

        present_consensusIDs = [
            consensusID
            for consensusID, consensus in crsContainerResult.consensus_dicts.items()
            if consensus
        ]
        if len(present_consensusIDs) == 0:
            log.warning(
                f"no consensus objects present in this single crsContainerResult {crsContainerResult.consensus_dicts.keys()}. Unused reads written to database for crIDs: {crIDs_unused_reads}"
            )
            continue

        # log.info(f"writing {[consensusID for consensusID,consensus in crsContainerResult.consensus_dicts.items() if consensus]} consensus objects to the database")
        c = conn.cursor()
        # write all consensus objects to the database
        c.executemany(
            "INSERT OR REPLACE INTO consensuses (id,consensus) VALUES (?,?)",
            [
                (consensusID, pickle.dumps(cattrs.unstructure(consensus)))
                for consensusID, consensus in crsContainerResult.consensus_dicts.items()
            ],
        )
        # write all unused reads to the database
        conn.commit()
    c.close()
    conn.close()


def load_all_consensusIDs_from_database(path_database: Path) -> list[str]:
    conn = sqlite3.connect("file:" + str(path_database) + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT id FROM consensuses")
    rows = c.fetchall()
    c.close()
    conn.close()
    return [row[0] for row in rows]


# def load_consensus_from_db(consensusID:str,path_database:Path) -> consensus_class.Consensus:
#     conn = sqlite3.connect('file:'+str(path_database)+'?mode=ro',uri=True)
#     c = conn.cursor()
#     c.execute("SELECT consensus FROM consensuses WHERE id = ?",(consensusID,))
#     row = c.fetchone()
#     c.close()
#     conn.close()
#     if row is None:
#         raise ValueError(f"consensus with ID {consensusID} not found in database {path_database}")
#     return cattrs.structure(pickle.loads(row[0]),consensus_class.Consensus)

# def load_all_consensuses_from_db(path_database:Path) -> list[consensus_class.Consensus]:
#     conn = sqlite3.connect('file:'+str(path_database)+'?mode=ro',uri=True)
#     c = conn.cursor()
#     try:
#         c.execute("SELECT consensus FROM consensuses")
#         rows = c.fetchall()
#         c.close()
#         conn.close()
#     except:
#         # list al tables in the database and print them alongside the error message
#         c.execute("SELECT name FROM sqlite_master WHERE type='table';")
#         tables = c.fetchall()
#         c.close()
#         conn.close()
#         raise ValueError(f"error loading consensuses from database {path_database}. Tables in database: {tables}")
#     return [cattrs.structure(pickle.loads(row[0]),consensus_class.Consensus) for row in tqdm(rows)]

# redundant function. its already preent in util
# def yield_all_consensuses_from_db(path_database:Path,n_rows:int=100) -> Generator[consensus_class.Consensus,None,None]:
#     conn = sqlite3.connect('file:'+str(path_database)+'?mode=ro',uri=True)
#     c = conn.cursor()
#     c.execute("SELECT consensus FROM consensuses")
#     while True:
#         rows = c.fetchmany(n_rows)
#         if not rows:
#             break
#         for row in rows:
#             yield cattrs.structure(pickle.loads(row[0]),consensus_class.Consensus)
#     c.close()
#     conn.close()


def load_unused_reads_from_db(
    crID: int, path_database: Path
) -> list[datatypes.SequenceObject]:
    conn = sqlite3.connect("file:" + str(path_database) + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT sequenceObject FROM unused_reads WHERE crID = ?", (crID,))
    rows = c.fetchall()
    c.close()
    conn.close()
    return [
        cattrs.structure(pickle.loads(row[0]), datatypes.SequenceObject) for row in rows
    ]


# %%

# path_cat_CrsContainerResults = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.tiles.txt")
# path_target_reference = Path("/home/vinzenz/development/LRSV-detection/development/reference/hs37d5/hs37d5.fa")
# path_initial_reference = Path("/home/vinzenz/development/LRSV-detection/development/reference/hs37d5/hs37d5.fa")
# path_consensus_fasta = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.consensus.fasta")
# consensus_to_target_alignments = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.consensus.target.bam")
# consensus_to_initial_alignments = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.consensus.initial.bam")
# path_database=Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/parametertuning/d03/HG002.final_db.sqlite")
# threads=8
# aln_args=''
# tmp_dir_path:Path|None=None


def consensus_containers_to_db(
    path_cat_CrsContainerResults: Path, path_database: Path
) -> None:
    log.info(f"constructing final database {path_database}")
    construct_final_database(path_database=path_database)

    counter = 0
    batch_size: int = 1000
    results: list[consensus_class.CrsContainerResult] = []
    for i, crsContainerResult in tqdm(
        enumerate(parse_crs_container_results(path_cat_CrsContainerResults))
    ):
        results.append(crsContainerResult)
        if i % batch_size == batch_size:
            counter += batch_size
            write_consensus_crsContainerResult_to_database(
                path_database=path_database, crsContainerResults=results
            )
            results = []
    if len(results) > 0:
        counter += len(results)
        write_consensus_crsContainerResult_to_database(
            path_database=path_database, crsContainerResults=results
        )
        results = []
    log.info(f"written {counter} annotated consensus objects to the database")


# put everything together

# %%


def run(args, **kwargs):
    consensus_containers_to_db(
        path_cat_CrsContainerResults=args.input, path_database=args.output
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Reads crs container objects and creates annotated Consensus objects that are written to the output database."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the concatenated crsContainerResults file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output path to the final database with tables consensuses and unused_reads.",
    )
    parser.add_argument(
        "--logfile",
        type=Path,
        required=False,
        default=None,
        help="Path to the logfile.",
    )
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
