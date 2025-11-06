# %%
import argparse
import json
import sqlite3

# from tqdm import tqdm
from itertools import combinations
from pathlib import Path, PosixPath

# from sklearn.cluster import AgglomerativeClustering
# import numpy as np
# import pandas as pd
from logzero import logger as log
from tqdm import tqdm

from . import datastructures, datatypes, signalstrength_to_crs

# %%


def serialize_crs_container(crs_container: dict) -> str:
    # debugging block, find the numpy datatypes
    # for cr in crs_container['crs']:
    #     try:
    #         json.dumps(cr.unstructure())
    #     except:
    #         raise(ValueError(f"Could not serialize cr: {cr}"))
    # try:
    #     json.dumps(crs_container['best_clusterings'])
    # except:
    #     raise(ValueError(f"Could not serialize best_clusterings: {crs_container['best_clusterings']}"))
    # try:
    #     json.dumps(crs_container['connecting_reads'])
    # except:
    #     raise(ValueError(f"Could not serialize connecting_reads: {crs_container['connecting_reads']}"))

    # each cr can be serialized
    return json.dumps({
        "crs": [cr.unstructure() for cr in crs_container["crs"]],
        "connecting_reads": crs_container["connecting_reads"],
    })


def create_containers_db(path_db: Path, timeout: float):
    assert timeout > 0.0, f"timeout must be > 0.0. It is {timeout}"
    assert type(path_db) == PosixPath or type(path_db) == str, (
        f"path_db must be a Path or str. It is {type(path_db)}"
    )

    if path_db.exists():
        log.warning(f"Database {path_db} exists. Overwriting it.")
        path_db.unlink()

    conn = sqlite3.connect(path_db)
    conn.execute(
        """CREATE TABLE IF NOT EXISTS containers
        (crID INTEGER PRIMARY KEY, data TEXT)"""
    )
    conn.execute(f"pragma busy_timeout={str(int(timeout * 1000))}")
    conn.commit()
    conn.close()


def write_containers_to_db(path_db: Path, crs_containers: list[dict[str, object]]):
    assert type(path_db) == PosixPath or type(path_db) == str, (
        f"path_db must be a Path or str. It is {type(path_db)}"
    )
    assert path_db.exists(), f"Database {path_db} does not exist"
    assert path_db.is_file(), f"Database {path_db} is not a file"
    assert type(crs_containers) == list, "crs_containers is not a list"
    assert all([type(crs_container) == dict for crs_container in crs_containers]), (
        "crs_containers is not a list of dicts"
    )

    conn = sqlite3.connect(path_db)
    c = conn.cursor()
    pairs = [
        (int(crs_container["crID"]), serialize_crs_container(crs_container))
        for crs_container in tqdm(crs_containers)
    ]
    c.executemany("INSERT INTO containers VALUES (?,?)", pairs)
    conn.commit()
    c.close()
    conn.close()


# %%


def get_crs_containers(
    crs_dict: dict[int, datatypes.CandidateRegion],
    min_connections: int = 3,
    bnd_border_tolerance: int = 60,
    bnd_min_size: int = 500,
) -> dict[str, list[int]]:
    """computes a dict of the form readname -> [crID] for all reads that connect two or more crs via BND signals"""
    """candidate regions are only connected if they share at least [min_connections] reads"""

    dict_connectors: dict[str, list[int]] = {}
    for cr in crs_dict.values():
        for signal in cr.sv_signals:
            # BND left
            if (
                signal.sv_type == 3
                and signal.ref_start >= cr.referenceStart + bnd_border_tolerance
                and signal.size >= bnd_min_size
            ):
                if signal.readname not in dict_connectors:
                    dict_connectors[signal.readname] = []
                dict_connectors[signal.readname].append(cr.crID)
            # BND right
            if (
                signal.sv_type == 4
                and signal.ref_start <= cr.referenceEnd - bnd_border_tolerance
                and signal.size >= bnd_min_size
            ):
                if signal.readname not in dict_connectors:
                    dict_connectors[signal.readname] = []
                dict_connectors[signal.readname].append(cr.crID)

    dict_connections: dict[tuple[int, int], list[str]] = {}
    for readname, crIDs in dict_connectors.items():
        if len(set(crIDs)) > 1:
            for crIDa, crIDb in combinations(crIDs, 2):
                if (crIDa, crIDb) not in dict_connections:
                    dict_connections[(crIDa, crIDb)] = []
                dict_connections[(crIDa, crIDb)].append(readname)
    dict_connections_filtered = {
        key: val for key, val in dict_connections.items() if len(val) >= min_connections
    }
    # all keys (crID tuples) in dict_connections_filtered are connected crIDs.
    # flat set of all connected crIDs:
    set_connected_crIDs: set[int] = set([
        int(crID) for crIDs in dict_connections_filtered.keys() for crID in crIDs
    ])
    set_singleton_crIDs: set[int] = set(list(crs_dict.keys())) - set_connected_crIDs

    # create connected components
    UF = datastructures.UnionFind(sorted(set_connected_crIDs))
    for crIDa, crIDb in dict_connections_filtered.keys():
        UF.union_by_name(crIDa, crIDb)
    CC = UF.get_connected_components()

    log.info("Creating crs containers")

    crs_containers: list[dict[str, object]] = []
    for crID in set_singleton_crIDs:
        crs_containers.append({
            "crID": crID,
            "crs": [crs_dict[crID]],
            "connecting_reads": dict(),
        })

    for cc in CC:
        # find the connecting reads. They need to be in dict_connections_filtered
        connecting_reads: dict[str, list[int]] = dict()
        for crIDa, crIDb in combinations(cc, 2):
            if (crIDa, crIDb) in dict_connections_filtered:
                for readname in dict_connections_filtered[(crIDa, crIDb)]:
                    if readname not in connecting_reads:
                        connecting_reads[readname] = []
                    connecting_reads[readname].append(crIDa)
                    connecting_reads[readname].append(crIDb)
        # make each value in connecting_reads unique
        for readname in connecting_reads.keys():
            connecting_reads[readname] = list(set(connecting_reads[readname]))
        crID_representative = min(cc)
        crs_containers.append({
            "crID": crID_representative,
            "crs": [crs_dict[crID] for crID in cc],
            "connecting_reads": connecting_reads,
        })

    return crs_containers


def crs_containers_QC(
    crs_containers: list[dict[str, object]], path_histogram: Path
) -> None:
    if not str(path_histogram).endswith(".png"):
        raise ValueError(f"path_histogram must end with .png. It is {path_histogram}")
    # validate crs_containers
    assert type(crs_containers) == list, "crs_containers is not a list"
    assert all([type(crs_container) == dict for crs_container in crs_containers]), (
        "crs_containers is not a list of dicts"
    )
    # check if all crs_containers have the keys 'crID', 'crs' and 'connecting_reads'
    assert all(["crs" in crs_container for crs_container in crs_containers]), (
        "crs_containers does not contain the key 'crs'"
    )

    # create a histogram of the number of crs per crs container
    data = [len(cr["crs"]) for cr in crs_containers]
    # import seaborn as sns
    # # # create histogram with seaborn
    # sns.histplot(data,bins=range(1,max(data)+1),discrete=True)
    # import matplotlib.pyplot as plt
    # plt.title('Histogram of crs per container')
    # plt.xlabel('Number of crs')
    # plt.ylabel('Number of containers')
    # plt.grid(axis='y')
    # plt.savefig(path_histogram)


def crs_to_crs_containers(
    path_crs: Path,
    path_db: Path,
    QC_histogram: Path | None = None,
    crIDs_file: Path | None = None,
    min_connections: int = 3,
    bnd_border_tolerance: int = 60,
    bnd_min_size: int = 500,
) -> list[dict]:
    crs_dict = {
        cr.crID: cr for cr in signalstrength_to_crs.load_crs_from_db(path_db=path_crs)
    }
    crs_containers = get_crs_containers(
        crs_dict=crs_dict,
        min_connections=min_connections,
        bnd_border_tolerance=bnd_border_tolerance,
        bnd_min_size=bnd_min_size,
    )
    # write crs_containers to sqlite3 database
    # to do so, the containers need to be json serialized
    log.info(f"Creating database {path_db}")
    create_containers_db(path_db=path_db, timeout=10.0)
    log.info(f"writing crs containers to {path_db}")
    write_containers_to_db(path_db=path_db, crs_containers=crs_containers)
    crIDs = [cr["crID"] for cr in crs_containers]
    if QC_histogram:
        log.info(f"Creating QC histogram at {QC_histogram}")
        crs_containers_QC(crs_containers=crs_containers, path_histogram=QC_histogram)
    # write all crIDs of all containers to newline separated text file
    log.info(f"Writing crIDs to {crIDs_file}")
    if crIDs_file:
        with open(crIDs_file, "w") as f:
            for crID in crIDs:
                print(crID, file=f)


# %%


def run(args, **kwargs):
    crs_to_crs_containers(
        path_crs=args.input,
        path_db=args.database,
        QC_histogram=args.QC_histogram,
        crIDs_file=args.crIDs,
        min_connections=args.min_connections,
        bnd_border_tolerance=args.bnd_border_tolerance,
        bnd_min_size=args.bnd_min_size,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="creates a bed file with columns chr, start, end, crID."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to candidate regions database file.",
    )
    parser.add_argument(
        "-o",
        "--database",
        type=Path,
        required=True,
        help="Path to sqlite3 database that will contain all crs containers.",
    )
    parser.add_argument(
        "--QC_histogram",
        type=Path,
        required=False,
        default=None,
        help="Path to QC histogram file. Can be html or png.",
    )
    parser.add_argument(
        "--crIDs",
        type=Path,
        required=False,
        default=None,
        help="Optional path to txt file that lists all minimal crIDs. That is, if a container contains more than one candidate region, its representaative (the min) is reported. Can be used in snakemake checkpoints.",
    )
    parser.add_argument(
        "--min-connections",
        type=int,
        required=False,
        default=3,
        help="Minimum number of connections between two crs to be considered connected.",
    )
    parser.add_argument(
        "--bnd-border-tolerance",
        type=int,
        required=False,
        default=60,
        help="Minimum distance from the border of a cr to be considered a BND signal.",
    )
    parser.add_argument(
        "--bnd-min-size",
        type=int,
        required=False,
        default=500,
        help="Minimum size of a BND signal to be considered.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
