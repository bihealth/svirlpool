# %%
import argparse
import csv
import json
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from . import datatypes, util


# %%
# --- use for optimization --- #
def get_readnames_of_CandidateRegion(cr: datatypes.CandidateRegion) -> set:
    return set([sn[8] for sn in cr.sv_signals])


# iterate all crs and save their read names to a dict of the form
# dict[str,{int}] or dict[readname:{crIDs}]
def create_readname_to_crID_dict(path_crs: Path) -> dict[str, set[int]]:
    dict_readname_to_crID = {}
    for cr in util.yield_from_crs(input=path_crs):
        for readname in get_readnames_of_CandidateRegion(cr):
            if readname in dict_readname_to_crID:
                dict_readname_to_crID[readname].add(cr.crID)
            else:
                dict_readname_to_crID[readname] = {cr.crID}
    return dict_readname_to_crID


# for each readname in a candidate region, find all other candidate regions that have the same readname
# then check if they are connected by BND or open and closed DELs


# %%


def find_open_DELLs(a: datatypes.CandidateRegion) -> set:
    """returns all readIDs that have an open DELL, which means that they are an open left end of a large DEL."""
    expected_reads = set([sn[8] for sn in a.sv_signals if sn[5] == 1])
    open_reads = []
    used_indeces = []
    for i in range(len(a.sv_signals) - 1, -1, -1):
        svtype = a.sv_signals[i][5]
        readID = a.sv_signals[i][8]
        if svtype == 1:  # DELL
            if readID in expected_reads:
                open_reads.append(readID)
                used_indeces.append(i)
        if svtype == 2:  # DELR
            if readID in expected_reads:
                expected_reads.remove(readID)
    return set(open_reads)


def find_open_DELRs(a: datatypes.CandidateRegion) -> set:
    """returns all readIDs that have an open DELR, which means that they are an open right end of a large DEL."""
    expected_reads = set([sn[8] for sn in a.sv_signals if sn[5] == 2])
    open_reads = []
    used_indeces = []
    for i in range(len(a.sv_signals)):
        sv_type = a.sv_signals[i][5]
        readID = a.sv_signals[i][8]
        if sv_type == 2:  # DELR
            if readID in expected_reads:
                open_reads.append(readID)
                used_indeces.append(i)
        if sv_type == 1:  # DELL
            if readID in expected_reads:
                expected_reads.remove(readID)
    return set(open_reads)


def find_open_BNDs(a: datatypes.CandidateRegion, minimum_size: int = 100) -> set:
    """returns all readIDs that have a BND signal in a."""
    # sn[13] == -1 excludes all BNDs that are in a repeat
    return set(
        [
            sn[8]
            for sn in a.sv_signals
            if sn[5] >= 3 and sn[6] >= minimum_size and sn[13] == -1
        ]
    )


def indeces_open_DELLs(a: datatypes.CandidateRegion) -> list:
    """returns all readIDs that have an open DELL, which means that they are an open left end of a large DEL."""
    expected_reads = set([sn[8] for sn in a.sv_signals if sn[5] == 1])
    used_indeces = []
    for i in range(len(a.sv_signals) - 1, -1, -1):
        svtype = a.sv_signals[i][5]
        readID = a.sv_signals[i][8]
        if svtype == 1:  # DELL
            if readID in expected_reads:
                used_indeces.append(i)
        if svtype == 2:  # DELR
            if readID in expected_reads:
                expected_reads.remove(readID)
    return used_indeces


def indeces_open_DELRs(a: datatypes.CandidateRegion) -> list:
    """returns all readIDs that have an open DELR, which means that they are an open right end of a large DEL."""
    expected_reads = set([sn[8] for sn in a.sv_signals if sn[5] == 2])
    used_indeces = []
    for i in range(len(a.sv_signals)):
        sv_type = a.sv_signals[i][5]
        readID = a.sv_signals[i][8]
        if sv_type == 2:  # DELR
            if readID in expected_reads:
                used_indeces.append(i)
        if sv_type == 1:  # DELL
            if readID in expected_reads:
                expected_reads.remove(readID)
    return used_indeces


def indeces_open_BNDs(a: datatypes.CandidateRegion) -> list:
    """returns all readIDs that have a BND signal in a."""
    open_bnds = find_open_BNDs(a)
    return [i for i, sn in enumerate(a.sv_signals) if sn[5] >= 4 and sn[8] in open_bnds]


def connect_crs(
    input: Path,
    minimum_connecting_reads: int,
    minimum_connecting_coverage: float,
    minimum_abs_coverage: int = 3,
) -> dict:
    crs = list(util.yield_from_crs(input=input))
    # dict_readIDs_to_sampleIDs = {}
    # for cr in crs:
    #     for svc in cr.sv_signals:
    #         dict_readIDs_to_sampleIDs[svc[7]] = svc[8]
    # for each cr in crs, create three sets of open_BNDs, open_DELLs, open_DELRs
    open_svs = {}
    for cr in crs:
        x = [find_open_BNDs(cr), find_open_DELLs(cr), find_open_DELRs(cr)]
        if len(x[0]) + len(x[1]) + len(x[2]) >= minimum_connecting_reads:
            open_svs[cr.crID] = x
    # find all connecting readIDs between all pairs of crs in open_svs.keys()
    # a connection is a key tuple (crID_a,crID_b) with a set of readIDs as value
    connections = {}
    lk = list(open_svs.keys())
    for i in tqdm(range(len(open_svs.keys()))):
        ki = lk[i]
        for j in range(i + 1, len(open_svs.keys())):
            kj = lk[j]
            # connect BND with BND
            # connect DELL in i with DELR in j
            connection = [
                open_svs[ki][0].intersection(open_svs[kj][0]),
                open_svs[ki][1].intersection(open_svs[kj][2]),
            ]
            sum_connection = len(set(connection[0]).union(set(connection[1])))
            if sum_connection > 0:
                # filter out connections that do not have enough coverage or
                # alt reads compared to total reads
                df_ki = pd.DataFrame(crs[ki].sv_signals)
                df_kj = pd.DataFrame(crs[kj].sv_signals)
                ind_ki = indeces_open_DELLs(crs[ki])
                ind_ki.extend(indeces_open_BNDs(crs[ki]))
                ind_kj = indeces_open_DELRs(crs[kj])
                ind_kj.extend(indeces_open_BNDs(crs[kj]))
                coverage_ki = (
                    df_ki.iloc[ind_ki, :]
                    .loc[:, [5, 9, 7]]
                    .groupby([9])
                    .mean()
                    .iloc[:, 1]
                    .to_numpy()
                )
                coverage_alt_ki = (
                    df_ki.iloc[ind_ki, [9, 0]]
                    .groupby([9])
                    .count()
                    .iloc[:, 0]
                    .to_numpy()
                )
                coverage_kj = (
                    df_kj.iloc[ind_kj, :]
                    .loc[:, [5, 9, 7]]
                    .groupby([9])
                    .mean()
                    .iloc[:, 1]
                    .to_numpy()
                )
                coverage_alt_kj = (
                    df_kj.iloc[ind_kj, [9, 0]]
                    .groupby([9])
                    .count()
                    .iloc[:, 0]
                    .to_numpy()
                )
                ki_relative_cov = coverage_alt_ki / coverage_ki
                kj_relative_cov = coverage_alt_kj / coverage_kj
                is_valid = (
                    (
                        (ki_relative_cov > minimum_connecting_coverage).any()
                        or (kj_relative_cov > minimum_connecting_coverage).any()
                    )
                    and coverage_ki.sum() > minimum_abs_coverage
                    and coverage_kj.sum() > minimum_abs_coverage
                )
                if is_valid:
                    connections[(ki, kj)] = connection
    return connections


# %%
def crs_to_connections(
    input: Path,
    output: Path,
    minimum_connecting_reads: int,
    minimum_connecting_coverage: float,
    minimum_abs_coverage: int,
) -> None:
    """generates and saves to output a dictionary of connections between crs in input."""
    connections = connect_crs(
        input=input,
        minimum_abs_coverage=minimum_abs_coverage,
        minimum_connecting_reads=minimum_connecting_reads,
        minimum_connecting_coverage=minimum_connecting_coverage,
    )
    # save connections to file
    connections_remapped = [[k, [list(w) for w in v]] for k, v in connections.items()]
    with open(output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for key, value in connections_remapped:
            writer.writerow([json.dumps(key), json.dumps(value)])


def run(args, **kwargs):
    crs_to_connections(
        input=args.input,
        output=args.output,
        minimum_abs_coverage=args.minimum_abs_coverage,
        minimum_connecting_coverage=args.minimum_connecting_coverage,
        minimum_connecting_reads=args.minimum_connecting_reads,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Computes MCRs from CRs and does QC both for pre-split MCRs and post-split."
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Path to crs file."
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to output tsv mcrs file."
    )
    parser.add_argument(
        "--minimum_connecting_reads",
        type=int,
        required=False,
        default=1,
        help="Minimum number of connecting reads between two crs.",
    )
    parser.add_argument(
        "--minimum_connecting_coverage",
        type=float,
        required=False,
        default=0.15,
        help="Minimum relative coverage of connecting reads between two crs.",
    )
    parser.add_argument(
        "--minimum_abs_coverage",
        type=int,
        required=False,
        default=3,
        help="Minimum absolute coverage of all reads reads between two crs.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
