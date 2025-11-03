# %%
# this script loads a signals file and a rafs file and a crs file
# subsets all data to the provided interval and prints the results
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import util


# %%
def print_signals(signals, region):
    # subset signals
    # sort signals by columns 0,1,2
    mask = signals[0] == region[0]
    df = signals.loc[mask, :]
    mask = df[1] >= region[1]
    mask &= df[1] <= region[2]
    df = df.loc[mask, :]
    # print each row of signals to terminal with a separating '\t'
    samples = df[7].unique()
    # for _,row in df.iterrows():
    #     print(*row,sep='\t')
    # # print list of samples
    print("unique samplenames")
    for sample in samples:
        print(sample)
    return


def print_rafs(rafs, region):
    rafs = sorted(
        rafs,
        key=lambda x: (
            x.referenceID,
            x.referenceAlignmentStart,
            x.referenceAlignmentEnd,
        ),
    )
    np_raf_regions = np.array(
        list(
            map(
                lambda x: [
                    x.referenceID,
                    x.referenceAlignmentStart,
                    x.referenceAlignmentEnd,
                ],
                rafs,
            )
        )
    )
    mask = np_raf_regions[:, 0] == region[0]
    mask &= np_raf_regions[:, 1] >= region[1]
    mask &= np_raf_regions[:, 2] <= region[2]
    # add rafs that overlap with the region
    mask2 = np_raf_regions[:, 0] == region[0]
    mask2 &= np_raf_regions[:, 1] <= region[1]
    mask2 &= np_raf_regions[:, 2] >= region[2]
    samplenames = []
    for index in np.where(mask | mask2)[0]:
        raf = rafs[index]
        samplenames.append(raf.read_name)
        # print(raf.unstructure())
        # print("----------")
    # print unique sample anmes
    print("unique samplenames")
    for sample in set(samplenames):
        print(sample)


def print_crs(crs, region):
    crs = sorted(crs, key=lambda x: (x.referenceID, x.referenceStart, x.referenceEnd))
    np_crs_regions = np.array(
        list(map(lambda x: [x.referenceID, x.referenceStart, x.referenceEnd], crs))
    )
    # also check for overlapping crs, not only for crs that are completely in the region
    # 1) check for crs that are completely in the region
    mask = np_crs_regions[:, 0] == region[0]
    mask &= np_crs_regions[:, 1] >= region[1]
    mask &= np_crs_regions[:, 2] <= region[2]
    # 2) check for crs that are overlapping with the region
    mask2 = np_crs_regions[:, 0] == region[0]
    mask2 &= np_crs_regions[:, 1] <= region[1]
    mask2 &= np_crs_regions[:, 2] >= region[2]
    for index in np.where(mask | mask2)[0]:
        cr = crs[index]
        print(
            cr.chr,
            cr.referenceStart,
            cr.referenceEnd,
            cr.referenceEnd - cr.referenceStart,
            cr.crID,
            sep="\t",
        )
        for samplename in set([row[7] for row in cr.sv_signals]):
            print(samplename)
        print("----------")


def print_signal_processing(
    region: str,
    input_signals: Path,
    input_raf: Path,
    input_crs: Path,
    path_reference: Path,
):
    signals = pd.read_csv(input_signals, sep="\t", header=None, dtype={9: str})
    rafs = list(util.yield_from_raf(input_raf))
    crs = list(util.yield_from_crs(input_crs))
    refdict = util.create_ref_dict(reference=path_reference)
    rrefdict = {v: k for k, v in refdict.items()}
    region = region.replace(",", "").replace(":", "-").split("-")
    chr = region[0]
    region[0] = rrefdict[region[0]]
    region[1] = int(region[1])
    region[2] = int(region[2])
    print("printing signal processing for region", chr, region[1], region[2])
    print(
        "--------------------------------------------------------------------------------"
    )
    print("signals")
    print_signals(signals, region)
    print(
        "--------------------------------------------------------------------------------"
    )
    print("rafs")
    print_rafs(rafs, region)
    print(
        "--------------------------------------------------------------------------------"
    )
    print("crs")
    print_crs(crs, region)


# %%


def run(args, **kwargs):
    print_signal_processing(
        region=args.region,
        input_signals=args.signals,
        input_raf=args.rafs,
        input_crs=args.crs,
        path_reference=args.reference,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Computes the consensus objects of a given fastq file of cut reads for a certain candidate region."
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        required=True,
        help="Region in the format chr:start-end.",
    )
    parser.add_argument(
        "-s", "--signals", type=Path, required=True, help="Path to the signals file."
    )
    parser.add_argument(
        "-a", "--rafs", type=Path, required=True, help="Path to the rafs file."
    )
    parser.add_argument(
        "-c", "--crs", type=Path, required=True, help="Path to the crs file."
    )
    parser.add_argument(
        "-f",
        "--reference",
        type=Path,
        required=True,
        help="Path to the reference file.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()

# # %%
# input_signals = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.signal")
# input_raf = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.raf")
# input_crs = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.crs")
# path_reference = Path("/home/vinzenz/development/LRSV-detection/development/test/hs37d5.fa")

# region = "1:82,718,748-82,723,465"
# signals = pd.read_csv(input_signals,sep='\t',header=None,dtype={9:str})
# rafs = list(util.yield_from_raf(input_raf))
# crs = list(util.yield_from_crs(input_crs))
# refdict = util.create_ref_dict(reference=path_reference)
# rrefdict = {v:k for k,v in refdict.items()}
# region = region.replace(',','').replace(':','-').split('-')
# chr = region[0]
# region[0] = rrefdict[region[0]]
# region[1] = int(region[1])
# region[2] = int(region[2])
# # %%
# # subset signals
# # sort signals by columns 0,1,2
# print_signals(signals,region)
# print_rafs(rafs,region)
# print_crs(crs,region)
# %%
