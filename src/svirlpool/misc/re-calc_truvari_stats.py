# %%
import argparse
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path


# %%
def calc_new_stats(args, **kwargs) -> None:
    # overlap=0.8
    # wdir=Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/truvari/ont-r10.35x")
    with open(args.wdir / "diff.bed", "w") as f:
        cmd_intersect = f"bedtools intersect -f {args.overlap} -e -a {str(args.wdir / 'fp.bed')} -b {str(args.wdir / 'fn.bed')}"
        subprocess.check_call(shlex.split(cmd_intersect), stdout=f)

    D = defaultdict(int)
    for stat, metric in zip(
        ["diff", "tp-call", "fn", "fp"], ["delta", "tp", "fn", "fp"]
    ):
        cmd_wcl = f"wc -l {str(args.wdir / (stat+'.bed'))}"
        output = subprocess.check_output(shlex.split(cmd_wcl))
        n = int(str(output)[2:].split(" ")[0])
        D[metric] = n

    # re-calc metrics
    TP = D["tp"] + D["delta"]
    FP = D["fp"] - D["delta"]
    FN = D["fn"] - D["delta"]

    prec = TP / (TP + FP)
    sens = TP / (TP + FN)
    F1 = (2 * prec * sens) / (prec + sens)

    oPrec = D["tp"] / (D["tp"] + D["fp"])
    oSens = D["tp"] / (D["tp"] + D["fn"])
    oF1 = (2 * oPrec * oSens) / (oPrec + oSens)

    results = {
        n: [m, o]
        for m, n, o in zip(
            [F1, sens, prec, TP, FP, FN],
            ["F1", "sens", "prec", "TP", "FP", "FN"],
            [oF1, oSens, oPrec, D["tp"], D["fp"], D["fn"]],
        )
    }

    for k in results.keys():
        print(f"{k}\t{results[k][0]}\t{results[k][1]}")
    if args.outfile:
        with open(args.outfile, "w") as f:
            for k in results.keys():
                print(f"{k}\t{results[k][0]}\t{results[k][1]}", file=f)


# %%


def get_parser():
    parser = argparse.ArgumentParser(
        description="re-calculates the performance statistics of truvari output"
    )
    parser.add_argument(
        "--wdir", type=Path, required=True, help="Path to the truvari directory"
    )
    parser.add_argument(
        "--script",
        type=Path,
        required=False,
        default="/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/scripts/result_vcf_to_beds.sh",
        help="Path to the 'result_vcf_to_beds.sh' script.",
    )
    parser.add_argument(
        "--overlap",
        type=float,
        required=False,
        default=0.8,
        help="Necessary fraction of mutual overlap",
    )
    parser.add_argument(
        "--outfile",
        type=Path,
        required=False,
        help="Prints the new metrics to the provided file.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    calc_new_stats(args)
    return


if __name__ == "__main__":
    main()
