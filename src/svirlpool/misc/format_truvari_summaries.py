# this script takes as input a list of truvari summary.json files
# and a list of the same length of names for each summary
# then creates a table from the summary files
# and outputs it as a tsv file

import json

import pandas as pd


def format_truvari_summaries(
    summary_files: list[str],
    names: list[str],
    output_full: str | None = None,
    output_mini: str | None = None,
    output_gt_matrix: str | None = None,
    print_to_console: bool = True,
) -> None:
    # read in all summary files
    if len(summary_files) != len(names):
        raise ValueError("summary_files and names must be of the same length")
    summaries = []
    for f, name in zip(summary_files, names):
        with open(f) as fh:
            summary = json.load(fh)
            summary["name"] = name
            summaries.append(summary)
    # create a dataframe from the summaries
    df = pd.DataFrame(summaries)
    # reorder columns to have name first
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index("name")))
    df = df.reindex(columns=cols)
    # print the dataframe to console
    if print_to_console:
        print(df.to_string(index=False))
    # output the full dataframe to a tsv file
    if output_full is not None:
        df.to_csv(output_full, sep="\t", index=False)
    # output a mini dataframe to a tsv file (only name, precision, recall, f1)
    if output_mini is not None:
        df_mini = df[
            [
                "name",
                "precision",
                "recall",
                "f1",
                "FP",
                "FN",
                "TP-base",
                "gt_concordance",
            ]
        ]
        df_mini.to_csv(output_mini, sep="\t", index=False)
    # output the gt_matrix to a tsv file
    if output_gt_matrix is not None:
        with open(output_gt_matrix, "w") as fh:
            for i, row in df.iterrows():
                fh.write(f">{row['name']}\n")
                gt_matrix = row["gt_matrix"]
                # get all unique keys
                keys = set()
                for k1 in gt_matrix.keys():
                    keys.add(k1)
                    for k2 in gt_matrix[k1].keys():
                        keys.add(k2)
                keys = sorted(keys)
                # write header
                fh.write("\t" + "\t".join(keys) + "\n")
                # write rows
                for k1 in keys:
                    fh.write(k1)
                    for k2 in keys:
                        count = gt_matrix.get(k1, {}).get(k2, 0)
                        fh.write(f"\t{count}")
                    fh.write("\n")
                fh.write("\n")


import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Format truvari summary.json files into a table"
    )
    parser.add_argument(
        "-i", "--summary-files", nargs="+", help="List of truvari summary.json files"
    )
    parser.add_argument(
        "-n",
        "--names",
        nargs="+",
        help="List of names for each summary file (same length as summary_files)",
    )
    parser.add_argument(
        "--output-full", type=str, help="Output full table to this tsv file"
    )
    parser.add_argument(
        "-o",
        "--output-mini",
        type=str,
        help="Output mini table (name, precision, recall, f1) to this tsv file",
    )
    parser.add_argument(
        "--output-gt-matrix", type=str, help="Output gt_matrix to this tsv file"
    )
    parser.add_argument(
        "--no-console", action="store_true", help="Do not print the table to console"
    )
    args = parser.parse_args()

    format_truvari_summaries(
        summary_files=args.summary_files,
        names=args.names,
        output_full=args.output_full,
        output_mini=args.output_mini,
        output_gt_matrix=args.output_gt_matrix,
        print_to_console=not args.no_console,
    )
