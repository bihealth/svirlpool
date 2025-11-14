# # %%

# import argparse
# from pathlib import Path

# import numpy as np
# import pandas as pd
# import plotly.express as px

# from . import datatypes, util

# # %%


# def cr_size(cr: datatypes.CandidateRegion) -> int:
#     return cr.referenceEnd - cr.referenceStart


# def n_signals(cr: datatypes.CandidateRegion) -> int:
#     return len(cr.sv_signals)


# def n_reads(cr: datatypes.CandidateRegion) -> int:
#     return len({s.aug_readname() for s in cr.sv_signals})


# def avg_coverage(cr: datatypes.CandidateRegion) -> float:
#     return float(np.mean([s.coverage for s in cr.sv_signals]))


# def n_repeats(cr: datatypes.CandidateRegion) -> int:
#     return len({s.repeatID for s in cr.sv_signals if s.repeatID != -1})


# # iterate crs and extract stats
# def create_df(input: Path) -> pd.DataFrame:
#     return pd.DataFrame(
#         [
#             [cr_size(cr), n_signals(cr), n_reads(cr), avg_coverage(cr), n_repeats(cr)]
#             for cr in util.yield_from_crs(input)
#         ],
#         columns=["cr_size", "n_signals", "n_reads", "avg_coverage", "n_repeats"],
#     )


# def plots(df: pd.DataFrame, output: Path) -> None:
#     # plot histograms of n_signals, n_reads, avg_coverage, n_repeats
#     for col in df.columns:
#         fig = px.histogram(df, x=col, nbins=100)
#         fig.write_html(output / f"QC_crs.{col}.html")


# def QC_crs(input: Path, output: Path) -> None:
#     df = create_df(input)
#     plots(df, output)


# def run(args, **kwargs):
#     QC_crs(input=args.input, output=args.output)


# def get_parser():
#     parser = argparse.ArgumentParser(
#         description="Convolves over the singal to get the localized signalstrength to be used to create candidate regions."
#     )
#     parser.add_argument(
#         "-i", "--input", type=Path, required=True, help="tsv file of candidate regions."
#     )
#     parser.add_argument(
#         "-o",
#         "--output",
#         type=Path,
#         required=True,
#         help="output directory for the plots of name scheme '[output]/QC_crs.{{colname}}.html.'",
#     )
#     return parser


# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return


# if __name__ == "__main__":
#     main()
