# #%%
# # load unfiltered signaldepths
# # load mosdepth output

# from pathlib import Path
# import argparse
# import subprocess
# from shlex import split
# import pandas as pd
# import numpy as np
# import plotly.express as px
# from . import util
# import tempfile
# import csv
# from tqdm import tqdm
# from collections import defaultdict

# from logzero import logger as log

# #%%
# # =============================================================================
# #  sketch the steps
# # =============================================================================

# def print_histograms(dict_all_counts:dict, dict_mean_cov:dict, path_figures:Path, cutoff_factor:float):
#     # create histograms from count dicts
#     # e.g. create df from counts dict and merge successive n cells
#     df_histogram_all = pd.DataFrame(dict_all_counts)
#     # merge successive n rows for each column
#     bin_width = int(round(np.log(df_histogram_all.sum().sum())*2))
#     df_histogram_all = df_histogram_all.groupby(np.arange(len(df_histogram_all))//bin_width).sum()
#     # change index to names of intervals (e.g. 0-29, 30-59, ...)
#     df_histogram_all.index = [f"{i*bin_width}-{(i+1)*bin_width-1}" for i in df_histogram_all.index]
#     # stack columns into one column
#     df_histogram_all = df_histogram_all.stack().reset_index()
#     df_histogram_all.columns=['bin','sample','count']
#     fig_bar_all_total = px.bar(df_histogram_all,x='bin',y='count',facet_row='sample')
#     fig_bar_all_total.write_html(path_figures / "histogram_signaldepth_cutoff_all.html")

#     # --- plot a zoomed in version without binning between 1 and maximum of (2*mean depth)
#     df_histogram_small = pd.DataFrame(dict_all_counts).iloc[:max([round(mean_cov*cutoff_factor) for mean_cov in dict_mean_cov.values()]),:]
#     df_histogram_small = df_histogram_small.stack().reset_index()
#     df_histogram_small.columns=['bin','sample','count']
#     fig_bar_all_small = px.bar(df_histogram_small,x='bin',y='count',facet_row='sample')
#     fig_bar_all_small.write_html(path_figures / "histogram_signaldepth_cutoff_passed.html")

# def get_mean_covs(path_sampledicts:Path,threads:int,override_depth:int) -> dict:
#     dict_ID_bam = util.load_sampledicts(path_sampledicts)[2]
#     # for each provided sample, calc depth and filter signaldepths
#     # if mosdepth summary is provided, use it.
#     # if not, then use a provided alignments (bam) file
#     # to call mosdepth and use its summary.
#     dict_mean_cov = {}
#     for ID,path_bam in dict_ID_bam.items():
#         path_mosdepth = Path(path_bam).with_suffix('.mosdepth.summary.txt')
#         if path_mosdepth.exists():
#             log.info(f"Found mosdepth summary file for {ID}: {str(path_mosdepth)}")
#             df_mosdepth = pd.read_csv(path_mosdepth, sep='\t')
#             mean_cov = float(df_mosdepth.loc[:,'mean'].iloc[-1])
#         elif override_depth > 0:
#             log.info(f"Did not find mosdepth summary file for {ID}. Using override depth: {override_depth}")
#             mean_cov = override_depth
#         else:
#             log.info(f"Did not find mosdepth summary file for {ID}. Creating new one. Calling mosdepth..")
#             #tmp_dir = tempfile.TemporaryDirectory()
#             #prefix = Path(tmp_dir.name) / Path(path_bam).stem
#             cmd_mosdepth = split(f"mosdepth -x -t {threads} -m {Path(path_bam).with_suffix('')} {path_bam}")
#             log.debug(f"cmd mosdepth: {' '.join(cmd_mosdepth)}")
#             subprocess.run(cmd_mosdepth)
#             # read mosdepth summary
#             df_mosdepth = pd.read_csv(path_mosdepth, sep='\t')
#             mean_cov = float(df_mosdepth.loc[:,'mean'].iloc[-1])
#             log.debug(f"found mean coverage: {mean_cov} of {ID}")
#         dict_mean_cov[ID] = mean_cov
#     return dict_mean_cov

# def filter_signaldepths(input:Path,output:Path,dict_mean_cov:dict,path_sampledicts:Path,cutoff_factor:float,min_coverage:int) -> tuple:
#     dict_ID_bam = util.load_sampledicts(path_sampledicts)[2]
#     dict_all_counts = {ID:defaultdict(int) for ID in dict_ID_bam.keys()}
#     log.info("filtering signaldepths..")
#     with open(input) as f:
#         with open(output, 'w') as f_out:
#             reader = csv.reader(f, delimiter='\t')
#             for line in tqdm(reader):
#                 ID = int(line[8])
#                 depth = int(line[11])
#                 dict_all_counts[ID][depth] += 1
#                 passed = depth <= round(cutoff_factor*dict_mean_cov[ID])
#                 passed = passed and depth >= min_coverage
#                 if passed:
#                     # write line to output file
#                     print(*line, sep='\t', file=f_out)
#     return dict_all_counts

# def filter_and_QC_signaldepths(
#             input:Path,
#             output:Path,
#             figures:Path,
#             sampledicts:Path,
#             threads:int,
#             cutoff_factor:float,
#             min_coverage:int,
#             override_depth:int):
#     dict_mean_cov = get_mean_covs(
#         path_sampledicts=sampledicts,
#         threads=threads,
#         override_depth=override_depth)
#     log.info(f"mean coverages: {dict_mean_cov}")
#     dict_all_counts = filter_signaldepths(
#         input=input,
#         output=output,
#         dict_mean_cov=dict_mean_cov,
#         path_sampledicts=sampledicts,
#         cutoff_factor=cutoff_factor,
#         min_coverage=min_coverage)
#     print_histograms(
#         dict_all_counts=dict_all_counts,
#         dict_mean_cov=dict_mean_cov,
#         path_figures=figures,
#         cutoff_factor=cutoff_factor)

# def run(args,**kwargs):
#     filter_and_QC_signaldepths(
#         input=args.input,
#         output=args.output,
#         sampledicts=args.sampledicts,
#         cutoff_factor=args.cutoff_factor,
#         threads=args.threads,
#         figures=args.figures,
#         min_coverage=args.min_coverage,
#         override_depth=args.override_depth)

# def get_parser():
#     parser = argparse.ArgumentParser(description="Filters signaldepths for maximum depth of coverage on a locus and creates QC plots.")
#     parser.add_argument('-i','--input',type=Path,required=True,help='Input file with signaldepths.')
#     parser.add_argument('-o','--output',type=Path,required=True,help='Output file with filtered signaldepths.')
#     parser.add_argument('-s','--sampledicts',type=Path,required=True,help='Sampledicts file.')
#     parser.add_argument('-t','--threads',type=int,required=True,help='Number of threads to use. Especially important if mosdepth requires to create the summary files in the location of the bam files.')
#     parser.add_argument('-c','--cutoff_factor',type=float,required=False,default=6.0,help='Factor to multiply the mean coverage with to get the cutoff for the maximum depth of coverage.')
#     parser.add_argument('-f','--figures',type=Path,required=True,help='basepath to save QC figures.')
#     parser.add_argument('-m','--min_coverage',type=int,required=False,default=2,help='Minimum coverage per sample to keep a SV signal.')
#     parser.add_argument('--override-depth',type=int,default=0,help='Override the depth of coverage with a fixed value. Useful for testing.')
#     return parser

# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()
