# # run alignment_to_raf for each bam file and cat all raf files to output
# import argparse
# import tempfile
# import typing
# from shlex import split
# import subprocess
# import logzero
# from logzero import logger
# from pathlib import Path
# from glob import glob

# from . import alignments_to_rafs, util

# def bams_to_raf(
#         reference:Path,
#         sampledicts:typing.List[dict],
#         min_singal_size:int,
#         min_mapq:int,
#         regions:Path,
#         threads:int,
#         output:Path,
#         blacklist:Path):
#     # create temporary directory in which all temp raf files are stored
#     tmp_dir = tempfile.TemporaryDirectory()
#     # for each bam file, create a tmp raf file in tmp_dir and run
#     # collect_signal_ReadAlignmentFragments on it
#     raf_tmp_files = []
#     dict_samplename_to_path,dict_samplename_to_id = util.load_sampledicts(sampledicts)[:2]
#     logger.info("collecting for all samples")
#     for samplename,bam in dict_samplename_to_path.items():
#         logger.info("  collecting for %s -- %s", samplename, bam)
#         tmp_raf = tempfile.NamedTemporaryFile(dir=tmp_dir.name,delete=False,suffix=".raf",prefix=samplename+'.')
#         raf_tmp_files.append(tmp_raf.name)
#         alignments_to_rafs.collect_signal_ReadAlignmentFragments(
#             reference=reference,
#             alignments=bam,
#             sampleID=dict_samplename_to_id[samplename],
#             min_singal_size=min_singal_size,
#             min_mapq=min_mapq,
#             regions=regions,
#             blacklist=blacklist,
#             threads=threads,
#             output=Path(tmp_raf.name))
#     # cat all tmp raf files to output
#     util.concat_txt_files(raf_tmp_files,output)
#     tmp_dir.cleanup()

# def run(args,**kwargs):
#     bams_to_raf(
#         reference=args.reference,
#         sampledicts=args.sampledicts,
#         min_singal_size=args.min_singal_size,
#         min_mapq=args.min_mapq,
#         regions=args.regions,
#         threads=args.threads,
#         output=args.output,
#         blacklist=args.blacklist)

# def get_parser():
#     parser = argparse.ArgumentParser(description='Extract read alignment fragment objects from provided bam files.')
#     parser.add_argument("-r","--reference", type=Path,required=True, help="Path to indexed reference.fa (+.fai)")
#     parser.add_argument("-d","--sampledicts", type=Path,required=True, help="Path to sample dictionaries in json format.")
#     parser.add_argument("-o","--output", type=Path,required=True, help="Path to output table to be sorted")
#     parser.add_argument("-m","--min_singal_size", type=int,required=False, default=6, help="Minimum intra alignment even (ins,del) to be reportet as a signal.")
#     parser.add_argument("-c","--min_clipped_length", type=int,required=False, default=150, help="Minimum bp of a clipped tail to be considered a break end.")
#     parser.add_argument("-q","--min_mapq", type=int,required=False, default=30, help="Minimum map quality to be used for raf extraction.")
#     parser.add_argument("-t","--threads", type=int,required=False, default=3, help="Number of threads for parallel processing.")
#     parser.add_argument("-g","--regions", type=Path,required=False, default=Path(""), help="Path to regions file (BED). Defaults to the whole reference genome.")
#     parser.add_argument("-b","--blacklist", type=Path,required=False, default=Path(""), help="Path to blacklist file (BED).")
#     return parser

# def main():
#     logzero.loglevel(logzero.DEBUG)
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()
