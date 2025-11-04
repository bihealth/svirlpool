# import argparse
# from pathlib import Path
# from tqdm import tqdm
# import gzip
# import tempfile
# from bisect import insort_left

# from . import alignments_to_rafs



# def split_DELs(
#         path_signaldepths: Path,
#         path_output: Path,
#         threads:int
#         ) -> None:
#     cache = []
#     last_chrID = None
#     tmp_output = tempfile.NamedTemporaryFile(suffix=".tmp.split_dels.bed",delete=True)
#     with open(tmp_output.name, 'w') as out:
#         for line in tqdm(gzip.open(path_signaldepths, 'rt')):
#             current_signal = line.strip().split('\t')
#             if last_chrID and last_chrID != current_signal[0]:
#                 # save all items in cache to output
#                 for cached_signal in cache:
#                     print(*cached_signal,sep='\t',file=out)
#                 cache = []
#             while len(cache)>0 and int(current_signal[1]) > int(cache[0][1]):
#                 cached_signal = cache.pop(0)
#                 print(*cached_signal,sep='\t',file=out)
#             if current_signal[3] == 1:
#                 # create DELL
#                 dell = current_signal.copy()
#                 dell[2] = str(int(dell[1])+1)
#                 dell[3] = str(1)
#                 print(*dell,sep='\t',file=out)
#                 delr = current_signal.copy()
#                 delr[1] = str(int(delr[2])-1)
#                 delr[3] = str(2)
#                 insort_left(a=cache,x=delr,key=lambda x: (int(x[0]),int(x[1])))
#                 # cache.append(delr)
#                 # cache = sorted(cache, key=lambda x: int(x[1]))
#             else:
#                 print(*current_signal,sep='\t',file=out)
#             last_chrID = current_signal[0]
#         if len(cache)>0:
#             for cached_signal in cache:
#                 print(*cached_signal,sep='\t',file=out)
#     alignments_to_rafs.compress_and_index_bedlike(
#         sort_numerically=True,
#         input=tmp_output.name,
#         output=path_output,
#         threads=threads)

# # %%


# def run(args,**kwargs):
#     split_DELs(
#         path_signaldepths=args.input,
#         path_output=args.output,
#         threads=args.threads,
#     )

# def get_parser():
#     parser = argparse.ArgumentParser(description='Split DELs (1) to DELLs (1) and DELRs (2) and save to output.')
#     parser.add_argument("-i","--input", type=Path,required=True, help="Path to the signaldepths tsv file.")
#     parser.add_argument("-o","--output", type=Path,required=True, help="Path to the output signals tsv file.")
#     parser.add_argument("-t","--threads", type=int,required=False, default=1, help="Number of threads to use.")
#     return parser

# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()
