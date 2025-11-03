# #%%
# # this script loads a table of datatypes.InDelBndCandidate objects that have been unstructured and
# # written via json dump to a tsv with columns refChr, refStart, InDelBndCandidate
# # then it merges the InDelBndCandidate objects into svCandidate objects (if they can be merged)
# # and loads ref and alt sequences from the provided consensus fasta file and
# # the provided reference fasta file
# # finally it writes the svCandidate objects to a tsv file
# #%%
# import argparse
# import typing
# import csv
# import json
# import cattrs
# from pathlib import Path
# from Bio import SeqIO
# import tempfile
# import subprocess
# from logzero import logger as log
# import numpy as np
# from tqdm import tqdm

# from . import datatypes, util

# # TODO: unit test
# def ref_sequences_from_svs(
#         IDBCs:typing.List[datatypes.InDelBndCandidate],
#         path_reference:Path,
#         max_region_size:int) -> typing.List[str]:
#     log.info("Adding ref sequences to SVs..")
#     # bed table from svs
#     bed_table = np.array([[idbc.refChr,int(idbc.refStart),int(idbc.refEnd) if idbc.svType=='DEL' else int(idbc.refEnd)+1] for idbc in IDBCs],dtype=object)
#     # adjust bed table for regions that are too large. If a region is too large, the end coordinate is set to start+1
#     # WARNING! Any ref==alt field need to be processed later. alt and ref can't remain the same according to VCF standard.
#     region_sizes = bed_table[:,2] - bed_table[:,1]
#     mask = region_sizes > max_region_size
#     bed_table[mask,2] = bed_table[mask,1]
#     log.info(f"{sum(mask)} regions were too large. The ref sequences of these regions are limited to size 1.")
#     # save bed_table as temporary bed file
#     tmp_bed = tempfile.NamedTemporaryFile(mode='w',delete=False,suffix='.bed')
#     with open(tmp_bed.name,'w') as f:
#         for line in bed_table:
#             print('\t'.join(map(str,line)),sep='\t',file=f)
#     # run bedtools getfasta
#     cmd_getfasta = ['bedtools','getfasta','-fi',str(path_reference),'-bed',tmp_bed.name]
#     tmp_alt_fasta = tempfile.NamedTemporaryFile(mode='w',delete=False,suffix='.fasta')
#     with open(tmp_alt_fasta.name,'w') as f:
#         log.info(f"extracting ref sequences into {tmp_alt_fasta.name}...")
#         subprocess.check_call(cmd_getfasta,stdout=f)
#     # load fasta file
#     ref_sequences = list(map(lambda x: str(x.seq),list(SeqIO.parse(tmp_alt_fasta.name, "fasta"))))
#     return ref_sequences

# # TODO: unit test
# def get_alt_sequences(
#         IDBCs:typing.List[datatypes.InDelBndCandidate],
#         path_consensus_sequences:Path) -> typing.List[str]:
#     consensus_sequences = {record.id:str(record.seq) for record in SeqIO.parse(path_consensus_sequences,"fasta")}
#     sequences = []
#     for idbc in IDBCs:
#         if idbc.consensus_name not in consensus_sequences:
#             raise ValueError(f"consensus name {idbc.consensus_name} not found in consensus fasta file {str(path_consensus_sequences)}.")
#         if idbc.svType == 'INS':
#             if bool(idbc.is_reverse):
#                 sequences.append(''.join(util.reverse_complement(list(consensus_sequences[idbc.consensus_name][idbc.consensus_start:idbc.consensus_end]))))
#             else:
#                 sequences.append(consensus_sequences[idbc.consensus_name][idbc.consensus_start:idbc.consensus_end])
#         else:
#             sequences.append(consensus_sequences[idbc.consensus_name][idbc.consensus_start])
#     return sequences


# # TODO: unit test
# def preprocess_IDBCs(
#         path_idbc_table:Path,
#         path_consensus_sequences:Path,
#         path_reference:Path,
#         path_output:Path,) -> None:

#     # 1) load the table of InDelBndCandidate objects
#     log.info("loading IDBCs..")
#     IDBCs = list(util.yield_IDBCs(path_idbc_table))

#     # TODO: optimize to avoid loading all sequences to memory.
#     log.info("loading consensus sequences..")
#     alt_sequences = get_alt_sequences(IDBCs=IDBCs,path_consensus_sequences=path_consensus_sequences)

#     log.info(f"extracting ref sequences..")
#     ref_sequences = ref_sequences_from_svs(IDBCs,path_reference,1_000_000)

#     log.info(f"writing to {path_output}..")
#     # write to a new tsv file of the format refChr, refStart, idbc, ref, alt
#     with open(path_output,"w") as f:
#         writer = csv.writer(f,delimiter='\t',quotechar='"',quoting=csv.QUOTE_MINIMAL)
#         for idbc,alt,ref in zip(IDBCs,alt_sequences,ref_sequences):
#             writer.writerow([idbc.refChr,idbc.refStart,json.dumps(cattrs.unstructure(idbc)),ref,alt])
#     log.info("done.")
# # %%

# # argparse block
# def run(args,**kwargs):
#     preprocess_IDBCs(
#         path_idbc_table=args.input,
#         path_consensus_sequences=args.consensuses,
#         path_reference=args.reference,
#         path_output=args.output)

# def get_parser():
#     parser = argparse.ArgumentParser(description='creates a vcf file from a provided svs file.')
#     parser.add_argument("-i", "--input", type=Path,required=True, help="Path to .IDBCs.tsv file.")
#     parser.add_argument("-c", "--consensuses", type=Path,required=True, help="Path to consensus fasta file.")
#     parser.add_argument("-r", "--reference", type=Path,required=True, help="Path to reference fasta file.")
#     parser.add_argument("-o", "--output", type=Path,required=True, help="Path to output annotated aIDBCs.tsv file.")
#     return parser

# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()
