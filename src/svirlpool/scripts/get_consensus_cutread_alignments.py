# import argparse
# from pathlib import Path
# from Bio import SeqIO, SeqRecord
# from Bio.Seq import Seq
# import pysam
# import json
# import subprocess
# import shlex

# from logzero import logger as log

# from . import datatypes, consensus_containers_to_db, sv_calling, consensus_lib, util


# def parse_SequenceObject_to_BioSeq(sequenceObject:datatypes.SequenceObject) -> SeqRecord:
#     sr = SeqRecord.SeqRecord(
#         Seq(sequenceObject.sequence),
#         id=sequenceObject.name,
#         description="" if not sequenceObject.description else sequenceObject.description)
#     if sequenceObject.qualities and len(sequenceObject.qualities) == len(sequenceObject.sequence):
#         sr.letter_annotations["phred_quality"] = sequenceObject.qualities
#     return sr

# def write_consensus_sequence_and_alignments_to_files(
#         path_db:Path,
#         crID:int,
#         path_base_output:Path) -> list[str]:
#     """reads the consensus objects with matching crIDs from the database and gets all reconstructible_reads and writes those alignments to a bam file \
# the bam file is written to path_base_output+.crID+.bam and indexed. The consensus is written to path_base_output+crID+.fasta and indexed"""
#     path_base_output = Path(path_base_output)
#     # load all consensusIDs from the database that start with crID
#     log.info(f"loading all consensusIDs from the database {path_db}")
#     consensusIDs = consensus_containers_to_db.load_all_consensusIDs_from_database(path_database=path_db)
#     log.info(f"found {len(consensusIDs)} consensus objects in the database, names: {consensusIDs}")
#     # filter all consensusIDs that start with crID
#     log.info(f"filtering all consensusIDs that start with {crID}")
#     consensusIDs = [consensusID for consensusID in consensusIDs if consensusID.startswith(str(crID)+'.')]
#     # a consensus object has a consensusID and a list of reconstructible_reads
#     # load all consensus objects from the database
#     # create a dict consensusID:[reconstructed_reads]
#     log.info(f"processing all consensus objects with IDs: {consensusIDs} from the database {path_db}")
#     dict_alignments:dict[str,list[pysam.AlignedSegment]] = dict()
#     dict_consensus:dict[str,consensus_class.Consensus] = dict()
#     for consensusID in consensusIDs:
#         consensus:consensus_class.Consensus = next(util.yield_consensus_objects(consensusIDs=[consensusID], path_db=path_db,silent=True))
#         log.info(f"loaded consensus object {consensusID} with {len(consensus.reconstructible_reads)} reconstructible reads")
#         consensus_alignments:list[pysam.AlignedSegment] = []
#         for recr in consensus.reconstructible_reads:
#             aln_pysam = recr.alignment.to_pysam()
#             seq = consensus_lib.reconstruct_ReconstructibleSequence(reconstructibleSequence=recr,reference=consensus.consensus_sequence)
#             if aln_pysam.is_reverse:
#                 seq.sequence = str(Seq(seq.sequence).reverse_complement())
#             aln_pysam.query_sequence = seq.sequence
#             if len(seq.sequence) != aln_pysam.infer_query_length():
#                 log.error(f"reconstructed sequence length {len(seq.sequence)} does not match alignment length {aln_pysam.infer_query_length()}")
#                 # print the erroneous recr as string to terminal
#                 print(json.dumps(recr.unstructure()))
#                 raise f"reconstructed sequence length {len(seq.sequence)} does not match alignment length {aln_pysam.infer_query_length()}"

#         #     consensus_alignments.append(aln_pysam)
#         # for each consensus_alignment, reconstruct the sequence and add it to the pysam alignment object
#         dict_alignments[consensusID] = consensus_alignments
#         dict_consensus[consensusID] = consensus

#     # open an output bam file and write all alignmments to it
#     path_bam = str(path_base_output) + ".alignments.bam"
#     # load all alignments

#     # write all pysam alignments to the bam file
#     log.info(f"writing all alignments to {path_bam}")
#     with pysam.AlignmentFile(path_bam,"wb",header=dict_alignments[consensusIDs[0]][0].header) as f:
#         for consensusID in consensusIDs:
#             for alignment in dict_alignments[consensusID]:
#                 f.write(alignment)
#     # with pysam.AlignmentFile(path_bam.replace('bam', 'sam'),"w",header=consensus_alignments[0].header) as f:
#     #     for alignment in consensus_alignments:
#     #         f.write(alignment)
#     # index the bam file
#     log.info(f"indexing {path_bam}")
#     cmd_index = f"samtools index {path_bam}"
#     subprocess.run(shlex.split(cmd_index),check=True)
#     # write the consensus to a fasta file
#     path_fasta = str(path_base_output) + ".consensuses.fasta"
#     log.info(f"writing consensus sequence to path_fasta")
#     with open(path_fasta,"w") as f:
#         for consensus in dict_consensus.values():
#             seq = datatypes.SequenceObject(
#                 name=consensus.ID,
#                 id=consensus.ID,
#                 sequence=consensus.consensus_sequence)
#             sr = parse_SequenceObject_to_BioSeq(seq)
#             SeqIO.write(sr,f,"fasta")
#     # index the fasta file
#     log.info(f"indexing {path_fasta}")
#     cmd_index = f"samtools faidx {path_fasta}"
#     subprocess.run(shlex.split(cmd_index),check=True)
#     return consensusIDs


# def svPrimitives_annotations_for_consensus_sequence(
#         path_db:Path,
#         consensusID:str,
#         reference_name:str):
#     """reads the svPrimitives from the database that have matching consensusID and writes a BED file of consensusName, start, end, infostring (svtype, ref coords)"""
#     svPrimitives:list[SVprimitive_class.SVprimitive] = sv_calling.read_svPrimitives_from_db(
#         database=path_db,consensusIDs=[consensusID],reference_name=reference_name)[consensusID]
#     # for each svPrimitive, calculate the start and end position on the consensus sequence
#     # and fill a list with tuples (start,end,svtype)
#     dict_sv_types = {0:'INS',1:'DEL',3:'BND',4:'BND'}
#     buffer={0:5,1:5,3:30,4:30}
#     intervals:list[tuple[str,int,int,str]] = []
#     for svPrimitive in svPrimitives:
#         infostring = f"SVTYPE={dict_sv_types[svPrimitive.sv_type]};REF={svPrimitive.chr}:{svPrimitive.ref_start}-{svPrimitive.ref_end}"
#         intervals.append((consensusID,svPrimitive.read_start-buffer[svPrimitive.sv_type], svPrimitive.read_end+buffer[svPrimitive.sv_type], infostring))
#     return intervals


# def get_cut_reads_from_consensus(
#         path_db:Path,
#         crID:int,
#         path_base_output:Path,
#         reference_name:str) -> None:
#     """reads the consensus objects with matching crIDs from the database and gets all reconstructible_reads and writes those alignments to a bam file \
# the bam file is written to path_base_output+.crID+.bam and indexed. The consensus is written to path_base_output+crID+.fasta and indexed"""
#     consensusIDs = write_consensus_sequence_and_alignments_to_files(
#         path_db=path_db,
#         crID=crID,
#         path_base_output=path_base_output)
#     log.info(f"found {str(consensusIDs)} consensus sequences. Retrieving annotations for them..")
#     all_annotations = []
#     for consensusID in consensusIDs:
#         all_annotations.append(svPrimitives_annotations_for_consensus_sequence(
#             path_db=path_db,
#             consensusID=str(consensusID),
#             reference_name=reference_name))
#     # all_annotations is a list of lists. flatten it and then write it to a bed file
#     flat_annotations = [item for sublist in all_annotations for item in sublist]
#     with open(str(path_base_output) + ".annotations.bed","w") as f:
#         for annotation in flat_annotations:
#             print(*annotation,sep="\t",file=f)


# def run(args,**kwargs):
#     get_cut_reads_from_consensus(
#         path_db=args.database,
#         crID=args.crID,
#         path_base_output=args.output,
#         reference_name=args.reference_name)


# def get_parser():
#     parser = argparse.ArgumentParser(description="Get the cut reads for a given crID from the database and write them to a bam file.")
#     parser.add_argument("-c","--crID",type=int,required=True,help="The crID to get the cut reads for.")
#     parser.add_argument("-d","--database",type=Path,required=True,help="The path to the database.")
#     parser.add_argument("-o","--output",type=Path,required=True,help="The path to the output directory with base filename.")
#     parser.add_argument("-r","--reference-name",type=str,required=True,help="The name of the genome reference, e.g. hs37d5 or GRCh38.")

#     return parser


# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()

# # %%
