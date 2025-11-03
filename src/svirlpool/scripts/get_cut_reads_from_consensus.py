# import argparse
# from pathlib import Path
# from Bio import SeqIO, SeqRecord
# from Bio.Seq import Seq

# from logzero import logger as log

# from . import datatypes, consensus_containers_to_db, consensus_lib, util


# def parse_SequenceObject_to_BioSeq(sequenceObject:datatypes.SequenceObject) -> SeqRecord:
#     sr = SeqRecord.SeqRecord(
#         Seq(sequenceObject.sequence),
#         id=sequenceObject.name,
#         description="" if not sequenceObject.description else sequenceObject.description)
#     if sequenceObject.qualities and len(sequenceObject.qualities) == len(sequenceObject.sequence):
#         sr.letter_annotations["phred_quality"] = sequenceObject.qualities
#     return sr


# def get_cut_reads_from_consensus(
#         path_db:Path,
#         path_reference_index:Path,
#         threads:int,
#         crID:int,
#         path_fasta:Path,
#         path_bam:Path) -> None:
#     # load all consensusIDs from the database that start with crID
#     log.info(f"loading all consensusIDs from the database {path_db}")
#     consensusIDs = consensus_containers_to_db.load_all_consensusIDs_from_database(path_database=path_db)
#     # filter all consensusIDs that start with crID
#     log.info(f"filtering all consensusIDs that start with {crID}")
#     consensusIDs = [consensusID for consensusID in consensusIDs if consensusID.startswith(str(crID)+'.')]
#     # a consensus object has a consensusID and a list of reconstructible_reads
#     # load all consensus objects from the database
#     # create a dict consensusID:[reconstructed_reads]
#     log.info(f"loading all consensus objects with IDs: {consensusIDs} from the database {path_db}")
#     dict_reconstructed_reads = dict()
#     for consensusID in consensusIDs:
#         consensus:consensus_class.Consensus = next(util.yield_consensus_objects(
#             consensusIDs=[consensusID], path_db=path_db))
#         # reconstruct all reads and save them to dict_reconstructed_reads
#         dict_reconstructed_reads[consensusID] = [
#             parse_SequenceObject_to_BioSeq(consensus_lib.reconstruct_ReconstructibleSequence(
#                 reconstructibleSequence=reconstructible,reference=consensus.consensus_sequence)) \
#                     for reconstructible in consensus.reconstructible_reads]
#         # append the consensusID to each readname
#         for sr in dict_reconstructed_reads[consensusID]:
#             sr.id = f"{sr.id}.{consensusID}"

#     # now write all reads to a fasta file (path_fasta)
#     log.info(f"writing all reads to {path_fasta}")
#     with open(path_fasta,"w") as f:
#         for consensusID in consensusIDs:
#             for sr in dict_reconstructed_reads[consensusID]:
#                 SeqIO.write(sr,f,"fasta")

#     # align all reads to the reference index
#     log.info(f"aligning all reads to the reference index {path_reference_index}")
#     util.align_reads_with_minimap(
#         reference=path_reference_index,
#         bamout=path_bam,
#         reads=path_fasta,
#         tech="map-ont",
#         threads=threads)


# def run(args,**kwargs):
#     get_cut_reads_from_consensus(
#         path_db=args.input,
#         path_reference_index=args.reference_index,
#         threads=args.threads,
#         crID=args.consensusID,
#         path_fasta=args.output_fasta,
#         path_bam=args.output_bam)


# def get_parser():
#     parser = argparse.ArgumentParser(description="Call SVs from the final database and write them to a vcf file.")
#     parser.add_argument("-c","--consensusID", type=int, required=True, help="The consensusID that should be used to extract the reads.")
#     parser.add_argument("-i","--input", type=Path, required=True, help="Path to the database that contains consensus objects.")
#     parser.add_argument("-r","--reference-index", type=Path, required=True, help="Path to the minimap2 index")
#     parser.add_argument("-t","--threads", type=int, default=1, help="Number of threads to use.")
#     parser.add_argument("-of","--output-fasta", type=Path, required=True, help="Path to the output fasta file.")
#     parser.add_argument("-ob","--output-bam", type=Path, required=True, help="Path to the output bam file.")
#     return parser


# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()
