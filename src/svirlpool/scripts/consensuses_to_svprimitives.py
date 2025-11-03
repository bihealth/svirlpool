# # align consensus sequences to a given reference with a provided trf file
# # 1) iterate all consensus containers and their consensus sequences in database and write the sequences to fasta file
# # 2) align the consensus sequences to the reference with minimap2 (map-ont)
# # 3) build intervaltree of the trf annotations
# # 4) load the alignments and sort by query_name (==consensusID)
# # 5)
# # 6) iterate consensus DB and

# #%%
# import argparse
# from pathlib import Path,PosixPath
# import subprocess
# import tempfile
# import json
# import sqlite3
# import cattrs
# import pickle
# import sys
# import itertools

# from logzero import logger as log
# from logzero import logfile
# from intervaltree import IntervalTree
# import pysam
# from tqdm import tqdm
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord, Seq

# from . import consensus, datatypes, util, consensus_align_lib, sv_calling, SVComplexes, consensus_class

# #%%

# def parse_crs_container_results(path: PosixPath):
#     with open(path, "r") as f:
#         for line in f:
#             crs_container_result = json.loads(line)
#             yield cattrs.structure(crs_container_result, consensus_class.CrsContainerResult)

# def align_consensus_sequences_to_reference(
#         path_reference:Path,
#         path_consensus:Path,
#         bamout:Path,
#         threads:int,
#         aln_args:str='',
#         tech:str='asm5',
#         tmp_dir_path:Path|None=None) -> None:
#     with tempfile.TemporaryDirectory(dir=tmp_dir_path,delete=False if tmp_dir_path else True) as tmp_dir:
#         #tmp_bamout = tempfile.NamedTemporaryFile(dir=tmp_dir, suffix=".bam", delete=False if tmp_dir_path else True)
#         log.info("aligning consensus sequences to reference with minimap2")
#         # check if a minimap2 index exists
#         if Path(f"{path_reference}.mmi").exists():
#             log.info("using existing minimap2 index")
#             path_reference = Path(f"{path_reference}.mmi")
#         util.align_reads_with_minimap(
#             reference=path_reference,
#             bamout=bamout,
#             reads=path_consensus,
#             tech=tech,
#             threads=threads,
#             aln_args=aln_args)
#             # bamout=tmp_bamout.name,
#         # refine alignments by aligning the clipped tails of the consensus sequences to the reference
#         # TODO: fully implement
#         # util.add_clipped_tails_to_alignments(
#         #     alignments=tmp_bamout.name,
#         #     reference=path_reference,
#         #     threads=threads,
#         #     tech=tech,
#         #     output=bamout,
#         #     aln_args=aln_args,
#         #     min_clipped_length=300)


# def load_alignments_with_sequences(
#         consensus_to_target_alignments:Path,
#         path_target_reference:Path,
#         tmp_dir_path:Path|None=None) -> dict[str,tuple[list[datatypes.Alignment],list[datatypes.SequenceObject]]]:
#     # load all alignments and parse from pysam.alignedSegment to datatypes.Alignment
#     log.info(f"loading alignments from {consensus_to_target_alignments}...")
#     dict_consensus_to_target_alignments:dict[str,list[datatypes.Alignment]] = {}
#     dict_regions:dict[tuple[str,int],str] = {}
#     with open(consensus_to_target_alignments, "rb") as f:
#         tmp_bed = tempfile.NamedTemporaryFile(dir=tmp_dir_path,delete=False if tmp_dir_path else True, suffix=".bed")
#         for aln in pysam.AlignmentFile(f):
#             if aln.is_mapped: # no mapped consensuses -> empty list
#                 if aln.query_name not in dict_consensus_to_target_alignments:
#                     dict_consensus_to_target_alignments[aln.query_name] = []
#                 dict_consensus_to_target_alignments[aln.query_name].append(
#                         util.parse_alignment(aln,parse_qualities=False))
#                 region = f"{aln.reference_name}:{aln.reference_start}-{aln.reference_end}"
#                 index = len(dict_consensus_to_target_alignments[aln.query_name])-1
#                 dict_regions[(aln.query_name,index)] = region

#     # --- fetch reference sequences --- #
#     # construct bed file of all regions (chr/tstart/tend) of all consensus to ref alignments ( if aligned to the reference )
#     tmp_bed = tempfile.NamedTemporaryFile(dir=tmp_dir_path,delete=False if tmp_dir_path else True, suffix=".bed", prefix="consensus_to_target_alignments.reference_regions.")
#     with open(tmp_bed.name, "w") as f:
#         for alignments in dict_consensus_to_target_alignments.values():
#             for aln in alignments:
#                 pysam_aln:pysam.AlignedSegment = aln.to_pysam()
#                 print(pysam_aln.reference_name, pysam_aln.reference_start, pysam_aln.reference_end, pysam_aln.query_name,sep="\t", file=f)

#     # fetch reference sequences
#     cmd_getfasta = ['bedtools','getfasta','-fi',str(path_target_reference),'-bed',tmp_bed.name]
#     tmp_alt_fasta = tempfile.NamedTemporaryFile(mode='w',delete=False,suffix='.fasta')
#     with open(tmp_alt_fasta.name,'w') as f:
#         log.info(f"extracting ref sequences into {tmp_alt_fasta.name}...")
#         subprocess.check_call(cmd_getfasta,stdout=f)
#     # load fasta file
#     ref_sequences = {record.name:str(record.seq) for record in SeqIO.parse(tmp_alt_fasta.name, "fasta")}
#     # create a dict consensusID:([alignments],[ref_sequences])
#     results = {}
#     for consensusID in dict_consensus_to_target_alignments.keys():
#         alignments = dict_consensus_to_target_alignments[consensusID]
#         ref_seqs = [datatypes.SequenceObject(name=f"{consensusID}.{i}",id=f"{consensusID}.{i}",description=f"consensusID={consensusID},alignmentID={i}",sequence=ref_sequences[dict_regions[(consensusID,i)]]) for i,aln in enumerate(alignments)]
#         results[consensusID] = (alignments,ref_seqs)
#     return results

# # def construct_final_database(
# #         path_database:Path) -> None:
# #     # use sqlite3 to create a database with tables:
# #     # - consensusess [ID,consensus object] (from crsContainerResult)
# #     # - unused_reads [readname,crIDs,SequenceObject] (from crsContainerResult) # crIDs is a list of INTs
# #     conn = sqlite3.connect('file:'+str(path_database)+'?mode=rwc',uri=True)
# #     c = conn.cursor()
# #     c.execute("""CREATE TABLE IF NOT EXISTS consensuses
# #         (id VARCHAR(200) PRIMARY KEY,
# #         consensus TEXT)""")
# #     conn.commit()
# #     c.execute("""CREATE TABLE IF NOT EXISTS unused_reads
# #         (aug_aug_readname VARCHAR(200) PRIMARY KEY,
# #         crID INTEGER,
# #         sequenceObject TEXT)""")
# #     conn.commit()
# #     conn.close()

# # def write_consensus_crsContainerResult_to_final_database(
# #         path_database:Path,
# #         crsContainerResults:list[consensus_class.CrsContainerResult]):
# #     # iterate over crsContainerResults
# #     # and write all consensus objects to the database
# #     # and write all unused reads to the database. append to the aug_name of each read its crID
# #     conn = sqlite3.connect('file:'+str(path_database)+'?mode=rwc',uri=True)
# #     for crsContainerResult in tqdm(crsContainerResults):
# #         c = conn.cursor()
# #         # write all consensus objects to the database
# #         c.executemany("INSERT OR REPLACE INTO consensuses (id,consensus) VALUES (?,?)",
# #             [(consensusID,pickle.dumps(cattrs.unstructure(consensus))) for consensusID,consensus in crsContainerResult.consensus_dicts.items()])
# #         # write all unused reads to the database
# #         conn.commit()
# #         c.executemany("INSERT OR REPLACE INTO unused_reads (aug_aug_readname,crID,sequenceObject) VALUES (?,?,?)",
# #             [[f"{sequenceObject.name}.{sequenceObject.description}",int(sequenceObject.description),pickle.dumps(sequenceObject.unstructure())] for sequenceObject in crsContainerResult.unused_reads])
# #         conn.commit()
# #     conn.close()

# #%%


# def load_alignments(path_alignments:Path,samplename:str|None=None,parse_DNA:bool=True) -> dict[str,list[datatypes.Alignment]]:
#     """produces a dict consensusID:list[Alignment]. If a consensus has no alignments, the list is empty"""
#     # load all alignments and parse from pysam.alignedSegment to datatypes.Alignment
#     log.info(f"loading alignments from {path_alignments}...")
#     dict_consensus_to_alignments:dict[str,list[datatypes.Alignment]] = {}
#     with open(path_alignments, "rb") as f:
#         for aln in pysam.AlignmentFile(f):
#             if aln.is_mapped: # no mapped consensuses -> empty list
#                 if aln.query_name not in dict_consensus_to_alignments:
#                     dict_consensus_to_alignments[aln.query_name] = []
#                 if not parse_DNA:
#                     aln.query_sequence = None
#                 dict_consensus_to_alignments[aln.query_name].append(
#                         util.parse_alignment(aln=aln,samplename=samplename,parse_qualities=False,parse_sequence=False))
#     return dict_consensus_to_alignments

# # test!
# def trf_to_interval_tree(
#         path_trf:Path) -> dict[str,IntervalTree]:
#     """produces a dict chrom:IntervalTree"""
#     # load the trf intervals to an intervaltree and set their index as their names
#     log.info(f"loading trf annotations from {path_trf}...")
#     trf_intervals = {}
#     with open(path_trf, "r") as f:
#         for i,line in enumerate(f):
#             chrom,start,end = line.strip().split()
#             if chrom not in trf_intervals:
#                 trf_intervals[chrom] = IntervalTree()
#             trf_intervals[chrom][int(start):int(end)] = i
#     return trf_intervals

# def add_trf_annotations_to_alignments(
#         alignments:dict[str,list[datatypes.Alignment]],
#         trf_intervals:dict[str,IntervalTree],
#         reference_name:str) -> dict[str,list[datatypes.ConsensusAlignment]]:
#     """Finds the trf annotations for each alignment and produces a dict consensusID:list[tuple[Alignment,(chrom,start,end)]
#     WARNING: This function is destructuve and modifies the alignments in place.

#     Args:
#         alignments (dict[str,list[datatypes.Alignment]]): _description_
#         trf_intervals (dict[str,IntervalTree]): _description_
#         reference_name (str): the name of the reference, e.g. GRCh38 or HS1

#     Returns:
#         dict[str,list[datatypes.ConsensusAlignment]]: dict reference_name:dict consensusID:list[ConsensusAlignment]
#     """
#     dict_results:dict[str,list[datatypes.ConsensusAlignment]] = {} # dict reference_name:dict consensusID:list[ConsensusAlignment]
#     for consensusID,alns in alignments.items():
#         if consensusID not in dict_results:
#             dict_results[consensusID] = []
#         consensusAlignments:dict[str,list] = {reference_name:[]}
#         for aln in alns:
#             assert type(reference_name) == str, f"reference_name is not a str, but {type(reference_name)}, data={reference_name}"
#             consensusAlignment:datatypes.ConsensusAlignment = datatypes.ConsensusAlignment(alignment=aln, trf_intervals=[],reference_name=reference_name)
#             pysam_aln = aln.to_pysam() #type: ignore
#             chrom = aln.reference_name #type: ignore
#             if chrom in trf_intervals:
#                 for interval in trf_intervals[chrom][pysam_aln.reference_start:pysam_aln.reference_end]:
#                     assert type(interval.data) == int, f"interval.data is not an int (repeatID), but {type(interval.data)}, data={interval.data}"
#                     consensusAlignment.trf_intervals.append((interval.begin,interval.end,interval.data))
#             consensusAlignments[reference_name].append(consensusAlignment)
#             aln=None # clear memory
#         dict_results[consensusID] = consensusAlignments[reference_name]
#     return dict_results


# def create_database(output_db:Path,name_reference:str) -> None:
#     # use sqlite3 to create a database with tables: consensusAlignments+'__'+reference_name of form consensusID:list[ConsensusAlignment]
#     # each list[ConsensusAlignment] is unstructured and pickled
#     # the types are consensusID:TEXT, consensusAlignments:BLOB
#     # only create a new file if the file does not exist, else just add the table
#     # '.' in name_reference is illegal in sqlite3
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.','_')
#         log.warning(f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}")
#     conn = sqlite3.connect('file:'+str(output_db)+'?mode=rwc',uri=True)
#     conn.execute(f"""CREATE TABLE IF NOT EXISTS consensusAlignments__{name_reference}
#         (consensusID VARCHAR(100) PRIMARY KEY,
#         consensusAlignments BLOB)""")
#     conn.commit()
#     conn.close()


# # TODO: speed up
# # this is extremely slow
# def write_consensus_alignments_to_database(
#         path_db:Path,
#         dict_alignments:dict[str,list[datatypes.ConsensusAlignment]],
#         name_reference:str,
#         is_initial_reference:bool,
#         cache_size:int=1000) -> None:
#     # iterate over dict_alignments and write all alignments to the database
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.','_')
#         log.warning(f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}")

#     # witht he is_initial_reference annotation, we can check if the alignments in this database are suitable for lift-over operations
#     if is_initial_reference: # add is_initial_reference annotation to the database
#         util.insert_table_metadata(path_database=path_db, table_name=f"consensusAlignments__{name_reference}", tag="is_initial_reference")
#     else:
#         util.insert_table_metadata(path_database=path_db, table_name=f"consensusAlignments__{name_reference}", tag="is_switch_reference")

#     # write all alignments to the database
#     conn = sqlite3.connect('file:'+str(path_db)+'?mode=rwc',uri=True)
#     cache = []
#     query = f"INSERT OR REPLACE INTO consensusAlignments__{name_reference} (consensusID,consensusAlignments) VALUES (?,?)"
#     for consensusID,consensusAlignments in tqdm(dict_alignments.items()):
#         data = [consensusID,pickle.dumps(cattrs.unstructure(consensusAlignments))]
#         cache.append(data)
#         # c.execute(f"INSERT OR REPLACE INTO consensusAlignments__{name_reference} (consensusID,consensusAlignments) VALUES (?,?)", data)
#         if len(cache) >= cache_size:
#             c = conn.cursor()
#             c.executemany(query, cache)
#             conn.commit()
#             cache = []
#             c.close()
#         else:
#             continue
#     if len(cache) > 0:
#         c = conn.cursor()
#         c.executemany(query, cache)
#         conn.commit()
#         c.close()
#     conn.close()

# import typing

# def yield_consensus_alignments_from_database(
#         path_db:Path,
#         name_reference:str) -> typing.Generator[dict[str,list[datatypes.ConsensusAlignment]],None,None]:
#     if "." in name_reference:
#         name_reference = name_reference.replace(".","_")
#         log.warning(f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}")
#     conn = sqlite3.connect('file:'+str(path_db)+'?mode=ro',uri=True)
#     c = conn.cursor()
#     c.execute(f"SELECT consensusID,consensusAlignments FROM consensusAlignments__{name_reference}")
#     for row in c:
#         consensus_alignment = cattrs.structure(pickle.loads(row[1]),list[datatypes.ConsensusAlignment])
#         yield {row[0]:consensus_alignment}
#     c.close()
#     conn.close()

# def load_consensus_alignments_from_database(
#         path_db:Path,
#         name_reference:str) -> dict[str,list[datatypes.ConsensusAlignment]]:
#     # iterate over dict_alignments and write all alignments to the database
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.','_')
#         log.warning(f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}")
#     conn = sqlite3.connect('file:'+str(path_db)+'?mode=ro',uri=True)
#     c = conn.cursor()
#     c.execute(f"SELECT consensusID,consensusAlignments FROM consensusAlignments__{name_reference}")
#     results = {}
#     for row in c:
#         results[row[0]] = cattrs.structure(pickle.loads(row[1]),list[datatypes.ConsensusAlignment])
#     c.close()
#     conn.close()
#     return results

# def load_specific_consensus_alignments_from_database(
#         path_db:Path,
#         name_reference:str,
#         consensusIDs:set[str]) -> dict[str,list[datatypes.ConsensusAlignment]]:
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.','_')
#         log.warning(f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}")
#     conn = sqlite3.connect('file:'+str(path_db)+'?mode=ro',uri=True)
#     c = conn.cursor()
#     c.execute(f"SELECT consensusID,consensusAlignments FROM consensusAlignments__{name_reference} WHERE consensusID IN ({','.join(['?']*len(consensusIDs))})",tuple(consensusIDs))
#     results = {}
#     for row in c:
#         results[row[0]] = cattrs.structure(pickle.loads(row[1]),list[datatypes.ConsensusAlignment])
#     c.close()
#     conn.close()
#     return results

# # faster than querying the database fro Consensus objects (in load_consensus_sequences_from_db)
# # def load_consensus_sequences_from_alignments(
# #         path_alignments:Path) -> dict[str,str]:
# #     result:dict[str,str] = dict()
# #     for aln in pysam.AlignmentFile(path_alignments):
# #         # if aln has hard clips, continue
# #         if aln.is_unmapped:
# #             result[aln.query_name] = aln.query_sequence
# #         elif aln.cigartuples[0][0] != 5 and aln.cigartuples[-1][0] != 5:
# #             result[aln.query_name] = aln.query_sequence
# #         else:
# #             continue
# #     return result


# def simplified_consensus_sequence(consensus:consensus_class.Consensus) -> bytes:
#     s,e = consensus.consensus_padding.consensus_interval_on_sequence_with_padding
#     core = consensus.consensus_padding.sequence[s:e]
#     left = "X" * consensus.consensus_padding.padding_size_left
#     right = "X" * consensus.consensus_padding.padding_size_right
#     return pickle.dumps(left + core + right)


# def align_consensuses_to_reference(
#         path_db:Path,
#         path_reference:Path,
#         consensus_to_reference_alignments:Path,
#         path_consensus_fasta:Path|None=None,
#         threads:int=8,
#         tmp_dir_path:Path|None=None) -> tuple[dict[str,bytes], dict[str,tuple[int,int]]]:

#     core_sequences:dict[str,bytes] = {} # save key = consensusID, value = sequence (compressed with pickle and simplified padding parts)
#     core_intervals:dict[str,tuple[int,int]] = {} # save key = consensusID, value = (start,end) of the core sequence
#     with tempfile.TemporaryDirectory(dir=tmp_dir_path, delete=False if tmp_dir_path else True) as tmp_dir:
#         if path_consensus_fasta is None:
#             tmp_consensus_fasta = tempfile.NamedTemporaryFile(dir=tmp_dir,delete=False if tmp_dir_path else True, suffix=".fasta",prefix="consensus_sequences.")
#             path_consensus_fasta = Path(tmp_consensus_fasta.name)

#         with open(path_consensus_fasta, "w") as f:
#             for consensus_object in util.yield_consensus_objects(path_db=path_db):
#                 core_sequences[consensus_object.ID] = pickle.dumps(consensus_object.consensus_sequence) # compress for smaller memory footprint
#                 core_intervals[consensus_object.ID] = consensus_object.consensus_padding.consensus_interval_on_sequence_with_padding
#                 seq = consensus_object.consensus_padding.sequence # write full padded sequence to the fasta file
#                 id = consensus_object.ID
#                 # description is for debugging purposes
#                 description = f"padding_size_left={str(consensus_object.consensus_padding.padding_size_left)},padding_size_right={str(consensus_object.consensus_padding.padding_size_right)},consensus_interval_on_sequence_with_padding={consensus_object.consensus_padding.consensus_interval_on_sequence_with_padding}"
#                 seqRec = SeqRecord(
#                     seq=Seq(seq),
#                     id=id,
#                     name=id,
#                     description=description)
#                 SeqIO.write(seqRec, f, "fasta")

#         log.info(f"aligning consensus sequences to the reference {path_reference}...")
#         #TODO(feature): enable secondary alignments to find alternative spots for consensus sequences on the final references when SV calling.
#         aln_args = " --secondary=no "
#         tmp_bam_out_unfiltered = tempfile.NamedTemporaryFile(dir=tmp_dir,delete=False if tmp_dir_path else True, suffix=".bam",prefix="consensus_to_reference_unfiltered.")
#         align_consensus_sequences_to_reference(
#             path_consensus=path_consensus_fasta,
#             path_reference=path_reference,
#             bamout=tmp_bam_out_unfiltered.name,
#             aln_args=aln_args,
#             threads=threads,
#             tech="map-ont")
#     # filter the consensus alignments for local alignments that are not covering the core intervals
#     log.info(f"filtering the consensus alignments for local alignments that are not covering the core intervals...")
#     with pysam.AlignmentFile(tmp_bam_out_unfiltered, 'rb') as f:
#         header = f.header
#         log.info(f"filtering the consensus alignments for local alignments that are not covering the core intervals...")
#         with pysam.AlignmentFile(consensus_to_reference_alignments, 'wb', header=header) as g:
#             for aln in f:
#                 if aln.is_mapped:
#                     # check if the alignment covers the core interval
#                     traced_back_ref_start, traced_back_ref_end = util.get_interval_on_ref_in_region(a=aln,start=core_intervals[aln.query_name][0],end=core_intervals[aln.query_name][1])
#                     if traced_back_ref_start != traced_back_ref_end:
#                         g.write(aln)
#     # index the bam file
#     cmd_index = ['samtools','index',consensus_to_reference_alignments]
#     subprocess.check_call(cmd_index)
#     return core_sequences, core_intervals


# def consensus_alignments_to_svPrimitives(
#         samplename:str,
#         path_db:Path,
#         path_reference:Path,
#         name_reference:str,
#         output_svprimitives:Path,
#         output_svcomplexes:Path,
#         path_trf:Path,
#         output_consensus_alignments_db:Path,
#         threads:int,
#         path_consensus_fasta:Path|None=None,
#         consensus_to_reference_alignments:Path|None=None,
#         chunksize:int=100,
#         tmp_dir_path:Path|None=None,
#         sqlite_timeout:float|None=None,
#         dont_merge_horizontally:bool=False,
#         is_initial_reference:bool=False) -> None:
#     # write all consensus sequences to a fasta file and align to the target reference
#     # write all alignments to a file
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.','_')
#         log.warning(f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}")
#     if dont_merge_horizontally:
#         log.warning("dont_merge_horizontally is set to True. This will not merge SVs within the same consensus alignment. This is not recommended for general use cases.")

#     with tempfile.TemporaryDirectory(dir=tmp_dir_path, delete=False if tmp_dir_path else True) as tmp_dir:
#         if consensus_to_reference_alignments is None:
#             consensus_to_reference_alignments = Path(tmp_dir) / "consensus_to_reference.bam"

#         sequences_core, intervals_core = align_consensuses_to_reference(
#             path_db=path_db,
#             path_reference=path_reference,
#             consensus_to_reference_alignments=consensus_to_reference_alignments,
#             path_consensus_fasta=path_consensus_fasta,
#             threads=threads,
#             tmp_dir_path=tmp_dir)

#         log.info(f"loading alignments from {consensus_to_reference_alignments}...")
#         consensus_alignments:dict[str,list[datatypes.Alignment]] = load_alignments(path_alignments=consensus_to_reference_alignments,parse_DNA=False)
#         n_alignments = sum([len(alignments) for alignments in consensus_alignments.values()])
#         log.info(f"{n_alignments} consensus alignments generated.")
#         # TODO: check correct number of alignments

#         log.info(f"adding tandem repeat annotations to the consensus-to-reference alignments...")
#         dict_alignments:dict[str,list[datatypes.ConsensusAlignment]] = add_trf_annotations_to_alignments(
#             alignments=consensus_alignments,
#             trf_intervals=trf_to_interval_tree(path_trf),
#             reference_name=name_reference)
#         consensus_alignments=dict() # has no real use anymore, free memory
#         # add samplename to consensusAlignments.alignment members
#         for consensusID,consensusAlignments in dict_alignments.items():
#             for consensusAlignment in consensusAlignments:
#                 consensusAlignment.alignment.samplename = samplename

#         unaligned_consensusIDs = set(sequences_core.keys()) - set(dict_alignments.keys())
#         if len(unaligned_consensusIDs) > 0:
#             if len(unaligned_consensusIDs) == len(dict_alignments):
#                 raise ValueError(f"all consensus sequences are unaligned to the reference. Either the path to the alignments or the fasta file is corrupted or you need to allocate more memory.")
#             log.info(f"{len(unaligned_consensusIDs)} consensus sequences were not aligned to the reference. They are: {unaligned_consensusIDs}")

#         for consensusID,consensusAlignments in tqdm(dict_alignments.items()):
#             for consensusAlignment in consensusAlignments:
#                 # add proto_svs (MergedSVSignal)s to all ConsensusAlignments in dict_alignments
#                 # contains the alt sequences
#                 sequence = pickle.loads(sequences_core[consensusID]) # is just forward or backward.
#                 merged_svs:list[datatypes.MergedSVSignal] = consensus_align_lib.parse_sv_signals_from_consensus(
#                     samplename=samplename,
#                     consensusAlignment=consensusAlignment,
#                     consensus_sequence=sequence,
#                     interval_core=intervals_core[consensusID])
#                 consensusAlignment.proto_svs=merged_svs

#         # then add the reference sequences to all MergedSVSignals in proto_svs in all ConsensusAlignments in dict_alignments
#         dict_alignments:dict[str,list[datatypes.ConsensusAlignment]] = consensus_align_lib.add_ref_sequences_to_dict_alignments(
#                 dict_alignments=dict_alignments,
#                 path_reference=path_reference)
#         log.info(f"dict_alignments object has size {sys.getsizeof(dict_alignments)} bytes")

#         # ========================== START DEBUG ========================== #
#         # pickle dict_alignments
#         # debug_path = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/merge_svs_in_dict_alignments_input.1.pkl")
#         # with open(debug_path, 'wt') as f:
#         #     pickle.dumps(dict_alignments, f)
#         # ========================== END DEBUG ========================== #

#         # then truly merge the merged sv signals
#         #kmer_letter_dict = util.compute_letter_dict("ACGTN")
#         log.info(f"merging SV signals in all consensus alignments...")
#         # merge SVs within the same consensus alignment (horizontal merging)
#         if not dont_merge_horizontally:
#             consensus_align_lib.merge_svs_in_dict_alignments(
#                 dict_alignments=dict_alignments,
#                 merge_large_dels_with_small_gaps=(50,500),
#                 max_low_complexity_ignored_size=15)

#         # generate SVcomplexes from the MergedSVSignals in the ConsensusAlignments


#         # create the output database
#         create_database(output_db=output_consensus_alignments_db,name_reference=name_reference)
#         # write the annotated SVs to output_db
#         log.info(f"writing consensus alignments to {output_consensus_alignments_db}...")
#         write_consensus_alignments_to_database(
#             path_db=output_consensus_alignments_db,
#             dict_alignments=dict_alignments,
#             name_reference=name_reference,
#             is_initial_reference=is_initial_reference)

#         # =====================================================================
#         #  generate all SV primitives
#         # sv primitives have their alignments as members
#         svPrimitives:list[SVprimitive_class.SVprimitive] = sv_calling.generate_svPrimitives(
#             samplename=samplename,
#             name_reference=name_reference,
#             consensus_alignments=dict_alignments,
#             path_consensuses=path_db,
#             threads=threads,
#             chunksize=chunksize,
#             DB_timeout=sqlite_timeout)

#         # =====================================================================
#         #  generate all SV complexes
#         # sv primitives at this point are fully garnished: their adjacencies are set, they have ref and alt sequences
#         # consensus alignments are present here
#         svComplexes:list[datatypes.SVcomplex] = SVComplexes.SVprimitives_to_SVcomplexes(svPrimitives=svPrimitives)
#         svPrimitives = SVComplexes.subtract_SVcomplexes_from_SVprimitives(
#             svPrimitives=svPrimitives,
#             svComplexes=svComplexes)

#         SVComplexes.add_ref_sequences_to_svcomplexes_inplace(SVcomplexes=svComplexes, path_reference=path_reference)
#         SVComplexes.add_alt_sequences_to_svcomplexes_inplace(SVcomplexes=svComplexes, intervals_core=intervals_core, sequences_core=sequences_core)

#         # create the output database for SVcomplexes
#         log.info(f"creating output database {output_svcomplexes} table..")
#         SVComplexes.create_svComplexes_db(
#             database=output_svcomplexes,
#             name_reference=name_reference)

#         # then write SVcomplexes to DB
#         SVComplexes.write_svComplexes_to_db(
#             database=output_svcomplexes,
#             data=svComplexes,
#             name_reference=name_reference)


#         # =====================================================================
#         #  write all SV primitives to database

#         # create output database
#         log.info(f"creating output database {output_svprimitives} table..")
#         sv_calling.create_svPrimitves_db(
#             database=output_svprimitives,
#             name_reference=name_reference)

#         sv_calling.write_svPrimitives_to_db(
#             database=output_svprimitives,
#             data=svPrimitives,
#             name_reference=name_reference)

#         # log.info(f"all_svPrimitives object has size {sys.getsizeof(all_svPrimitives)} bytes")

#         # # if allow_list:
#         # #     log.info(f"filtering svPrimitives by allow list {str(allow_list)}...")
#         # #     all_svPrimitives,filtered_svPrimitives = sv_calling.filter_svPrimitives_by_allow_list(allow_list=allow_list,svPrimitives=all_svPrimitives)
#         # #     log.info(f"filtered {len(filtered_svPrimitives)} svPrimitives by allow list with remaining {len(all_svPrimitives)} svPrimitives.")

#         # #log.info(f"merging {len(all_svPrimitives)} svPrimitives...")
#         # #all_svPrimitives:list[SVprimitive_class.SVprimitive] = sv_calling.merge_svPrimitives_with_mp(svPrimitives=all_svPrimitives,threads=threads)

#         # log.debug(f"checking the results...")
#         # # check the results
#         # for svPrimitive in all_svPrimitives:
#         #     # if sv_type is 0, then alt_sequence must have len > 0
#         #     # if sv_type is 1, then ref_sequence must have len > 0
#         #     if svPrimitive.sv_type == 0 and len(svPrimitive.get_alt_sequence()) == 0:
#         #         raise ValueError(f"svPrimitive {svPrimitive} has sv_type {str(svPrimitive.sv_type)} but no alt_sequence")
#         #     if svPrimitive.sv_type == 1 and len(svPrimitive.get_ref_sequence()) == 0:
#         #         raise ValueError(f"svPrimitive {svPrimitive} has sv_type {str(svPrimitive.sv_type)} but no ref_sequence")
#         #     if svPrimitive.sv_type >= 3 and len(svPrimitive.get_alt_sequence()) == 0:
#         #         raise ValueError(f"svPrimitive {svPrimitive} has sv_type {str(svPrimitive.sv_type)} but no alt_sequence")

#         # log.info(f"writing svPrimitives to {output_db} for reference {name_reference}...")
#         # # write svPrimitves to input database
#         # sv_calling.write_svPrimitives_to_db(database=output_db,data=all_svPrimitives, name_reference=name_reference)


#     #%%

# # %%

# def run(args,**kwargs):
#     consensus_alignments_to_svPrimitives(
#         samplename=args.samplename,
#         path_db=args.input,
#         path_reference=args.reference,
#         name_reference=args.name_reference,
#         output_svprimitives=args.svprimitives,
#         output_svcomplexes=args.svcomplexes,
#         path_trf=args.repeats,
#         output_consensus_alignments_db=args.output,
#         threads=args.threads,
#         path_consensus_fasta=args.consensus_fasta,
#         consensus_to_reference_alignments=args.consensus_alignments,
#         chunksize=args.chunksize,
#         sqlite_timeout=args.sqlite_timeout,
#         tmp_dir_path=args.tmp_dir_path,
#         dont_merge_horizontally=args.dont_merge_horizontally,
#         is_initial_reference=args.is_initial_reference)


# # superseeded by __main__:main()
# # def get_parser():
# #     parser = argparse.ArgumentParser(description="Aligns consensus sequences to a reference and create consensus alignment objects that are annotated with ref and alt sequences and are saved to an output database.")
# #     parser.add_argument("-i","--input", type=Path, required=True, help="Path to the input database with the consensus objects.")
# #     parser.add_argument("-r","--reference", type=Path, required=True, help="Path to the reference genome fasta file.")
# #     parser.add_argument("-n","--name-reference", type=str, required=True, help="Name of the reference genome, e.g. hs37d5.")
# #     parser.add_argument("-o","--output", type=Path, required=True, help="Path to the output database that will contain the ConsensusAlignment objects. Can be the same as input, where a table will be added with the name consensusAlignment__[reference_name].")
# #     parser.add_argument("--repeats", type=Path, required=True, help="Path to the tandem repeat annotations file for the chosen reference.")
# #     parser.add_argument("--consensus-fasta", type=Path, required=False, default=None, help="Path to the consensus sequences in fasta format. If not provided, the consensus sequences will be written to a temporary file.")
# #     parser.add_argument("--consensus-alignments", type=Path, required=False, default=None, help="Path to the consensus to reference alignments in bam format. If not provided, the alignments will be written to a temporary file.")
# #     parser.add_argument("-t","--threads", type=int, default=8, help="Number of threads to use.")
# #     #parser.add_argument("--aln-args", type=str, required=False, default='', help="Additional arguments for minimap2.")
# #     parser.add_argument("--tmp-dir-path", type=Path, required=False, default=None, help="Path to the temporary directory.")
# #     parser.add_argument("--logfile", type=Path, required=False, default=None, help="Path to the logfile.")
# #     return parser


# # def main():
# #     parser = get_parser()
# #     args = parser.parse_args()
# #     if args.logfile:
# #         logfile(str(args.logfile))
# #     run(args)
# #     return

# # if __name__ == "__main__":
# #     main()
