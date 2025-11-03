# #%%
# import argparse
# from pathlib import Path
# from collections import defaultdict
# import typing
# import sqlite3
# import time
# from multiprocessing import Pool
# from tqdm import tqdm

# import pysam
# from . import datatypes,datastructures,util

# from logzero import logger as log

# #%%
# # =============================================================================
# #  helper functions
# def overlaps(a:pysam.AlignedSegment, refID:int, start:int, end:int) -> bool:
#     return a.reference_id == refID \
#         and max(0,min(a.reference_end,end) - max(a.reference_start,start)) > 0

# def consecutive_regions(
#         cra:datatypes.CandidateRegion,
#         crb:datatypes.CandidateRegion,
#         crs_fetch_merge_distance:int) -> bool:
#     """returns True if crb is consecutive to cra, False otherwise"""
#     return cra.referenceID == crb.referenceID \
#         and min(cra.referenceEnd,crb.referenceEnd) \
#         - max(cra.referenceStart,crb.referenceStart) \
#             > -crs_fetch_merge_distance

# # step: set of full_connecting_reads
# def read_connections(path_connections:Path,threshold:int) -> typing.Tuple[dict,datastructures.UnionFind]:
#     """Creates a dict of sets of full connecting reads and a UnionFind object from a connections file."""
#     connections = util.load_crs_connections(input=path_connections)
#     if len(connections) == 0:
#         return defaultdict(set), datastructures.UnionFind([]), set()
#     a,b = list(zip(*list(connections.keys())))
#     UF = datastructures.UnionFind([*a,*b])
#     # connect all connections with a number of connecting reads >= threshold
#     for (crA,crB),(bnds,dels) in connections.items():
#         if len(bnds) + len(dels) >= threshold:
#             UF.union_by_name(crA,crB)
#     # create sets of shared reads for representative crs
#     full_reads = defaultdict(set)
#     for (crA,crB),(bnds,dels) in connections.items():
#         for bnd in bnds:
#             full_reads[UF.name_of_id(UF.find_by_name(crA))].add(bnd)
#             full_reads[UF.name_of_id(UF.find_by_name(crB))].add(bnd)
#     connected_crs = set(UF._name_to_id.keys())
#     return full_reads, UF, connected_crs

# # generate a dict of crIDs that map either to themselves or to their representative
# # from a given UnionFind object
# def representative_crs(UF:datastructures.UnionFind,crIDs:set) -> dict:
#     dict_representatives = dict()
#     crIDs_in_UF = set(UF._name_to_id.keys())
#     for crID in crIDs:
#         if crID in crIDs_in_UF:
#             dict_representatives[crID] = UF.name_of_id(UF.find_by_name(crID))
#         else:
#             dict_representatives[crID] = crID
#     return dict_representatives

# # =============================================================================
# #  functions
# # process_crs_on_bams fills two queues:
# #  1. queue_reads: (readname,sampleID,sequence,qualities)
# #  2. queue_intervals: (readname,start,end,refID,refStart,refEnd,sampleID)

# def write_reads_to_reads_db(
#         path_reads_database:Path,
#         buffer_reads:list,
#         timeout:float=300.0) -> None:
#     log.debug(f"writing {len(buffer_reads)} reads to database..")
#     """Writes reads to reads database."""
#     con_reads_db = sqlite3.connect(path_reads_database, timeout=timeout)
#     cur_reads_db = con_reads_db.cursor()
#     # write to db if readname does not exist
#     cur_reads_db.executemany(
#         "INSERT OR IGNORE INTO reads VALUES (?,?,?,?)",
#         buffer_reads
#     )
#     con_reads_db.commit()
#     con_reads_db.close()

# def write_intervals_to_intervals_db(
#         path_intervals_database:Path,
#         buffer_intervals:list,
#         timeout:float) -> None:
#     log.debug(f"writing {len(buffer_intervals)} intervals to database..")
#     """Writes intervals to intervals database."""
#     con_intervals_db = sqlite3.connect(path_intervals_database, timeout=timeout)
#     cur_intervals_db = con_intervals_db.cursor()
#     cur_intervals_db.executemany(
#         "INSERT INTO intervals VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
#         buffer_intervals)
#     con_intervals_db.commit()
#     con_intervals_db.close()

# # this function should be run in parallel
# def process_crs_on_bams(
#         dict_sampleID_bam:dict,
#         crs:list,
#         dict_representatives:dict,
#         path_reads_database:Path,
#         path_intervals_database:Path,
#         max_chunk_size:int,
#         tolerated_clipped_tails_length:int,
#         timeout:float) -> None:
#     # save all reads and intervals in memory
#     # then open connection to both databases
#     # and write to databases
#     #log.debug(f"I work on crs: {str([cr.crID for cr in crs])}")
#     region = (crs[0].chr,crs[0].referenceStart,crs[-1].referenceEnd)
#     used_reads = set([f"{row[8]}.{row[9]}" for cr in crs for row in cr.sv_signals])
#     reads_in_db = set()
#     buffer_reads = []
#     buffer_intervals = []
#     index_buffer_intervals = {}
#     missed_reads = dict()
#     n_index = 0
#     for sampleID,bampath in dict_sampleID_bam.items():
#         missed_reads[sampleID] = set()
#         bamfile = pysam.AlignmentFile(bampath,'rb')
#         # fetch all reads and keep in memory
#         #log.debug(f"fetching reads in region {region}")
#         # filter reads for all readnames in all crs
#         alignments:typing.List[pysam.AlignedSegment] = [aln for aln in bamfile.fetch(*region) if aln.query_name+'.'+str(sampleID) in used_reads]
#         for alignment in alignments:
#             if alignment.cigartuples[0][0] == 5 or alignment.cigartuples[-1][0] == 5:
#                 raise ValueError(f"Read {alignment.query_name} has hard clipped ends. This is not allowed. Use the -Y flag when aligning with minimap2.")
#         # TODO: implement a solution
#         # Problem: hard clipped queries don't bring the full read sequence with them. If relevant subsequences are to be extracted from
#         # the query sequence, then some sequence can be missing, since the hard clipped ends are not part of the query sequence.
#         # to avoid this, so far, minimap2 is provided the -Y flag to only produce soft clipped tails.
#         # But we need another solution, where all alignments of a read are searched, until the full read sequence is found.
#         # this full read sequence can then be used to extract the relevant subsequences.
#         # How to find all alignments of a read?
#         # iterate all alignments and build a dict that has all alignments and their chr:start-end if the alignment has no hard clipped ends
#         # then this dict can be queried for the alignment region of the read to get the alignment from the bam file and fetch its DNA sequence

#         for cr in crs:
#             reads_in_cr = set([f"{row[8]}.{row[9]}" for row in cr.sv_signals])
#             #log.debug(f"working on cr {cr.crID} in region {region} with {len(reads_in_cr)} reads")
#             reads_in_region = set()
#             for aln in alignments:
#                 qname = aln.query_name+'.'+str(sampleID)
#                 if (not qname in reads_in_cr) or aln.is_secondary or aln.is_duplicate or not overlaps(
#                             a=aln,
#                             refID=cr.referenceID,
#                             start=cr.referenceStart,
#                             end=cr.referenceEnd):
#                     #log.debug(f"skipping read {qname} in region {region}")
#                     continue
#                 # must be a number. if not, the read is not in the cr
#                 interval = util.get_interval_on_read_in_region(
#                     a=aln,
#                     start =cr.referenceStart,
#                     end   =cr.referenceEnd)
#                 # check types of interval.
#                 if type(interval[0]) != int or type(interval[1]) != int:
#                     raise ValueError(f"interval[0] and interval[1] must be of type int. Got {type(interval[0])} and {type(interval[1])} instead.")
#                 if qname in reads_in_region:
#                     # extend existing interval
#                     # find read with same name in buffer_intervals
#                     # TODO: adjust to full read length coordinates (?)
#                     other = buffer_intervals[index_buffer_intervals[qname]]
#                     other[3] = min(interval[0],other[3])
#                     other[4] = max(interval[1],other[4])
#                     start_aln_on_read,end_aln_on_read = util.query_start_end_on_read(aln)
#                     other[11] = min(start_aln_on_read,other[11])
#                     other[12] = max(end_aln_on_read,other[12])
#                     other[13] = min(max(aln.reference_start,cr.referenceStart),other[13])
#                     other[14] = max(min(aln.reference_end,cr.referenceEnd),other[14])
#                     #log.debug(f"extended interval for read {qname} in region {region}")
#                     continue
#                 else:
#                     # add existing interval and read to buffers
#                     reads_in_region.add(qname)
#                     # if the alignment does not have hard clipped ends, it is a full read
#                     # and can be written to the database. If the orientation of the alignment is is_reverse,
#                     # then the reverse complement of the sequence and the reverse qualities must be written to the database
#                     # DEBUG
#                     # if qname in reads_in_db:
#                     #     log.debug(f"skipping read {qname}, because it is already in the database")
#                     # if aln.cigartuples is None:
#                     #     log.debug(f"skipping read {qname}, because it has no cigar tuples")
#                     # if aln.cigartuples[0][0] == 5 or aln.cigartuples[-1][0] == 5:
#                     #     log.debug(f"skipping read {qname}, because it has hard clipped ends")
#                     # DEBUG END
#                     if qname not in reads_in_db and \
#                             aln.cigartuples[0][0] != 5 and aln.cigartuples[-1][0] != 5:
#                         reads_in_db.add(qname)
#                         if aln.is_reverse:
#                             seq = ''.join(util.reverse_complement(list(aln.query_sequence)))
#                             qual = aln.query_qualities[::-1]
#                             #log.debug(f"add to buffer reads: {qname},{sampleID}..")
#                             buffer_reads.append(
#                                 (qname,
#                                 sampleID,
#                                 seq,
#                                 qual))
#                             #log.debug(f"added read {qname} to buffer_reads (reverse alignment)")
#                         else:
#                             #log.debug(f"add to buffer reads: {qname},{sampleID}..")
#                             buffer_reads.append(
#                                 (qname,
#                                 sampleID,
#                                 aln.query_sequence,
#                                 aln.query_qualities))
#                             #log.debug(f"added read {qname} to buffer_reads")
#                     start_aln_on_read,end_aln_on_read = util.query_start_end_on_read(aln)
#                     if aln.is_reverse:
#                         full_read_length = end_aln_on_read + aln.cigartuples[0][1] if aln.cigartuples[0][0] in (4,5) else 0
#                     else:
#                         full_read_length = end_aln_on_read + aln.cigartuples[-1][1] if aln.cigartuples[-1][0] in (4,5) else 0
#                     #log.debug(f"add to buffer_intervals: {qname},{cr.crID}..")
#                     buffer_intervals.append(
#                         [qname,
#                         cr.crID,
#                         dict_representatives[cr.crID],
#                         int(interval[0]),
#                         int(interval[1]),
#                         cr.referenceID,#''.join(map(lambda x: chr( x+33 ), read.query_qualities))))
#                         cr.referenceStart,
#                         cr.referenceEnd,
#                         sampleID,
#                         not aln.is_reverse,
#                         int(full_read_length),
#                         start_aln_on_read,
#                         end_aln_on_read,
#                         max(aln.reference_start,cr.referenceStart),
#                         min(aln.reference_end,cr.referenceEnd)])
#                     #log.debug(f"added interval for read {qname} to buffer_intervals")
#                     index_buffer_intervals[qname] = n_index
#                     n_index += 1
#             # if buffer exceeds max_chunk_size, write to database
#             if len(buffer_reads) >= max_chunk_size:
#                 write_reads_to_reads_db(
#                     path_reads_database=path_reads_database,
#                     buffer_reads=buffer_reads,
#                     timeout=timeout)
#                 #log.debug(f"written {len(buffer_reads)} reads to database")
#                 buffer_reads = []
#             if len(buffer_intervals) >= max_chunk_size*100:
#                 write_intervals_to_intervals_db(
#                     path_intervals_database=path_intervals_database,
#                     buffer_intervals=buffer_intervals,
#                     timeout=timeout)
#                 #log.debug(f"written {len(buffer_intervals)} intervals to database")
#                 buffer_intervals = []
#                 index_buffer_intervals = {}
#                 n_index = 0
#         # write remaining reads and intervals to database
#         if len(buffer_reads) > 0:
#             write_reads_to_reads_db(
#                 path_reads_database=path_reads_database,
#                 buffer_reads=buffer_reads,
#                 timeout=timeout)
#             #log.debug(f"written {len(buffer_reads)} reads to database")
#             buffer_reads = []
#         if len(buffer_intervals) > 0:
#             write_intervals_to_intervals_db(
#                 path_intervals_database=path_intervals_database,
#                 buffer_intervals=buffer_intervals,
#                 timeout=timeout)
#             #log.debug(f"written {len(buffer_intervals)} intervals to database")
#             buffer_intervals = []
#             index_buffer_intervals = {}
#             n_index = 0
#         bamfile.close()
#         # report any reads that are missing from the db
#         missing_local_reads = reads_in_region - reads_in_db
#         if len(missing_local_reads) > 0:
#             missed_reads[sampleID].update(missing_local_reads)
#             #log.debug(f"missed {len(missing_local_reads)} reads in region {region}")
#     # log.info(f"Finished working on crs: {str([cr.crID for cr in crs])}")
#     log.debug(f"Finished working on region {region} with crs={str([cr.crID for cr in crs])}")
#     return missed_reads
#     # open conenction to reads database and intervals database

# def process_func(kwargs) -> dict:
#     process_crs_on_bams(**kwargs)

# # =============================================================================
# #  DEBUG INPUT
# # first, create a set of all used readnames by iterating all crs and saving the readnames to a set (in memory)
# # path_crs = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.crs")
# # path_connections = Path('/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.connections')
# # path_reads_db = Path('/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/reads_database.db')
# # path_read_intervals_db = Path('/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/read_intervals_database.db')
# # path_sampledicts = Path("/home/vinzenz/development/LRSV-detection/development/test/HG002HG003/HG002HG003.sampledicts.json")

# # max_threads = 4 # each subprocess gets a set of adjacent crs and reads the intervals from the bam files
# # threshold = 1
# # crs_fetch_merge_distance = 50_000
# # tolerated_clipped_tails_length = 100
# # chunk_size = 1000 # number of reads to be written to database at once

# def create_reads_db(path_reads_db:Path,timeout:float):
#     con_reads = sqlite3.connect(path_reads_db)
#     cur_reads = con_reads.cursor()
#     cur_reads.execute("""CREATE TABLE IF NOT EXISTS reads
#         (readname VARCHAR(200) PRIMARY KEY,
#         sampleID INTEGER,
#         sequence TEXT,
#         qualities TEXT)""")
#     con_reads.execute(f'pragma busy_timeout={str(int(timeout*1000))}')
#     con_reads.commit()
#     con_reads.close()

# def create_intervals_db(path_intervals_db:Path,timeout:float):
#     con_intervals = sqlite3.connect(path_intervals_db)
#     cur_intervals = con_intervals.cursor()
#     cur_intervals.execute("""CREATE TABLE IF NOT EXISTS intervals
#         (readname VARCHAR(200),
#         crID INTEGER,
#         mcrID INTEGER,
#         start INTEGER,
#         end INTEGER,
#         refID INTEGER,
#         refStart INTEGER,
#         refEnd INTEGER,
#         sampleID INTEGER,
#         forward BOOLEAN,
#         qlen INTEGER,
#         qstart INTEGER,
#         qend INTEGER,
#         rstart INTEGER,
#         rend INTEGER)""")
#     con_intervals.execute(f'pragma busy_timeout={str(int(timeout*1000))}')
#     con_intervals.commit()
#     con_intervals.close()

# def write_missed_reads_to_db(
#         dict_sampleID_bam:dict,
#         path_reads_db:Path,
#         missed_reads:dict,
#         timeout:float,
#         chunk_size:int) -> set:
#     totally_missed_reads = set()
#     for sampleID,readnames in missed_reads.items():
#         # fill this buffer with reads that are not in the database until there are chunk_size reads
#         # then write all reads to the database
#         saved_reads = set()
#         buffer = []
#         # open the bam file of sampleID and iterate the alignments there.
#         for aln in tqdm(pysam.AlignmentFile(dict_sampleID_bam[sampleID],'rb')):
#             if aln.query_name in readnames and not aln.query_name in saved_reads:
#                 if aln.cigartuples[0][0] != 5 and aln.cigartuples[-1][0] != 5:
#                     saved_reads.add(aln.query_name)
#                     if aln.is_reverse:
#                         seq = ''.join(util.reverse_complement(list(aln.query_sequence)))
#                         qual = aln.query_qualities[::-1]
#                         buffer.append(
#                             (aln.query_name,
#                             sampleID,
#                             seq,
#                             qual))
#                     else:
#                         buffer.append(
#                             (aln.query_name,
#                             sampleID,
#                             aln.query_sequence,
#                             aln.query_qualities))
#             if len(buffer) >= chunk_size:
#                 write_reads_to_reads_db(
#                     path_reads_database=path_reads_db,
#                     buffer_reads=buffer,
#                     timeout=timeout)
#                 buffer = []
#         if len(buffer) > 0:
#             write_reads_to_reads_db(
#                 path_reads_database=path_reads_db,
#                 buffer_reads=buffer,
#                 timeout=timeout)
#         totally_missed_reads.update(readnames - saved_reads)
#     return totally_missed_reads


# def crs_to_dbs(
#         path_crs:Path,
#         path_connections:Path,
#         path_reads_db:Path,
#         path_read_intervals_db:Path,
#         path_sampledicts:Path,
#         threads:int,
#         threshold_connections:int,
#         crs_fetch_merge_distance:int,
#         tolerated_clipped_tails_length:int,
#         chunk_size:int,
#         output_mcrIDs:Path,
#         timeout:float) -> None:
#     # =============================================================================
#     # prepare
#     dict_sampleID_bam = util.load_sampledicts(input=path_sampledicts)[2]
#     _, UF, __ = read_connections(path_connections=path_connections,threshold=threshold_connections)
#     # create reads database with columns: readname, sampleID, sequence, qualities
#     create_reads_db(path_reads_db=path_reads_db,timeout=timeout)
#     # create intervals database with columns: readname, crID, representative, start, end, refID, refStart, refEnd, sampleID
#     create_intervals_db(path_intervals_db=path_read_intervals_db,timeout=timeout)
#     # =============================================================================
#     #  main loop to write reads and intervals to databases
#     log.info("collecting stretches of adjacent candidate regions..")
#     crs_iter = util.yield_from_crs(input=path_crs)
#     collected_crs = [next(crs_iter)]
#     jobs = []
#     while True:
#         next_cr = next(crs_iter,None)
#         # process existing container of crs
#         if not next_cr or not consecutive_regions(next_cr,collected_crs[-1],crs_fetch_merge_distance):
#             args = {
#                 'dict_sampleID_bam': dict_sampleID_bam,
#                 'crs': collected_crs,
#                 'dict_representatives': representative_crs(UF,[cr.crID for cr in collected_crs]),
#                 'path_reads_database': path_reads_db,
#                 'path_intervals_database': path_read_intervals_db,
#                 'max_chunk_size': chunk_size,
#                 'tolerated_clipped_tails_length': tolerated_clipped_tails_length,
#                 'timeout': timeout
#             }
#             jobs.append(args)

#             # reset containers to next cr
#             if next_cr:
#                 collected_crs = [next_cr]
#             else:
#                 collected_crs = []
#                 break
#         # next_cr is consecutive -> add to containers
#         elif next_cr:
#             collected_crs.append(next_cr)

#     # --- do the jobs in parallel --- #
#     time_start = time.time()
#     if threads > 1:
#         with Pool(threads) as pool:
#             missed_reads_list = list(pool.map(process_func, jobs))
#     else:
#         missed_reads_list = [process_crs_on_bams(**job) for job in jobs]

#     log.debug(f"finished writing to dbs jobs.")

#     # missed_reads_list is of the format [{sampleID: {readname,..},..},..]
#     # at first, all dicts are merged into one dict of the form {sampleID: {readname,..},..}
#     missed_reads = dict()
#     for missed_reads_dict in missed_reads_list:
#         if missed_reads_dict:
#             for sampleID,readnames in missed_reads_dict.items():
#                 if readnames:
#                     if sampleID in missed_reads:
#                         missed_reads[sampleID].update(readnames)
#                     else:
#                         missed_reads[sampleID] = readnames

#     log.debug(f"finished flattening missed reads.")

#     if len(missed_reads) > 0:
#         total_missed_reads = sum([len(reads) for reads in missed_reads.values()])
#         log.warning(f"{total_missed_reads} reads could not be found on the first sighting in the bam file(s). Need to re-iterate the bam file(s).")
#         # then, the missed reads are iterated.
#         totally_missed_reads = write_missed_reads_to_db(
#             dict_sampleID_bam=dict_sampleID_bam,
#             missed_reads=missed_reads,
#             path_reads_db=path_reads_db,
#             chunk_size=chunk_size,
#             timeout=timeout)
#         if len(totally_missed_reads) > 0:
#             raise(f"Some read sequences could not be extracted from the bam file. \
#     This means that for all {len(totally_missed_reads)} reads there is no complete read sequence stored in the bam file. \
#     Make sure to use minimap2 and if you filter the bam file, make sure to always remove all alignments of reads that are filtered. \
#     The read names are: {totally_missed_reads}")

#     print(f"finished in {time.time()-time_start:.2f} seconds")
#     log.info("finished writing to databases")
#     # get unique mcrIDs from the intervals database
#     # and write them to a file (output_mcrIDs)
#     con_intervals = sqlite3.connect(path_read_intervals_db)
#     cur_intervals = con_intervals.cursor()
#     cur_intervals.execute("SELECT DISTINCT mcrID FROM intervals")
#     mcrIDs = [int(mcrID[0]) for mcrID in cur_intervals.fetchall()]
#     con_intervals.close()
#     with open(output_mcrIDs,'w') as f:
#         for mcrID in sorted(mcrIDs):
#             print(mcrID,file=f)
#     log.info(f"written mcrIDs to {output_mcrIDs}")
#     return

# # %%
# # argparse
# def run(args,**kwargs):
#     crs_to_dbs(
#         path_crs=args.crs,
#         path_connections=args.connections,
#         path_reads_db=args.reads_db,
#         path_read_intervals_db=args.intervals_db,
#         path_sampledicts=args.sampledicts,
#         threads=args.threads,
#         threshold_connections=args.threshold_connections,
#         crs_fetch_merge_distance=args.crs_fetch_merge_distance,
#         tolerated_clipped_tails_length=args.tolerated_clipped_tails_length,
#         chunk_size=args.chunk_size,
#         output_mcrIDs=args.output_mcrIDs,
#         timeout=args.timeout)

# def get_parser():
#     parser = argparse.ArgumentParser(description='Creates a database with all reads and their intervals in candidate regions.')
#     parser.add_argument('-c','--crs', type=Path, required=True, help='input path to candidate regions file')
#     parser.add_argument('-n','--connections', type=Path, required=True, help='input path to candidate regions connections file')
#     parser.add_argument('-s','--sampledicts', type=Path, required=True, help='input path to sampledicts file')
#     parser.add_argument('-r','--reads-db', type=Path, required=True, help='output path to reads database')
#     parser.add_argument('-i','--intervals-db', type=Path, required=True, help='output path to intervals database')
#     parser.add_argument('-t','--threads', type=int, required=False, default=8, help='param number of threads to use')
#     parser.add_argument('-o','--output-mcrIDs', type=Path, required=True, help='output path to a txt file with mcrIDs. Used in the snakemake workflow consensus rule.')
#     parser.add_argument('--threshold-connections', type=int, required=False, default=2, help='param minimum number of connectiong reads between two candidate regions to be considered as connected')
#     parser.add_argument('--crs-fetch-merge-distance', type=int, required=False, default=20_000, help='param maximum distance between two consecutive candidate regions to be merged into one fetch region')
#     parser.add_argument('--tolerated-clipped-tails-length', type=int, required=False, default=100, help='param maximum length of clipped tails to be considered as clipped')
#     parser.add_argument('--chunk-size', type=int, required=False, default=1_000, help='param number of reads to be collected in memory before writing to database')
#     parser.add_argument('--timeout', type=float, required=False, default=300.0, help='param timeout in seconds for sqlite3 connections')
#     return parser

# def main():
#     parser = get_parser()
#     args = parser.parse_args()
#     run(args)
#     return

# if __name__ == "__main__":
#     main()
