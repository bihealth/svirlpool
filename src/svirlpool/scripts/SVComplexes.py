# #%%
# import itertools
# from enum import Enum
# import sqlite3
# from pathlib import Path
# import cattrs
# import logging as log
# from typing import Generator
# import tempfile
# from Bio import SeqIO
# from shlex import split
# import subprocess
# import pickle

# from intervaltree import Interval


# from svirlpool.scripts import datatypes, SVprimitive_class

# SQLITE_TIMEOUT = 60  # seconds, adjust as needed

# #%%

# # def print_group_intervals(group:list[SVprimitive_class.SVprimitive]) -> None:
# #     """Prints the intervals of a group of SVprimitives."""
# #     for svp in group:
# #         print(f"{svp.chr}:{svp.ref_start}-{svp.ref_end} :: read={svp.read_start}-{svp.read_end} :: reverse={svp.aln_is_reverse}")

# def SVprimitives_to_SVcomplexes(svPrimitives:list[SVprimitive_class.SVprimitive]) -> list[datatypes.SVcomplex]:
#     """Generates SVcomplexes from the SVprimitives and the dict_consensus_alignments."""
#     # at first sort all svPrimitives by their consensusID and read_start
#     breakends:list[SVprimitive_class.SVprimitive] = [svp for svp in svPrimitives if (svp.sv_type == 3 or svp.sv_type == 4)]
#     breakends.sort(key=lambda x: (x.consensusID, x.read_start))
#     # then iterate chunks of grouped svPrimitives by their consensusID
#     svComplexes:list[datatypes.SVcomplex] = []
#     for consensusID,group in itertools.groupby(breakends, key=lambda x: x.consensusID):
#         group:SVprimitive_class.SVprimitive = sorted(group, key=lambda x: x.read_start)
#         tworelations:dict[tuple[int,int],set[TWORELATIONS]] = two_relations_of_group(group)
#         fourrelations:dict[tuple[int,int,int,int],set[FOURRELATIONS]] = four_relations_of_group(group)
#         # call inversions
#         svComplexes = _add_SVcomplexes_to_list(higher=call_inversions(group, fourrelations), lower=svComplexes)
#         # call deletions
#         svComplexes = _add_SVcomplexes_to_list(higher=call_deletions(group, tworelations), lower=svComplexes)
#         # call insertions
#         svComplexes = _add_SVcomplexes_to_list(higher=call_insertions(group, tworelations), lower=svComplexes)
#     return svComplexes

# # apply after calling SVcomplexes to a group of SVprimitives
# # this function removes all SVprimitives that are part of the given SVcomplexes from the group
# # this is useful to avoid calling break ends that are already part of a complex SV
# def subtract_SVcomplexes_from_SVprimitives(
#         svPrimitives:list[SVprimitive_class.SVprimitive],
#         svComplexes:list[datatypes.SVcomplex]) -> list[SVprimitive_class.SVprimitive]:
#     """returns a list of SVprimitives that are not part of the given SVcomplexes."""
#     used_svps:list[SVprimitive_class.SVprimitive] = [svp for svc in svComplexes for svp in svc.svprimitives]
#     result:list[SVprimitive_class.SVprimitive] = []
#     for svp in svPrimitives:
#         if svp not in used_svps:
#             result.append(svp)
#     return result


# def _add_SVcomplexes_to_list(
#         higher:list[datatypes.SVcomplex],
#         lower:list[datatypes.SVcomplex]) -> list[datatypes.SVcomplex]:
#     """Adds the lower SVcomplexes to the higher SVcomplexes and returns the combined list."""
#     """However, no lower svComplex can be added to a higher svComplex if it already contains a SVprimitive that is part of any higher one"""
#     used_svps:list[SVprimitive_class.SVprimitive] = [svp for svc in higher for svp in svc.svprimitives]
#     result = []
#     for sv_complex in lower:
#         # check if any SVprimitive of the lower sv_complex is already used in a higher sv_complex
#         if any(svp in used_svps for svp in sv_complex.svprimitives):
#             continue
#         # add the lower sv_complex to the higher sv_complexes
#         result.append(sv_complex)
#     # add the higher sv_complexes to the result
#     result.extend(higher)
#     return result


# class TWORELATIONS(Enum):
#     READCLOSED = 0 # read closed BNDs are directly adjacent without read sequence betwen them, e.g. deletions, translocations
#     REFCLOSED = 1 # ref closed BNDs are directly adjacent without a gap in the reference sequence, e.g. inversions, insertions
#     INVERSION = 2 # inversions are two BNDs that have different read orientations
#     TRANSLOCATION = 3 # two BNDs don't share the same chr
#     OVERLAP = 4 # the alignments of the two BNDs overlap on the reference
#     REFGAP = 5 # two BNDs are separated by a gap in the reference sequence, e.g. deletions
#     READGAP = 6 # two BNDs are separated by a gap in the read sequence, e.g. insertions
#     def __lt__(self, other):
#         return self.value < other.value
#     def __eq__(self, other):
#         return self.value == other.value
#     def __hash__(self):
#         return hash(self.value)

# class FOURRELATIONS(Enum):
#     HOP = 0 # breakends 2,3 are outside of the reference interval 1,4 given 1,4 share the same chr.
#     INVERSION = 1 # reverse aln is equal between breakends 1,4 and 2,3 bot not between 1,2 and 3,4


# # two-relations allow the detection of translocations, insertions and deletions
# def two_relations_of_group(group:list[SVprimitive_class.SVprimitive], prox:int=30) -> dict[tuple[int,int],set[TWORELATIONS]]:
#     if len(group) < 2:
#         return {}
#     assert all(x.consensusID == group[0].consensusID for x in group), "Not all SVprimitives have the same consensusID"
#     assert all(x.sv_type == 3 or x.sv_type == 4 for x in group), "Not all SVprimitives are of type BND"
#     assert all(group[i].read_start <= group[i+1].read_start for i in range(len(group)-1)), "SVprimitives are not sorted by read_start"
#     assert all(x.alignment_to_ref is not None for x in group), "Not all SVprimitives have an alignment_to_ref"

#     two_relations:dict[tuple[int,int],set[TWORELATIONS]]={}
#     for i in range(len(group)-1):
#         a:SVprimitive_class.SVprimitive = group[i]
#         b:SVprimitive_class.SVprimitive = group[i+1]
#         # both BNDs need to be on different alignments
#         if a.alignment_to_ref == b.alignment_to_ref:
#             continue
#         tags:set[TWORELATIONS] = set()
#         # READCLOSED - two adjacent bnds are no farther apart than 'prox' bp on the read
#         if abs(a.read_start - b.read_start) <= prox:
#             tags.add(TWORELATIONS.READCLOSED)
#         # REFCLOSED - two adjacent bnds are no farther apart than 'prox' bp on the reference
#         if abs(a.ref_start - b.ref_start) <= prox:
#             tags.add(TWORELATIONS.REFCLOSED)
#         # INVERSION - two adjacent bnds have different aln_is_reverse
#         if a.aln_is_reverse != b.aln_is_reverse:
#             tags.add(TWORELATIONS.INVERSION)
#         # TRANSLOCATION - two adjacent bnds are not on the same chromosome
#         if a.chr != b.chr:
#             tags.add(TWORELATIONS.TRANSLOCATION)

#         if a.chr == b.chr:
#             a_pysam = a.alignment_to_ref.to_pysam()
#             b_pysam = b.alignment_to_ref.to_pysam()
#             a_it = Interval(a_pysam.reference_start, a_pysam.reference_end)
#             b_it = Interval(b_pysam.reference_start, b_pysam.reference_end)
#             # OVERLAP - two adjacent bnds overlap on the reference
#             if a_it.overlap_size(b_it) >= prox:
#                 tags.add(TWORELATIONS.OVERLAP)
#             # REFGAP - two adjacent bnds are separated by a gap in the reference sequence
#             elif not a_it.overlaps(b_it) and abs(a.ref_start - b.ref_start) > prox:
#                 tags.add(TWORELATIONS.REFGAP)
#         # READGAP - two adjacent bnds are separated by a gap in the read sequence
#         if abs(a.read_start - b.read_start) > prox:
#             tags.add(TWORELATIONS.READGAP)

#         two_relations[(i,i+1)] = tags
#     return two_relations


# # four-relations allow the detection of complex SVs like inversions
# # some insertions, e.g. of Mobile Elements, can be aligned to a different locus
# # the generalization is a "hop" where a middle aigned segment between two adjacent segments is translocated.
# def four_relations_of_group(group:list[SVprimitive_class.SVprimitive], prox:int=5) -> dict[tuple[int,int,int,int],set[FOURRELATIONS]]:
#     if len(group) < 4:
#         return {}
#     assert all(x.consensusID == group[0].consensusID for x in group), "Not all SVprimitives have the same consensusID"
#     assert all(x.sv_type == 3 or x.sv_type == 4 for x in group), "Not all SVprimitives are of type BND"
#     assert all(group[i].read_start <= group[i+1].read_start for i in range(len(group)-1)), "SVprimitives are not sorted by read_start"

#     four_relations:dict[tuple[int,int,int,int],set[FOURRELATIONS]]={}
#     for i in range(len(group)-3):
#         a:SVprimitive_class.SVprimitive = group[i]
#         b:SVprimitive_class.SVprimitive = group[i+1]
#         c:SVprimitive_class.SVprimitive = group[i+2]
#         d:SVprimitive_class.SVprimitive = group[i+3]
#         # breakend pairs a,b and c,d need to be on different alignments
#         if a.alignment_to_ref == b.alignment_to_ref \
#                 or c.alignment_to_ref == d.alignment_to_ref:
#             continue
#         tags:set[FOURRELATIONS] = set()
#         # HOP - breakends b,c are outside of the reference interval a,d given a,d share the same chr.
#         if a.alignment_to_ref != b.alignment_to_ref \
#                 and b.alignment_to_ref == c.alignment_to_ref \
#                 and c.alignment_to_ref != d.alignment_to_ref:
#             interval_ad = Interval(a.ref_start, d.ref_start)
#             if a.chr != b.chr:
#                 tags.add(FOURRELATIONS.HOP)
#             elif not interval_ad.overlaps(b.ref_start, b.ref_start) and not interval_ad.overlaps(c.ref_start, c.ref_start):
#                 tags.add(FOURRELATIONS.HOP)
#         # INVERSION - reverse aln is equal between breakends a,d and b,c but not between a,b and c,d
#         if a.aln_is_reverse == d.aln_is_reverse and b.aln_is_reverse == c.aln_is_reverse \
#                 and a.aln_is_reverse != b.aln_is_reverse:
#             tags.add(FOURRELATIONS.INVERSION)

#         four_relations[(i,i+1,i+2,i+3)] = tags
#     return four_relations


# def call_inversions(
#         group:list[SVprimitive_class.SVprimitive],
#         fourrelations:dict[tuple[int,int,int,int],set[FOURRELATIONS]]) -> list[datatypes.SVcomplex]:
#     """Calls inversions from a group of SVprimitives and their four-relations."""
#     # an inversion is present if there exists a 4-hop and a 4-inv
#     if len(group) < 4:
#         return []
#     samplename = group[0].vcfID.split("_")[1]
#     results:list[datatypes.SVcomplex] = []
#     for (a,b,c,d),x in fourrelations.items():
#         if x == set([FOURRELATIONS.HOP, FOURRELATIONS.INVERSION]):
#             # create a SVcomplex from the group
#             sv_complex = datatypes.SVcomplex(
#                 svType=5, # inversion
#                 svprimitives=[group[a], group[b], group[c], group[d]],
#                 vcfID=f"{datatypes.SV_TYPE_DICT[5]}_{samplename}_{group[a].consensusID}"
#             )
#             results.append(sv_complex)
#     return results


# def call_deletions(
#         group:list[SVprimitive_class.SVprimitive],
#         tworelations:dict[tuple[int,int],set[TWORELATIONS]]) -> list[datatypes.SVcomplex]:
#     """Calls deletions from a group of SVprimitives and their two-relations."""
#     # a deletion is present if there exists a 2-readclosed and a 2-REFGAP
#     if len(group) < 2:
#         return []
#     samplename = group[0].vcfID.split("_")[1]
#     results:list[datatypes.SVcomplex] = []
#     for (a,b),x in tworelations.items():
#         if x == set([TWORELATIONS.READCLOSED, TWORELATIONS.REFGAP]):
#             # create a SVcomplex from the group
#             sv_complex = datatypes.SVcomplex(
#                 svType=1,
#                 svprimitives=[group[a], group[b]],
#                 vcfID=f"{datatypes.SV_TYPE_DICT[1]}_{samplename}_{group[a].consensusID}"
#             )
#             results.append(sv_complex)
#     return results


# def call_insertions(
#         group:list[SVprimitive_class.SVprimitive],
#         tworelations:dict[tuple[int,int],set[TWORELATIONS]]) -> list[datatypes.SVcomplex]:
#     """Calls insertions from a group of SVprimitives and their two-relations."""
#     # an insertion is present if there exists a 2-refclosed and a 2-readgap
#     if len(group) < 2:
#         return []
#     samplename = group[0].vcfID.split("_")[1]
#     results:list[datatypes.SVcomplex] = []
#     for (a,b),x in tworelations.items():
#         if x == set([TWORELATIONS.REFCLOSED, TWORELATIONS.READGAP]):
#             # create a SVcomplex from the group
#             sv_complex = datatypes.SVcomplex(
#                 svType=2, # insertion
#                 svprimitives=[group[a], group[b]],
#                 vcfID=f"{datatypes.SV_TYPE_DICT[2]}_{samplename}_{group[a].consensusID}"
#             )
#             results.append(sv_complex)
#     return results


# def create_svComplexes_db(
#         database: Path,
#         name_reference: str):
#     """Creates a database for SVcomplexes."""
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.', '_')
#     with sqlite3.connect('file:' + str(database) + '?mode=rwc', uri=True) as conn:
#         conn.execute(f"DROP TABLE IF EXISTS svComplexes__{name_reference}")
#         conn.execute(f"""CREATE TABLE svComplexes__{name_reference} (
#                         svComplexID VARCHAR(30) PRIMARY KEY,
#                         consensusID VARCHAR(30),
#                         name_reference VARCHAR(50),
#                         svComplex BLOB)""")
#         conn.execute(f"CREATE INDEX idx_svComplexes_consensusID ON svComplexes__{name_reference} (consensusID)")
#         conn.commit()

# def write_svComplexes_to_db(
#         database: str,
#         data: list[datatypes.SVcomplex],
#         name_reference: str) -> None:
#     """Writes SVcomplexes to the database."""
#     if '.' in name_reference:
#         name_reference = name_reference.replace('.', '_')
#     cache = []
#     for svComplex in data:
#         if not svComplex.vcfID or svComplex.vcfID == '':
#             raise ValueError(f"svComplex has no valid vcfID: {svComplex.vcfID}")
#         consensusID = svComplex.svprimitives[0].consensusID
#         svComplexID = svComplex.vcfID
#         cache.append((svComplexID, consensusID, name_reference, pickle.dumps(svComplex.unstructure())))

#     if cache:
#         query = f"INSERT INTO svComplexes__{name_reference} (svComplexID, consensusID, name_reference, svComplex) VALUES (?, ?, ?, ?)"
#         with sqlite3.connect(database, timeout=SQLITE_TIMEOUT) as conn:
#             c = conn.cursor()
#             c.executemany(query, cache)
#             conn.commit()
#             c.close()

# def read_svComplexes_from_db(
#         database: Path,
#         reference_name: str,
#         consensusIDs: list[str] | None = None) -> list[datatypes.SVcomplex]:
#     """Reads SVcomplexes from the database."""
#     if '.' in reference_name:
#         reference_name = reference_name.replace('.', '_')
#     svComplexes:list[datatypes.SVcomplex] = []
#     with sqlite3.connect(str(database)) as conn:
#         table_name = f"svComplexes__{reference_name}"
#         c = conn.cursor()
#         try:
#             if consensusIDs:
#                 placeholders = ','.join(['?' for _ in consensusIDs])
#                 query = f"SELECT svComplex FROM {table_name} WHERE consensusID IN ({placeholders})"
#                 c.execute(query, [*consensusIDs,])
#             else:
#                 query = f"SELECT svComplex FROM {table_name}"
#                 c.execute(query)

#             for row in c.fetchall():
#                 # Use cattrs.structure to deserialize the binary blob back into an SVcomplex object
#                 svComplexes.append(cattrs.structure(pickle.loads(row[0]), datatypes.SVcomplex))
#         except sqlite3.Error as e:
#             log.error(f"Error reading SVcomplexes from database: {e}")
#             # print all tables found in the database
#             c.execute("SELECT name FROM sqlite_master WHERE type='table';")
#             tables = c.fetchall()
#             log.error(f"Tables in database: {tables}")
#             c.close()
#             raise e
#         c.close()
#     return svComplexes

# def yield_svComplexes_from_db(
#         database: Path,
#         n_rows: int,
#         reference_name: str) -> Generator[datatypes.SVcomplex, None, None]:
#     """Yields SVcomplexes from the database in chunks."""
#     if '.' in reference_name:
#         reference_name = reference_name.replace('.', '_')
#     with sqlite3.connect("file:" + str(database) + "?mode=ro", uri=True) as conn:
#         c = conn.cursor()
#         try:
#             c.execute(f"SELECT svComplex FROM svComplexes__{reference_name}")
#             while True:
#                 rows = c.fetchmany(n_rows)
#                 if not rows:
#                     break
#                 for row in rows:
#                     yield cattrs.structure(pickle.loads(row[0]), datatypes.SVcomplex)
#         except sqlite3.Error as e:
#             log.error(f"Error reading SVcomplexes from database: {e}")
#             # print all tables found in the database
#             c.execute("SELECT name FROM sqlite_master WHERE type='table';")
#             tables = c.fetchall()
#             log.error(f"Tables in database: {tables}")
#             raise e
#         c.close()


# def add_alt_sequences_to_svcomplexes_inplace(
#         SVcomplexes:list[datatypes.SVcomplex],
#         sequences_core:dict[str,bytes],
#         intervals_core:dict[str,tuple[int,int]]) -> None:
#     """Adds alt sequences to the SVcomplexes."""
#     # thoughts
#     # 1) each svcomplex has oly one consensusID in all its constituting svprimitives
#     # 2) the alt sequence interval can be taken from the start positions of all svprimitives
#     for svcomplex in SVcomplexes:
#         consensusID = svcomplex.svprimitives[0].consensusID
#         if consensusID not in sequences_core:
#             raise(f"ConsensusID {consensusID} not found in sequences_core, cannot add alt sequence to SVcomplex {svcomplex.vcfID}")
#         if consensusID not in intervals_core:
#             raise(f"ConsensusID {consensusID} not found in intervals_core, cannot add alt sequence to SVcomplex {svcomplex.vcfID}")
#         if svcomplex.svType == 0: # insertion
#             start = min(svcomplex.svprimitives[0].read_start, svcomplex.svprimitives[1].read_start)
#             end = max(svcomplex.svprimitives[0].read_end, svcomplex.svprimitives[1].read_end)
#             sequence = pickle.loads(sequences_core[consensusID])[start:end]
#             svcomplex.alt_sequences = [sequence]
#         elif svcomplex.svType == 5: # inversion
#             start = min(svcomplex.svprimitives[1].read_start, svcomplex.svprimitives[2].read_start)
#             end = max(svcomplex.svprimitives[1].read_end, svcomplex.svprimitives[2].read_end)
#             sequence = pickle.loads(sequences_core[consensusID])[start:end]
#             svcomplex.alt_sequences = [sequence]


# def region_string_from_svcomplex(svcomplex:datatypes.SVcomplex) -> str:
#     """Returns a region string for the SVcomplex."""
#     if svcomplex.svType == 1: # deletion
#         chr = svcomplex.svprimitives[0].chr
#         start = min(svcomplex.svprimitives[0].ref_start, svcomplex.svprimitives[1].ref_start)
#         end = max(svcomplex.svprimitives[0].ref_start, svcomplex.svprimitives[1].ref_start)
#         return f"{chr}:{start}-{end}"
#     elif svcomplex.svType == 2: # insertion
#         chr = svcomplex.svprimitives[0].chr
#         start = min(svcomplex.svprimitives[0].ref_start, svcomplex.svprimitives[1].ref_start)
#         end = max(svcomplex.svprimitives[0].ref_start, svcomplex.svprimitives[1].ref_start)
#         return f"{chr}:{start}-{end}"
#     elif svcomplex.svType == 5: # inversion
#         chr = svcomplex.svprimitives[1].chr
#         start = min(svcomplex.svprimitives[1].ref_start, svcomplex.svprimitives[2].ref_start)
#         end = max(svcomplex.svprimitives[1].ref_start, svcomplex.svprimitives[2].ref_start)
#         return f"{chr}:{start}-{end}"
#     else:
#         raise ValueError(f"Unsupported SVtype {svcomplex.svType} for region string generation")

# def region_string_split(region_string:str) -> tuple[str,int,int]:
#     """Converts a region string to a BED line."""
#     chr, pos = region_string.split(':')
#     start, end = pos.split('-')
#     return (chr, int(start), int(end))


# def add_ref_sequences_to_svcomplexes_inplace(
#         SVcomplexes:list[datatypes.SVcomplex],
#         path_reference:Path|str) -> None:
#     """Adds ref sequences to the SVcomplexes."""
#     # generate regions
#     regions:list[str] = []
#     for svcomplex in SVcomplexes:
#         if svcomplex.svType == 1: # deletion
#             regions.append(region_string_from_svcomplex(svcomplex))
#     # save regions to tmp regions file
#     with tempfile.TemporaryDirectory() as tdir:
#         bed_path = Path(tdir) / "regions.bed"
#         with open(bed_path, 'w') as f:
#             for region in regions:
#                 bed_line:tuple[str,int,int] = region_string_split(region)
#                 print(*bed_line, sep='\t', file=f)
#         # use bedtools getfasta to get the sequences
#         cmd_getfasta = split(f"bedtools getfasta -fi {path_reference} -bed {str(bed_path)} -fo {tdir}/ref_sequences.fasta")
#         tmp_alt_fasta = Path(tdir) / "ref_sequences.fasta"
#         with open(tmp_alt_fasta, 'w') as f:
#             log.info(f"extracting ref sequences into {tmp_alt_fasta}...")
#             subprocess.check_call(cmd_getfasta, stdout=f)
#         # load fasta file # record name is the region string. ref_sequences is a dict of region:sequence
#         ref_sequences = {record.name: str(record.seq) for record in SeqIO.parse(tmp_alt_fasta, "fasta")}
#     # iterate over all svcomplexes and add the sequences.
#     for svcomplex in SVcomplexes:
#         if svcomplex.svType == 1:
#             # add ref sequences (ignoring original signals, just add for each svcomplex the corresponding ref sequence)
#             region = region_string_from_svcomplex(svcomplex)
#             svcomplex.ref_sequences = [ref_sequences[region].upper()]
#     # done


# # The functions below are examples from the sv primitive merging script.
# # We need to implement a similar logic for the merging of svComplexes.
# # They are composits of svprimitives, but they need their own merging logic,
# # that is: rules to determine of a pair of svComplexes can be merged.
# # and a function to merge a pair of svComplexes.

# def can_merge_svComplexes_pair(
#         first:datatypes.SVcomplex,
#         other:datatypes.SVcomplex,
#         far:int=10_000,
#         close:int=1_000,
#         very_close:int=50,
#         max_cohens_d:float=-2.0,
#         verbose:bool=False) -> bool:
#     """Checks if two SVcomplexes can be merged."""
#     assert isinstance(first, datatypes.SVcomplex), "first must be an instance of SVcomplex"
#     assert isinstance(other, datatypes.SVcomplex), "other must be an instance of SVcomplex"
#     # check if both are initialized by checking if their svprimitives are not empty lists
#     assert len(first.svprimitives) > 0, "first must be initialized with svprimitives"
#     assert len(other.svprimitives) > 0, "other must be initialized with svprimitives"
#     # two deletions can be transformed into a svprimitive of type 1.
#     # the catch is that both svprimitives' genotype measurements need to be merged somehow.


# from copy import deepcopy

# def merge_svPrimitives(svPrimitives:list[SVprimitive_class.SVprimitive],
#                        far:int=10_000,
#                        close:int=1_000,
#                        very_close:int=50,
#                        max_cohens_d:float=-2.0,
#                        verbose:bool=False) -> list[SVprimitive_class.SVprimitive]:
#     # --- merge svPrimitives that are most likely the same --- #
#     #svps:list[SVprimitive_class.SVprimitive|None] = deepcopy(sorted(svPrimitives,key=lambda x: (x.chr,x.ref_start,x.ref_end)))
#     SVPs = deepcopy(svPrimitives) # don't change the input
#     for i in range(len(SVPs)-1,-1,-1):
#         svPrimitive = SVPs[i]
#         if svPrimitive is None:
#             continue
#         j = i-1
#         # if svPrimitive.sv_type <= 1:
#         #     first_kmer_sketch = kmer_sketch_of_svPrimitive(svPrimitive,k=5)
#         while j >= 0:
#             other_svPrimitive = SVPs[j]
#             if other_svPrimitive is None:
#                 j -= 1
#                 continue
#             elif svPrimitive.chr != other_svPrimitive.chr:
#                 j -= 1
#                 continue
#             elif svPrimitive.sv_type != other_svPrimitive.sv_type:
#                 j -= 1
#                 continue

#             # TODO: don't merge if the k-mer sketches are very different.
#             # count kmers that are more frequent than 5% of the total kmers in both sketches
#             # if svPrimitive.sv_type <= 1:
#             #     other_kmer_sketch = kmer_sketch_of_svPrimitive(other_svPrimitive,k=5)
#             # else:
#             #     first_kmer_sketch = None
#             #     other_kmer_sketch = None
#             can_merge = can_merge_svPrimitives_pair(
#                 first=svPrimitive,
#                 other=other_svPrimitive,
#                 far=far,
#                 close=close,
#                 very_close=very_close,
#                 max_cohens_d=max_cohens_d,
#                 verbose=verbose)
#                 # first_kmer_sketch=first_kmer_sketch,
#                 # other_kmer_sketch=other_kmer_sketch)
#             if can_merge:
#                 SVPs[i] = merge_svPrimitives_pair(svPrimitive,other_svPrimitive)
#                 SVPs[j] = None
#                 svPrimitive = SVPs[i]
#             j -= 1
#     return [svPrimitive for svPrimitive in SVPs if svPrimitive is not None]


# def merge_svPrimitives_pair(first:SVprimitive_class.SVprimitive,second:SVprimitive_class.SVprimitive) -> SVprimitive_class.SVprimitive:
#     # check if first and second are fully initialized SVprimitives:

#     assert isinstance(first,SVprimitive_class.SVprimitive), "first must be an instance of SVprimitive"
#     assert isinstance(second,SVprimitive_class.SVprimitive), "second must be an instance of SVprimitive"
#     # check if both are initialized by checking if their original_signals are not empty lists
#     assert len(first.original_signals) > 0, "first must be initialized"
#     assert len(second.original_signals) > 0, "second must be initialized"
#     assert len(first.original_cr_signals) > 0, "first must be initialized with SV signals from the candidate region"
#     assert len(second.original_cr_signals) > 0, "second must be initialized with SV signals from the candidate region"
#     assert len(first.genotypeMeasurements) > 0, "first must be initialized with genotypeMeasurements"
#     assert len(second.genotypeMeasurements) > 0, "second must be initialized with genotypeMeasurements"

#     winner = first if max_support_of_any_samplename(first.genotypeMeasurements) >= max_support_of_any_samplename(second.genotypeMeasurements) else second
#     new_genotypeMeasurements = merge_genotypeMeasurements(first.genotypeMeasurements,second.genotypeMeasurements)
#     new_SV_primitive = deepcopy(first) if first == winner else second
#     new_SV_primitive.genotypeMeasurements = new_genotypeMeasurements
#     # merge original signals
#     new_SV_primitive.original_signals = first.original_signals + second.original_signals
#     # merge repeatIDs
#     new_SV_primitive.repeatIDs = sorted(set(first.repeatIDs + second.repeatIDs))
#     new_SV_primitive.reference_name = first.reference_name
#     new_SV_primitive.distortion = first.distortion + second.distortion
#     assert type(new_SV_primitive.reference_name) == str, "reference_name must be a string"
#     return new_SV_primitive


# def can_merge_svPrimitives_pair(
#         first:SVprimitive_class.SVprimitive,
#         other:SVprimitive_class.SVprimitive,
#         first_kmer_sketch:set|None=None,
#         other_kmer_sketch:set|None=None,
#         min_difference:int=8,
#         very_close:int=50,
#         close:int=1_000,
#         far:int=10_000,
#         max_cohens_d=-2.0,
#         verbose:bool=False) -> bool:
#     assert first.size >= 0 and other.size >= 0, "sizes must be >= 0"
#     if verbose:
#         print(f"trying to merge {first.sv_type} {first.size} {first.ref_start} with {other.sv_type} {other.size} {other.ref_start}")
#     if first.sv_type != other.sv_type:
#         if verbose:
#             print(f"can't merge because of different sv_types: {first.sv_type} and {other.sv_type}")
#         return False

#     ### --- precompute all features --- ###
#     first_is_bnd = first.sv_type > 1
#     other_is_bnd = other.sv_type > 1
#     distance = Interval(first.ref_start,first.ref_end).distance_to(Interval(other.ref_start,other.ref_end))
#     tolerance_a = int(round(np.mean(first.distortion)))
#     tolerance_b = int(round(np.mean(other.distortion)))
#     shared_reapeatIDs = set(first.repeatIDs).intersection(set(other.repeatIDs))
#     if verbose:
#         print(f"distance: {distance}, tolerance_a: {tolerance_a}, tolerance_b: {tolerance_b}, shared_reapeatIDs: {shared_reapeatIDs}, first_is_bnd: {first_is_bnd}, other_is_bnd: {other_is_bnd}")

#     ### --- RULES --- ###

#     def _bnd_merge() -> bool:
#         if first_is_bnd and other_is_bnd:
#             if distance <= very_close:
#                 if verbose:
#                     print(f"can merge BNDs because of distance: {distance} <= {very_close}")
#                 return True
#             else:
#                 if verbose:
#                     print(f"can't merge BNDs because of distance: {distance} > {very_close}")
#                 return False
#         return False

#     # ================================ 1) compare distances ================================ #

#     def _distance_is_ok() -> bool:
#         close_far = far if len(shared_reapeatIDs) > 0 else close
#         if distance > tolerance_a + tolerance_b + close_far + min_difference:
#             if verbose:
#                 print(f"can't merge because of distance: {distance} > {tolerance_a} + {tolerance_b} + {close_far} + {min_difference}")
#             return False
#         else:
#             if verbose:
#                 print(f"can merge because of distance: {distance} <= {tolerance_a} + {tolerance_b} + {close_far} + {min_difference}")
#             return True
#     # ================================ 2) compare sizes ================================ #

#     # adjust size estimates based on cohen's d if both distortion_sum_per_read have more than 2 values
#     # if only one has more than 2 values, test if the other is less than 2 stds away
#     # if both have less than 3 values, use the mean of the two and merge if their size difference is less than 8 or less than 2% of the smaller size
#     def _similar_size() -> bool:
#         if len(first.distortion_sum_per_read) > 2 and len(other.distortion_sum_per_read) > 2:
#             # if cohen's d of the sizes + values of distortion_sum_per_read is less than -1, they can't be merged
#             dist_a = np.array(list(first.distortion_sum_per_read.values())) + (-first.size if first.sv_type == 1 else first.size)
#             dist_b = np.array(list(other.distortion_sum_per_read.values())) + (-other.size if other.sv_type == 1 else other.size)
#             d = cohen_d(dist_a,dist_b)
#             if d < max_cohens_d:
#                 if verbose:
#                     print(f"can't merge because of cohen's d: {d} < {max_cohens_d}")
#                     print(f"distortion of A is {dist_a}, distortion of B is {dist_b}")
#                 return False
#             else:
#                 if verbose:
#                     print(f"can merge because of cohen's d: {d} >= {max_cohens_d}")
#                     print(f"distortion of A is {dist_a}, distortion of B is {dist_b}")
#                 return True
#         elif len(first.distortion_sum_per_read) > 2 or len(other.distortion_sum_per_read) > 2:
#             if len(first.distortion_sum_per_read) > 2:
#                 a = first
#                 b = other
#             else:
#                 a = other
#                 b = first
#             # check if b is in 2 stds of a
#             std_a = np.std(list(a.distortion_sum_per_read.values()))
#             if abs(b.size - a.size) > 2*std_a:
#                 if verbose:
#                     print(f"can't merge because of size: {abs(b.size - a.size)} > 2*{std_a}")
#                 return False
#             else:
#                 if verbose:
#                     print(f"can merge because of size: {abs(b.size - a.size)} <= 2*{std_a}")
#                 return True
#         else:
#             if not (abs(first.size - other.size) <= min_difference or abs(first.size - other.size) <= 0.02 * min(first.size,other.size)):
#                 if verbose:
#                     print(f"can't merge because of size: {abs(first.size - other.size)} > {min_difference} or {abs(first.size - other.size)} > 0.02 * {min(first.size,other.size)}")
#                 return False
#             else:
#                 if verbose:
#                     print(f"can merge because of size: {abs(first.size - other.size)} <= {min_difference} or {abs(first.size - other.size)} <= 0.02 * {min(first.size,other.size)}")
#                 return True

#     # ================================ 3) compare k-mers ================================ #
#     # ------ not used right now, because the giab benchmark worsened ----- #
#     # -- works not so well and is very slow so far -- #
#     # if kmer sketches are provided, first and other can only be merged, if they have similar kmer sketches
#     def _kmer_similarity() -> bool:
#         if first_kmer_sketch and other_kmer_sketch and len(first_kmer_sketch) > 0 and len(other_kmer_sketch) > 0:
#             # of one sketch is composed of only a few kmers, it is not reliable
#             sketch_ratio:float = 0.0
#             if len(first_kmer_sketch) < 10 and len(other_kmer_sketch) < 10:
#                 # both small
#                 sketch_ratio = util.intersection_ratio_of_smaller_set(first_kmer_sketch,other_kmer_sketch)
#                 MIN_RATIO = 0.3
#                 if sketch_ratio > MIN_RATIO:
#                     if verbose:
#                         print(f"can merge because of kmer sketch: {sketch_ratio} > {MIN_RATIO}")
#                     return True
#                 else:
#                     if verbose:
#                         print(f"can't merge because of kmer sketch: {sketch_ratio} <= {MIN_RATIO}")
#                     return False
#             if len(first_kmer_sketch) < 10 and len(other_kmer_sketch) > 10:
#                 sketch_ratio = util.intersection_ratio_of_smaller_set(first_kmer_sketch,other_kmer_sketch)
#                 MIN_RATIO = 0.45
#                 if sketch_ratio > MIN_RATIO:
#                     if verbose:
#                         print(f"can merge because of kmer sketch: {sketch_ratio} > {MIN_RATIO}")
#                     return True
#                 else:
#                     if verbose:
#                         print(f"can't merge because of kmer sketch: {sketch_ratio} <= {MIN_RATIO}")
#                     return False
#             else: # both large
#                 MIN_RATIO = 0.6
#                 sketch_ratio = util.intersection_ratio_of_smaller_set(first_kmer_sketch,other_kmer_sketch)
#                 if sketch_ratio > MIN_RATIO:
#                     if verbose:
#                         print(f"can merge because of kmer sketch: {sketch_ratio} > {MIN_RATIO}")
#                     return True
#                 else:
#                     if verbose:
#                         print(f"can't merge because of kmer sketch: {sketch_ratio} <= {MIN_RATIO}")
#                     return False
#         else:
#             if verbose:
#                 print(f"can merge: no kmer sketches provided")
#             return True

#     # ================================ CHECK CONDITIONS ================================ #

#     similar_size = _similar_size()
#     kmer_similarity = _kmer_similarity()
#     distance_ok = _distance_is_ok()
#     bnd_merge = _bnd_merge()


#     if bnd_merge:
#         if verbose:
#             print(f"can merge because of BND merge rules")
#         return True
#     if similar_size and distance_ok and kmer_similarity:
#         if verbose:
#             print(f"can merge because of size, distance and kmer similarity")
#         return True
#     else:
#         if verbose:
#             print(f"can't merge because of size, distance or kmer similarity")
#         return False
