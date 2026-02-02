# %%
import json
import pickle
from gzip import open as gzopen
from pathlib import Path

import cattrs
from intervaltree import Interval, IntervalTree

from svirlpool.svcalling.multisample_sv_calling import (
    SVcall, load_svCalls_from_json, load_svComposites_from_json,
    write_svCalls_to_vcf)
from svirlpool.svcalling.SVcomposite import SVcomposite

#%%
test_data_path = (
    Path(__file__).parent / "data" / "SVpatterns" / "svpattern_INV_parsing.json.gz"
)
#%%
path_merged_svComposites = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/tmp/all_svComposites.pkl"
with open(path_merged_svComposites, "rb") as f:
    merged_svComposites: list[SVcomposite] = pickle.load(f)

# %%
#def test_SVcalls_to_vcf() -> None:
path_SVcalls = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/tmp/svCalls.json.gz"
path_ref_bases_dict = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/tmp/ref_bases_dict.pkl"
path_covtrees = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/tmp/covtrees.pkl"

svCalls: list[SVcall] = load_svCalls_from_json(path_SVcalls)
with open(path_ref_bases_dict, "rb") as f:
    ref_bases_dict:dict[str, str] = pickle.load(f)
with open(path_covtrees, "rb") as f:
    covtrees:dict[str, dict[str, IntervalTree]] = pickle.load(f)
# %%
# e.g. SVcall(genotypes={'HG002': Genotype(samplename='HG002', genotype='0/0', gt_likelihood=1.0, genotype_quality=60, total_coverage=252, ref_reads=252, var_reads=0)}, passing=False, chrname='Y', end=58976098, start=58976097, svtype='INS', svlen=288, pass_altreads=False, pass_gq=True, precise=False, mateid='', consensusIDs=['HG002:1311.3'], ref_sequence=None, alt_sequence=b"\x80\x04\x95'\x01\x00\x00\x00\x00\x00\x00X \x01\x00\x00tctattccattccattccattccattccattccactgcattccattccattcctttccattccgttccatttcactccaccttctccgatccaatccatttcataccatccaattccatgattccactccattccattcataactttctacaagatctcactgtgtaacccaggctggtgagaggttctcactctgtcacccaggctggagtgtggtggcacaatatcacctgtcacacaggctggagtgcagcagcactatcttagctctggatgtcactctgtcac\x94.", sequence_id=None, description=None)
# check if a variant of the consensusID 193.0 is in svCalls
for svcall in svCalls:
    if "HG002:191.0" in svcall.consensusIDs:
        print(svcall)
# %%

from svirlpool.localassembly import SVpatterns

path_patterns_db = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/svpatterns.db"
# load sv patterns by cr ID
crID = 191
svps = SVpatterns.read_svPatterns_from_db(database=Path(path_patterns_db))

# %%
test = SVpatterns.read_svPatterns_from_db(database=Path(path_patterns_db), crIDs={crID})
# %%

# check the crs database "crs.db" in "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/crs.db"
# if the the candidate region 191 is present
# use:
# crs_dict = {
#         cr.crID: cr for cr in signalstrength_to_crs.load_crs_from_db(path_db=path_crs)
#     }
from svirlpool.candidateregions import signalstrength_to_crs

path_crs = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/HQ/svirlpool/HG002/crs.db"
crs_dict = {
        cr.crID: cr for cr in signalstrength_to_crs.load_crs_from_db(path_db=Path(path_crs))
    }
if 191 in crs_dict:
    print(f"crID 191 is present in crs.db: {crs_dict[191]}")
# %%
