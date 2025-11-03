# %%
# debug file to test consensus


# %%
# base_path = Path("/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/platinum/svirlpool/NA12889")

# original_cmd = """python3 -m svirlpool.scripts.consensus         -i NA12889.crs_containers.db         -s NA12889.sampledicts.json         -o consensus/0/consensus.4.txt         -t 2         -c 4         --logfile consensus/0/consensus.4.log"""
# # parse original command to a dict. if the word starts with "-" or with "--" then it is a key, else it is a value
# cmd_dict = {}
# key = None
# for word in original_cmd.split()[1:]:
#     if word.startswith("-"):
#         key = word
#     else:
#         cmd_dict[key] = word

# # test command
# crs_containers_to_consensus(
#     input=base_path / cmd_dict['-i'],
#     sampledicts=base_path / cmd_dict['-s'],
#     output=base_path / cmd_dict['-o'],
#     pairwise_distance_clustering_factor =2.0,
#     tolerance_factor_for_aligning_misfits=0.1,
#     threads=1,
#     crIDs=[int(cmd_dict['-c'])])

# %%
import pysam

from ..scripts import util


def test_get_interval_on_read_in_region__platinum_0():
    aln = None
    for a in pysam.AlignmentFile(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/tests/data/cut_reads/test.platinum.bam"
    ):
        if a.query_name == "886c55ec-a679-40af-9d8f-58e43f4c8a7a":
            aln = a
            break
    region_start = 21112695
    region_end = 21122675
    util.get_interval_on_read_in_region(a=aln, start=region_start, end=region_end)


test_get_interval_on_read_in_region__platinum_0()
# %%
