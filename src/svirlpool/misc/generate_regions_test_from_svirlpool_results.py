# load consensus objects and find the consensus objects that produced the most consensu sequences
# %%
from tqdm import tqdm

from ..scripts.consensus_align import parse_crs_container_results

# %%

path = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/HG002.30x/consensus_containers.txt"
path_reference = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa"
output = "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/test_regions/HG002_fp_regions.txt"

# at first, read the number of lines in the file
n_lines = sum(1 for line in tqdm(open(path)))
print(f"Number of lines in the file: {n_lines}")

threshold: int = 6

dict_counts = {}
dict_regions = {}

for crs_container_result in tqdm(parse_crs_container_results(path), total=n_lines):
    for consensusID, consensus_object in crs_container_result.consensus_dicts.items():
        cID = consensusID.split(".")[0]
        dict_counts[cID] = dict_counts.get(cID, 0) + 1
        if cID not in dict_regions:
            dict_regions[cID] = []
        dict_regions[cID].append(consensus_object.original_regions)
# %%
# sample the consensusIDs that produced 10 to 6 consensus sequences
list_selected = [cID for cID, count in dict_counts.items() if threshold >= count >= 6]
print(
    f"Number of consensusIDs with at least {threshold} consensus sequences: {len(list_selected)}"
)
# %%
import tempfile

# pick the first region of each consensusID in list_selected
# print to a bed file
tmp_bed = tempfile.NamedTemporaryFile(delete=True, suffix=".bed")
with open(tmp_bed.name, "w") as bed_file:
    for cID in list_selected:
        regions = dict_regions[cID]
        print(*regions[0][0], sep="\t", file=bed_file)
# %%
# sort it with bedtools
cmd_sort = f"bedtools sort -i {tmp_bed.name} -faidx {str(path_reference)}.fai"
import subprocess
from shlex import split

with open(output, "w") as sorted_bed_file:
    subprocess.run(split(cmd_sort), stdout=sorted_bed_file)
print(f"Sorted bed file written to {output}")
# %%
