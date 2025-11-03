# %%

import sys
from pathlib import Path

# Add the parent directory to sys.path to enable imports from scripts
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.consensus_align import align_padded_consensus_sequences

# %%


# def align_padded_consensus_sequences(
#         input_consensus_container_results:Path|str,
#         input_reference:Path|str,
#         path_bamout:Path|str,
#         path_fastaout:Path|str,
#         threads:int,
#         overwrite_fasta:bool=False,
#         tmp_dir_path:Path|None=None) -> tuple[dict[str,bytes], dict[str,tuple[int,int]]]:
def run():
    input_consensus_container_results = Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/bench_approx/consensus_containers.txt"
    )
    input_reference = Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/references/hs37d5/hs37d5.fa"
    )
    path_bamout = Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/bench_approx/consensus.bam"
    )
    path_fastaout = Path(
        "/data/cephfs-1/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/giab/test/bench_approx/consensus.fasta"
    )
    threads = 1
    overwrite_fasta = True
    tmp_dir_path = None

    align_padded_consensus_sequences(
        input_consensus_container_results=input_consensus_container_results,
        input_reference=input_reference,
        path_bamout=path_bamout,
        path_fastaout=path_fastaout,
        threads=threads,
        overwrite_fasta=overwrite_fasta,
        tmp_dir_path=tmp_dir_path,
    )


if __name__ == "__main__":
    run()
