# %%
# This script builds a project folder for a test case to use the LRSV detection workflow.
# Can be removed for launch.
# It allows to work on a subset of reads and alignments with given bed file but on several samples.

# =============================================================================
#  imports

import subprocess
from pathlib import Path
from shlex import split

import yaml
from logzero import logger as log

from . import util

# =============================================================================
#   user input
base_dir_name = "dels5k"
path_snakefile = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/src/workflows/main.smk"
)
path_project = Path(
    f"/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/{base_dir_name}"
)
path_bed = Path(
    f"/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/{base_dir_name}/targets.bed"
)
paths_bams = [
    Path(
        "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/alignments/minimap2.ont-r10.32x.bam"
    )
]
#    Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/alignments/minimap2.ont-r10.32x.alignments.bam")]
#    Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG003/alignments/ont.32x.minimap2.hs37d5.bam")]
paths_fastqs = [
    Path(
        "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/reads/ont-r10.32x.fastq"
    )
]
#    Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/reads/ont-r10.32x.reads.fastq")]
#    Path("/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG003/reads/GM24149.fastq")]
path_reference = Path(
    "/fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa"
)
path_reference_final = Path(
    "/fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa"
)
path_repeats = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/tools/pbsv/annotations/human_hs37d5.trf.bed"
)
path_lamassemble_matrix = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/tools/lamassemble/train/promethion.mat"
)
path_mappability = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/wgEncodeCrgMapabilityAlign50mer.lt_1.centered.bed"
)
benchmark = Path(
    "/fast/work/groups/cubi/projects/2022-10-18_May_LRSV-detection/development/HG/HG002/benchmarks/HG002_SVs_Tier1_v0.6.deletions.vcf.gz"
)

custom_names = ["HG002"]  # , "HG003"]
# used to sub-sample alignments (and reads)
subsamples = [1.0]  # ,1.0]

min_signal_size = 6
min_sv_size = 30

# =============================================================================
#   automatic paths

path_config_file = path_project / "wf_config.yml"
path_shell_script = path_project / "run_wf.sh"

# build dicts
dict_alias_bam = {custom_names[i]: p for i, p in enumerate(paths_bams)}
dict_alias_fastq = {custom_names[i]: p for i, p in enumerate(paths_fastqs)}
dict_subsamples = {n: subsamples[i] for i, n in enumerate(custom_names)}

paths_selected_bams = {
    k: path_project / "alignments" / f"{k}.bam" for k in custom_names
}
paths_selected_reads = {
    k: path_project / "reads" / f"{k}.readIDs.txt" for k in custom_names
}
# %%
# =============================================================================
#   setup directory structure

# create a folder for the project
# and create the subfolders: alignments, reads, QC, consensud, fastq
log.info(f"Creating project folder {path_project}")
path_project.mkdir(exist_ok=True)
for dirname in ["alignments", "reads", "QC", "consensus", "fastq"]:
    (path_project / dirname).mkdir(exist_ok=True, parents=True)

# copy bed file to project path
# check before if target file does not yet exist
if not (path_project / path_bed.name).exists():
    subprocess.check_call(split(f"cp {path_bed} {path_project}/{path_bed.name}"))

# %%
# # =============================================================================
# #   process reads and alignments

# use samtools to extract alignments from each bam file and copy these into the subfolder alignments
log.info("Extracting alignments from bam files")
for k, v in dict_alias_bam.items():
    if not paths_selected_bams[k].exists():
        log.info(f"Extracting alignments from {v}")
        cmd = f"samtools view -s {dict_subsamples[k]} -b -L {path_bed} {v}"
        with open(paths_selected_bams[k], "wb") as f:
            subprocess.check_call(split(cmd), stdout=f)
        # index the bam file
        if path_reference == path_reference_final:
            cmd_index = f"samtools index {paths_selected_bams[k]}"
            subprocess.check_call(split(cmd_index))

# then extract all unique read ids from the bam files and save them to txt files in the reads folder
log.info("Extracting read ids from bam files")
for k, p_bam in paths_selected_bams.items():
    if not paths_selected_reads[k].exists():
        log.info(f"Extracting read ids from {p_bam}")
        cmd = f"samtools view {p_bam}"
        cmd_cut = "cut -f1"
        cmd_sort = "sort"
        cmd_uniq = "uniq"
        with open(paths_selected_reads[k], "w") as f:
            p0 = subprocess.Popen(split(cmd), stdout=subprocess.PIPE)
            p1 = subprocess.Popen(
                split(cmd_cut), stdin=p0.stdout, stdout=subprocess.PIPE
            )
            p2 = subprocess.Popen(
                split(cmd_sort), stdin=p1.stdout, stdout=subprocess.PIPE
            )
            p3 = subprocess.Popen(split(cmd_uniq), stdin=p2.stdout, stdout=f)
            p3.communicate()
        # print number of lines in paths_selected_reads[k] to log
        log.info(
            f"Number of reads in {paths_selected_reads[k]}: {sum(1 for line in open(paths_selected_reads[k]))}"
        )

log.info("Extracting reads from fastq files")
# then extract the reads from the fastq files and save them to the reads folder
for k, v in dict_alias_fastq.items():
    if not (path_project / "reads" / f"{k}.fastq").exists():
        log.info(f"Extracting reads from {v}")
        p_out = path_project / "reads" / f"{k}.fastq"
        cmd = f"seqtk subseq {dict_alias_fastq[k]} {path_project / 'reads' / f'{k}.readIDs.txt'}"
        with open(p_out, "w") as f:
            subprocess.check_call(split(cmd), stdout=f)

# remove alignments and re-align fastq files to reference
if path_reference != path_reference_final:
    log.info("Re-aligning reads to reference")
    for k, v in dict_alias_fastq.items():
        if not (path_project / "alignments" / f"{k}.bam").exists():
            log.info(f"Re-aligning {v}")
            p_out = path_project / "alignments" / f"{k}.bam"
            p_in = path_project / "reads" / f"{k}.fastq"
            util.align_reads_with_minimap(
                reference=path_reference,
                reads=p_in,
                bamout=p_out,
                threads=8,
                tech="map-ont",
            )

# =============================================================================
#   process workflow script and config

# create config yml file for workflow
if not path_config_file.exists():
    log.info(f"Creating config file {path_config_file}")
    with open(path_config_file, "w") as f:
        yaml.dump(
            {
                "workdir": str(path_project),
                "reference": str(path_reference),
                "reference_final": str(path_reference_final),
                "prefix": path_project.name,
                "trf": str(path_repeats),
                "mat": str(path_lamassemble_matrix),
                "reads": [str(v) for k, v in dict_alias_fastq.items()],
                "alignments": [str(v) for k, v in paths_selected_bams.items()],
                "min_signal_size": min_signal_size,
                "min_sv_size": min_sv_size,
                "mappability": str(path_mappability),
                "regions": str(path_bed),
                "benchmark": str(benchmark),
            },
            f,
            default_flow_style=False,
        )

# create shell script to run workflow
if not path_shell_script.exists():
    log.info(f"Creating shell script {path_shell_script}")
    cmd = f"""#!/bin/bash
    snakemake \\
        --verbose \\
        --configfile={str(path_config_file)} \\
        --max-jobs-per-second=5 \\
        --snakefile={str(path_snakefile)} \\
        --use-conda \\
        --conda-frontend=mamba \\
        --jobs=8 \\
        --rerun-incomplete"""
    # save cmd to path_shell_script
    with open(path_shell_script, "w") as f:
        f.write(cmd)

log.info("Done")

# %%
