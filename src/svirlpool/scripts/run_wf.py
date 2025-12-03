# script to run the full workflow for a single sample

# define a function process_single_sample that is called from main
# it receives all input that is necessary to build a wf_config file

# define a config class that is parsed to a json as the example bewlow
# %%
import json
import subprocess
from pathlib import Path
from shlex import split

# %%


def check_annotations_file(
    input: Path,
    reference_fai: Path,
) -> None:
    if not input.exists() or input.stat().st_size == 0:
        raise FileNotFoundError(
            f"Annotation file {str(input)} does not exist or is empty."
        )
    # read chromosomes from reference_fai
    with open(reference_fai, "r") as f:
        reference_chroms = set()
        for line in f:
            reference_chroms.add(line.split("\t")[0])
    # read chromosomes from input annotation file
    with open(input, "r") as f:
        input_chroms = set()
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            input_chroms.add(line.split("\t")[0])
    # check if all chromosomes in input_chroms are in reference_chroms
    if not input_chroms.issubset(reference_chroms):
        missing_chroms = input_chroms - reference_chroms
        raise ValueError(
            f"""
Annotation file {str(input)} contains chromosomes not present in the reference fasta index {str(reference_fai)}.
Missing chromosomes: {", ".join(missing_chroms)}
To filter the annotation file, you can use the following bash command:
bedtools intersect -a {str(input)} -b <(awk '{{print $1"\t0\t"$2}}' {str(reference_fai)}) > {str(input).replace(".bed", ".filtered.bed")}
            """
        )


def validate_input(args) -> None:
    reference_fai = Path(str(args.reference) + ".fai")
    if not reference_fai.exists() or reference_fai.stat().st_size == 0:
        raise FileNotFoundError(
            f"Reference fasta index {str(reference_fai)} does not exist or is empty."
        )
    # check annotations files
    check_annotations_file(input=Path(args.trf), reference_fai=reference_fai)
    check_annotations_file(
        input=Path(args.mononucleotides), reference_fai=reference_fai
    )


def config_create_json(config: dict, path_json: Path):
    with open(path_json, "w") as f:
        json.dump(config, f)


def run_wf(args):
    validate_input(args)
    if not Path(args.workdir).exists():
        Path(args.workdir).mkdir(parents=True)
    dict_args = vars(args)
    # remove 'fast'and 'func' keys from dict_args
    dict_args.pop("func")
    dict_args.pop("fast")
    config_create_json(config=dict_args, path_json=Path(args.workdir) / "config.json")
    path_base_wf = Path(__file__).parent.parent / "workflows/main.smk"
    cmd_wf = split(
        f"snakemake \
        {'--unlock' if args.snakemake_unlock else ''} \
        {'--executor slurm --jobs ' + str(args.executor_slurm_jobs) if args.executor_slurm_jobs > 0 else ''} \
        --rerun-incomplete \
        --snakefile={str(path_base_wf)} \
        --configfile={Path(args.workdir) / 'config.json'} \
        --max-jobs-per-second=5 \
        --cores {args.threads}"
    )
    subprocess.check_call(cmd_wf)
