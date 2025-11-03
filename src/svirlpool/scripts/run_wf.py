# script to run the full workflow for a single sample

# define a function process_single_sample that is called from main
# it receives all input that is necessary to build a wf_config file

# define a config class that is parsed to a json as the example bewlow
# %%
import json
import subprocess
from pathlib import Path, PosixPath
from shlex import split

# %%


def config_create_json(config: dict, path_json: PosixPath):
    with open(path_json, "w") as f:
        json.dump(config, f)


def run_wf(args):
    # if len(args.alignments) != 1:
    #     raise ValueError("alignments can only accept one file right now.")
    # if workdir does not exist, create it
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
        {'--executor slurm --jobs '+str(args.executor_slurm_jobs) if args.executor_slurm_jobs > 0 else ''} \
        --rerun-incomplete \
        --snakefile={str(path_base_wf)} \
        --configfile={Path(args.workdir) / 'config.json'} \
        --max-jobs-per-second=5 \
        --cores {args.threads}"
    )
    subprocess.check_call(cmd_wf)
