# this script gte a bed file of mononucleotide stretches of a given length on a given reference genome
# and intersects it with a signals bed-like file.
# stats are then computed to determine if a sequencing sample is impacted by mononucleotide stretches and if so,
# the signal overlapping the mononucleotides is filtered.

# %%

import argparse
import subprocess
import tempfile
from pathlib import Path
from shlex import split

from logzero import logging as log

from ..signalprocessing import alignments_to_rafs
from ..util import util

# %%

# def yield_signaldepths_per_chr(chrID:int,filename:Path) -> typing.Generator[list, None, None]:
#     cmd_tabix = f"tabix -0 -f {str(filename)} {str(chrID)}"
#     process = subprocess.Popen(shlex.split(cmd_tabix), stdout=subprocess.PIPE)
#     for line in process.stdout:
#         raw = line.decode("utf-8").strip().split('\t')
#         yield raw


# def process_chr(mononucleotides_tree:IntervalTree,chrID:int,path_signaldepths:Path,output:Path,margin:int) -> None:
#     with open(output, 'w') as outf:
#         writer = csv.writer(outf, delimiter='\t')
#         for line in yield_signaldepths_per_chr(chrID=chrID,filename=path_signaldepths):
#             if len(line) < 7:
#                 log.warning(f"Skipping line: {line} because {len(line)} < 7")
#                 continue
#             sv_type = int(line[3])
#             sv_size = int(line[6])
#             start = max(0,int(line[1])-margin)
#             end = int(line[2])+margin
#             if sv_type != 1 or sv_size >= 100:
#                 writer.writerow(line)
#                 continue
#             overlapping_mononucleotides:set[Interval] = mononucleotides_tree.overlap(start,end)
#             # check for each mononucleotide if the signal is larger than 2x the mononucleotide
#             if len(overlapping_mononucleotides) == 0:
#                 writer.writerow(line)
#                 continue
#             for mononucleotide in overlapping_mononucleotides:
#                 if sv_size > 2*mononucleotide.data:
#                     writer.writerow(line)
#                     break


# def mp_process_function(args:dict) -> None:
#     return process_chr(**args)


# def lines_to_intervaltree(lines:list[list]) -> IntervalTree:
#     tree = IntervalTree()
#     for line in lines:
#         tree.addi(int(line[1]),int(line[2]),int(line[2])-int(line[1]))
#     return tree


# def mp_process_lines_to_intervaltree(args:dict) -> tuple[int,IntervalTree]:
#     return args['chrID'],lines_to_intervaltree(args['lines'])


# def filter_mononucleotide_exclusive_deletions(
#         mononucleotides:Path,
#         signaldepths:Path,
#         reference:Path,
#         output:Path,
#         margin:int,
#         threads:int,
#         tmp_dir_path:Path|None=None) -> None:
#     # build intervaltrees from mononucleotides (bed file)
#     # get chrIDs and chr Names

#     log.info("Building interval trees from mononucleotides")
#     ref_dict = util.create_ref_dict(reference)
#     ref_chr_to_ID_dict = {v:k for k,v in ref_dict.items()}
#     chrTrees = {i:IntervalTree() for i,_ in ref_dict.items()}
#     # load all mononucleotides lines
#     dict_lines_mononucleotides = {chrID:[] for chrID in chrTrees.keys()}
#     for line in csv.reader(open(mononucleotides, 'r'), delimiter='\t'):
#         chrID = ref_chr_to_ID_dict[line[0]]
#         dict_lines_mononucleotides[chrID].append(line)

#     jobs_args = [
#         {'lines':dict_lines_mononucleotides[chrID],
#          'chrID':chrID} for chrID in dict_lines_mononucleotides.keys()]

#     if threads > 1:
#         with mp.Pool(threads) as pool:
#             chrTrees = dict(pool.map(mp_process_lines_to_intervaltree, jobs_args, chunksize=1))
#             #chrTrees = {chrID:tree for chrID,tree in pool.map(mp_process_lines_to_intervaltree, jobs_args, chunksize=1)}
#     else:
#         chrTrees = {args['chrID']:lines_to_intervaltree(args['lines']) for args in jobs_args}

#     # create for each chrID a tmp output file
#     tmp_outputs = {i:tempfile.NamedTemporaryFile(dir=tmp_dir_path, prefix=f"tmp_out.{i}.",suffix=".tsv", delete=False if tmp_dir_path else True)
#                    for i in chrTrees.keys()}
#     log.info("creating jobs for multiprocessing")
#     jobs = [
#         {
#             'mononucleotides_tree':chrTrees[chrID],
#             'chrID':chrID,
#             'path_signaldepths':signaldepths,
#             'output':Path(tmp_outputs[chrID].name),
#             'margin':margin
#         }
#         for chrID in chrTrees.keys()]

#     if threads > 1:
#         with mp.Pool(threads) as pool:
#             pool.map(mp_process_function, jobs, chunksize=1)
#     else:
#         for job in jobs:
#             mp_process_function(job)

#     # cat all files together
#     log.info(f"Concatenating all files: {' '.join([outf.name for outf in tmp_outputs.values()])}")
#     tmp_output = tempfile.NamedTemporaryFile(prefix="tmp_all_output.",suffix=".tsv",delete=False if tmp_dir_path else True)
#     util.concat_txt_files(input_files=[v.name for v in tmp_outputs.values()],output=Path(tmp_output.name))

#     log.info(f"sorting, compressing and indexing {tmp_output.name} to {output}")
#     alignments_to_rafs.compress_and_index_bedlike(
#         sort_numerically=True,
#         input=tmp_output.name,
#         output=output,
#         threads=threads)
#     log.info("done")


def filter_mononucleotide_exclusive_deletions(
    mononucleotides: Path,
    signaldepths: Path,
    reference: Path,
    output: Path,
    margin: int,
    threads: int,
) -> None:
    # add margin of 2 bp to mononucleotides
    log.info(f"Adding margin of {margin} bp to mononucleotides")
    tmp_slop_mononucleotides = tempfile.NamedTemporaryFile(
        delete=True, suffix="slop_mononucleotides.bed"
    )
    cmd_slop = split(
        f"bedtools slop -b {str(margin)} -i {str(mononucleotides)} -g {str(reference)}.fai"
    )
    with open(tmp_slop_mononucleotides.name, "w") as f:
        subprocess.run(cmd_slop, stdout=f)

    # convert mononucleotides chr to chrID
    log.info("Converting mononucleotides chr to chrID")
    tmp_chrID_mononucleotides = tempfile.NamedTemporaryFile(
        delete=True, suffix="chrID_mononucleotides.bed"
    )
    util.bed_chr_to_chrID(
        input=mononucleotides,
        output=tmp_chrID_mononucleotides.name,
        reference=reference,
    )

    # intersect mononucleotides with signals so that a selection of signals is left
    tmp_out = tempfile.NamedTemporaryFile(delete=True, suffix=".out.signaldepths.tsv")
    log.info("Intersecting mononucleotides with signals")
    cmd_intersect = split(
        f"bedtools intersect -v -f 0.5 -r -a {str(signaldepths)} -b {tmp_chrID_mononucleotides.name}"
    )
    with open(tmp_out.name, "w") as f:
        subprocess.run(cmd_intersect, stdout=f)

    alignments_to_rafs.compress_and_index_bedlike(
        sort_numerically=True, input=tmp_out.name, output=output, threads=threads
    )


# %%


def run(args, **kwargs):
    filter_mononucleotide_exclusive_deletions(
        mononucleotides=args.mononucleotides,
        signaldepths=args.signaldepths,
        reference=args.reference,
        output=args.output,
        margin=args.margin,
        threads=args.threads,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        description="Filters signals for deletions that are mostly covered by mononucleotide stretches."
    )
    parser.add_argument(
        "-m",
        "--mononucleotides",
        type=Path,
        required=True,
        help="Path to mononucleotides bed file.",
    )
    parser.add_argument(
        "-s",
        "--signaldepths",
        type=Path,
        required=True,
        help="Path to bgzipped signaldepths tsv file.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Path to reference fasta file.",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to output file."
    )
    parser.add_argument(
        "-g",
        "--margin",
        type=int,
        required=False,
        default=5,
        help="Margin to add to mononucleotides.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use for compression.",
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    run(args)
    return


if __name__ == "__main__":
    main()
