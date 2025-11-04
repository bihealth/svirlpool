# --- SV calling module --- #
# queries a database with tables
# candidate_regions, consensuses, unused_reads, consensusAlignments__{name_reference}
# loads all candidate regions, and all alignments to dictionaries
# iterates the consensus objects
# - id VARCHAR(100) PRIMARY KEY
# - consensus BLOB

# loads all MergedSVSignals that correspond to name_reference and the consensusID
# takes from each MergedSVSignal the signal and generates a sv candidate from it.
# traces back the positon of the signal in the consensus via the cut read alignments
# and genotypes the SV candidate, subtracting the unused reads
# %%

import csv
import multiprocessing as mp
import pickle
import sqlite3
import subprocess
import tempfile
from copy import deepcopy
from datetime import datetime
from math import ceil
from pathlib import Path
from shlex import split
from typing import Generator

import attrs
import cattrs
import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree
from logzero import logger as log
from tqdm import tqdm

from . import SVprimitive_class, consensus_class, datatypes, genotyping, util

# %%


def find_closest_depth_estimator(
    depth_of_coverages: list[tuple[str, int, int]], position: int
) -> tuple[str, int, int]:
    if len(depth_of_coverages) == 0:
        raise ValueError("depth_of_coverages must not be empty")
    if len(depth_of_coverages) == 1:
        return depth_of_coverages[0]
    return min(depth_of_coverages, key=lambda x: abs(x[1] - position))


def print_it_cutreads_to_consensus_to_terminal(
    consensus: consensus_class.Consensus,
    it_cutreads_to_consensus: IntervalTree,
    line_width: int = 60,
) -> float:
    # calculate a scale factor for the consensus to fit its length to the terminal width
    scale_factor = line_width / len(consensus.consensus_sequence)
    # print a concatenation of '=' to represent the consensus sequence followed by a white space and the consensusID
    print("=" * line_width, consensus.ID, sep=" ")
    # then iterate the interval tree and create a new line to print the intervals of the cut reads
    # after position line_width+1 comes the sampleID (data of the interval)
    for interval in it_cutreads_to_consensus:
        s, e = interval.begin, interval.end
        # scale s,e
        s, e = int(ceil(s * scale_factor)), int(ceil(e * scale_factor))
        line = list("." * line_width)
        line[s:e] = ["|"] * (e - s)
        l = "".join(line)
        l += " " + str(interval.data)
        print(l)
    return scale_factor


def print_svPrimitives_to_terminal(
    svPrimitives: list[SVprimitive_class.SVprimitive],
    scale_factor_terminal: float,
    line_width: int = 60,
):
    # print each SVprimitive in a line where it starts on the consensus, followed by its type and size
    for svPrimitive in svPrimitives:
        start = int(ceil(svPrimitive.read_start * scale_factor_terminal))
        end = int(ceil(svPrimitive.read_end * scale_factor_terminal))
        line = list("." * line_width)
        len_event = max(1, end - start)
        line[start : start + len_event] = ["X"] * len_event
        l = "".join(line)
        l += " " + str(svPrimitive.sv_type) + " " + str(svPrimitive.size)
        print(l)


def print_reconstructed_alignments_to_terminal(
    consensus: consensus_class.Consensus, line_width: int = 60
):
    alignments = [rs.alignment.to_pysam() for rs in consensus.reconstructible_reads]
    util.display_ascii_alignments(
        alignments=alignments,
        max_ref_len=len(consensus.consensus_sequence),
        terminal_width=line_width,
    )


def svPrimitives_from_consensus(
    samplename: str,
    consensus: consensus_class.Consensus,
    consensus_alignments: list[consensus_class.ConsensusAlignment],
    min_signal_size: int = 30,
    verbose: bool = False,
) -> list[SVprimitive_class.SVprimitive]:
    results: list[SVprimitive_class.SVprimitive] = []
    depth_of_coverages = consensus.depth_of_coverages
    if len(depth_of_coverages) == 0:
        log.warning(f"depth_of_coverages is empty. Skipping consensus {consensus.ID}")
        return []
    intervals_cutread_alignments = consensus.intervals_cutread_alignments

    it_cutreads_to_consensus = IntervalTree()
    for interval in intervals_cutread_alignments:
        s, e = (
            (interval[0], interval[1])
            if interval[0] <= interval[1]
            else (interval[1], interval[0])
        )
        it_cutreads_to_consensus.addi(begin=s, end=e, data=interval[2])
    if verbose:
        print_reconstructed_alignments_to_terminal(consensus=consensus)
    # if proto_svs in consensus_alignment is empty, return empty list
    if all(
        [
            len(consensus_alignment.proto_svs) == 0
            for consensus_alignment in consensus_alignments
        ]
    ):
        return []

    # TODO: find and re-cluster and assemble all consensus sequences and their cut reads that have many proto_svs that are in none or different repeats
    # remove those from the set that is then passed to the alignment-based SV calling (below)
    # re-clustering based sv-calling
    # (...)
    # SV calling of other types of SVs
    # - deletions - connected BNDs without separating sequence
    # - duplication - inserted sequence that aligns to the consensus sequence
    # - inversion - three alignments with the middle part being reversed

    # alignment-based SV calling
    for idx_aln, consensus_alignment in enumerate(consensus_alignments):
        if len(consensus_alignment.proto_svs) == 0:
            continue

        # check if each consensus_alignment has at least one mergedSVSignal in proto_svs
        for idx_sv, mergedSVSignal in enumerate(consensus_alignment.proto_svs):
            if (
                mergedSVSignal.sv_type <= 2 and mergedSVSignal.size < min_signal_size
            ):  # filter very small indels
                continue
            # for debugging, print alt_sequence and ref_sequence of the mergedSVSignal
            # log.debug(f"alt_sequence: {mergedSVSignal.alt_sequence}, ref_sequence: {mergedSVSignal.ref_sequence}")
            # get all readnames that overlap read_start and read_end of the mergedSVSignal
            # if distance mergedSVSignal.read_start and mergedSVSignal.read_end is less than 30, only use start
            assert (
                mergedSVSignal.read_start <= mergedSVSignal.read_end
            ), "read_start must be smaller than or equal read_end"
            assert (
                len(mergedSVSignal.original_signals) > 0
            ), "original_signals must not be empty"
            assert (
                len(consensus.original_signals) > 0
            ), "original SV signals from consensus must not be empty"
            assert (
                type(consensus_alignment.reference_name) == str
            ), "reference_name must be a string"

            new_sv_primitive = SVprimitive_class.SVprimitive(
                **attrs.asdict(mergedSVSignal),
                original_cr_signals=[
                    signal.unstructure() for signal in consensus.original_signals
                ],
                consensusID=consensus.ID,
                alignment_to_ref=consensus_alignment.alignment,  # necessary to write down BNDs in the vcf correctly
                genotypeMeasurements=dict(),
                reference_name=consensus_alignment.reference_name,
            )
            new_sv_primitive.vcfID = f"{datatypes.SV_TYPE_DICT[new_sv_primitive.sv_type]}_{samplename}_{consensus.ID}-{str(idx_aln)}.{str(idx_sv)}"
            if verbose:
                print_svPrimitives_to_terminal(
                    [new_sv_primitive], 60 / len(consensus.consensus_sequence)
                )
            # only one genotype measurement per sampleID is generated
            overlapping_readnames_start = set(
                [
                    interval.data
                    for interval in it_cutreads_to_consensus[mergedSVSignal.read_start]
                ]
            )

            closest_depth_estimator: tuple[str, int, int] = (
                find_closest_depth_estimator(
                    depth_of_coverages=depth_of_coverages,
                    position=mergedSVSignal.read_start,
                )
            )
            new_sv_primitive.genotypeMeasurements[samplename] = (
                genotyping.GenotypeMeasurement(
                    samplename=samplename,
                    position_on_consensus=mergedSVSignal.read_start,
                    estimated_total_depth=closest_depth_estimator[2],
                    distance_to_closest_depth_estimator=abs(
                        closest_depth_estimator[1] - mergedSVSignal.read_start
                    ),
                    supporting_reads=list(overlapping_readnames_start),
                )
            )
            if mergedSVSignal.read_end - mergedSVSignal.read_start >= min_signal_size:
                # dict_closest_depth_estimators = find_closest_depth_estimators(depth_of_coverages=depth_of_coverages, position=mergedSVSignal.read_start)
                # add another genotype measurement for each sampleID for read_end
                overlapping_readnames_end = set(
                    [
                        interval.data
                        for interval in it_cutreads_to_consensus[
                            mergedSVSignal.read_end
                        ]
                    ]
                )
                # generate dict from readnames_start sampleID:readname (sampleID is the suffix after '.' in the readname)
                # dict_closest_depth_estimators = find_closest_depth_estimators(depth_of_coverages=depth_of_coverages, position=mergedSVSignal.read_end)
                new_genotype_measurement = genotyping.GenotypeMeasurement(
                    samplename=samplename,
                    position_on_consensus=mergedSVSignal.read_end,
                    estimated_total_depth=closest_depth_estimator[2],
                    distance_to_closest_depth_estimator=abs(
                        closest_depth_estimator[1] - mergedSVSignal.read_end
                    ),
                    supporting_reads=list(overlapping_readnames_end),
                )

                if len(
                    new_sv_primitive.genotypeMeasurements[samplename].supporting_reads
                ) < len(new_genotype_measurement.supporting_reads):
                    new_sv_primitive.genotypeMeasurements[samplename] = (
                        new_genotype_measurement
                    )

            # check if all sampleIDs have a GenotypeMeasurement. If not, raise an error
            # new_sv_primitive now should have at least one GenotypeMeasurement. check if that exists. If not, raise an error
            if samplename not in new_sv_primitive.genotypeMeasurements:
                raise ValueError(
                    f"sample {samplename} has no GenotypeMeasurement in new_sv_primitive"
                )

            # since we have the cut_read_alignment_signals on the consensus object, we can estimate a distortion
            # parameter that is given to the resulting sv_primitive. To allow the calculation of the distortion
            # at a lter point, here we just collect the necessary information from the cut_read_alignment_signals
            if len(consensus.reconstructible_reads) < len(
                consensus.cut_read_alignment_signals
            ):
                raise ValueError(
                    "len(consensus.reconstructible_reads) < len(consensus.cut_read_alignment_signals)"
                )
            sums_alignments = [0.0] * len(consensus.reconstructible_reads)
            for i, readAlignmentSignals in enumerate(
                consensus.cut_read_alignment_signals
            ):
                if (
                    readAlignmentSignals.read_name
                    not in new_sv_primitive.distortion_sum_per_read
                ):
                    new_sv_primitive.distortion_sum_per_read[
                        readAlignmentSignals.read_name
                    ] = 0.0
                sum_signals = 0
                for signal in readAlignmentSignals.SV_signals:
                    if signal.sv_type <= 1:  # insertion or deletion
                        sum_signals += signal.size
                    else:
                        sum_signals += (
                            0  # break ends don't contribute to the distortion
                        )
                new_sv_primitive.distortion_sum_per_read[
                    readAlignmentSignals.read_name
                ] += sum_signals
                sums_alignments[i] += sum_signals
            new_sv_primitive.distortion = sums_alignments
            new_sv_primitive.aln_is_reverse = (
                consensus_alignment.alignment.to_pysam().is_reverse
            )
            # add original signals to SVprimitive
            results.append(new_sv_primitive)
    # cluster the signals from ReadAlignmentSignals from Consensus with each sv_primitive from results as their centroids

    # check all results for completeness. each SVprimitive in results should have at least one GenotypeMeasurement
    svPrimitives_without_genotypeMeasurements = [
        svPrimitive
        for svPrimitive in results
        if len(svPrimitive.genotypeMeasurements) == 0
    ]
    if any(svPrimitives_without_genotypeMeasurements):
        raise ValueError(
            f"{len(svPrimitives_without_genotypeMeasurements)} SVprimitives without GenotypeMeasurements of a total of {len(results)}"
        )
    # else:
    #     log.info(f"All SVprimitives have GenotypeMeasurements of a total of {len(results)}")
    # write to database

    # break ends would occurr here! Add adjacencies..
    results = add_adjacencies_to_svPrimitives_breakends(results)
    return results


def add_adjacencies_to_svPrimitives_breakends(
    svPrimitives: list[SVprimitive_class.SVprimitive],
) -> list[SVprimitive_class.SVprimitive]:
    # sort all svPrimitives by consensusID, chr, and read_start
    svPrimitives = sorted(svPrimitives, key=lambda x: (x.consensusID, x.read_start))

    # two BNDs are adjacent if they are consecutive in the list and have the same consensusID
    # given a pair (A,B) of adjacent svPs:
    i = 0
    while i < len(svPrimitives) - 1:
        A: SVprimitive_class.SVprimitive = svPrimitives[i]
        B: SVprimitive_class.SVprimitive = svPrimitives[i + 1]
        # print(f"Checking {A.vcfID} and {B.vcfID} for adjacency")
        # print(f"A.consensusID == B.consensusID: {A.consensusID == B.consensusID} and A.alignment_to_ref != B.alignment_to_ref: {A.alignment_to_ref != B.alignment_to_ref} and A.sv_type >= 3 and B.sv_type >= 3: {A.sv_type >= 3 and B.sv_type >= 3}")
        if (
            A.consensusID == B.consensusID
            and A.alignment_to_ref != B.alignment_to_ref
            and A.sv_type >= 3
            and B.sv_type >= 3
        ):
            tp_str_A = ""
            tp_str_B = ""
            # if A and B are BNDs, add the adjacency to both
            if A.sv_type == 4:  # BNDR
                if B.sv_type == 3:  # t[p[, ]p]t
                    tp_str_A = "t[p["
                    tp_str_B = "]p]t"
                if B.sv_type == 4:  # t]p], t]p]
                    tp_str_A = "t]p]"
                    tp_str_B = "t]p]"
            if A.sv_type == 3:  # BNDL
                if B.sv_type == 3:  # [p[t, [p[t
                    tp_str_A = "[p[t"
                    tp_str_B = "[p[t"
                if B.sv_type == 4:  # ]p]t, t[p[
                    tp_str_A = "]p]t"
                    tp_str_B = "t[p["
            A.adjacent_bnd = datatypes.Adjacency(
                chr=B.chr,
                ref_pos=B.ref_start,
                gapsize=abs(B.read_start - A.read_start - 1),
                tp_str=tp_str_A,
                mate_id=B.vcfID,
            )
            B.adjacent_bnd = datatypes.Adjacency(
                chr=A.chr,
                ref_pos=A.ref_start,
                gapsize=abs(B.read_start - A.read_start - 1),
                tp_str=tp_str_B,
                mate_id=A.vcfID,
            )
            i += 2
        else:
            i += 1

    return svPrimitives


def generate_header(reference: Path, samplenames: list[str]) -> list[str]:
    reference = Path(reference)
    header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now().strftime('%Y%m%d')}",
        f"##reference=file://{str(reference.absolute())}",
    ]
    # add contigs to header
    # load reference index as pd dataframe
    ref_index = pd.read_csv(
        reference.with_suffix(reference.suffix + ".fai"), sep="\t", header=None
    )
    # compute lengths of contigs
    # add contigs to header
    for i in range(ref_index.shape[0]):
        header.append(
            f"##contig=<ID={str(ref_index.iloc[i,0])},length={ref_index.iloc[i,1]}>"
        )
    header.append(
        '##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of informative reads.">'
    )
    header.append(
        '##FILTER=<ID=PASS,Description="high quality and sufficient number of informative reads.">'
    )
    header.append(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">'
    )
    header.append(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion, DUP=Duplication, INV=Inversion">'
    )
    header.append(
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    header.append(
        '##INFO=<ID=PASS_ALTREADS,Number=1,Type=String,Description="Passed alt reads threshold">'
    )
    header.append(
        '##INFO=<ID=pass_GQ,Number=1,Type=String,Description="Passed Genotype precision threshold">'
    )
    header.append(
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">'
    )
    header.append(
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">'
    )
    header.append(
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">'
    )
    header.append(
        '##INFO=<ID=CONSENSUSID,Number=1,Type=String,Description="ID of the consensus that this SV originiates from. Other consensus sequences can also be involved.">'
    )
    header.append('##ALT=<ID=INS,Description="Insertion">')
    header.append('##ALT=<ID=DEL,Description="Deletion">')
    header.append('##ALT=<ID=DUP,Description="Duplication">')
    header.append('##ALT=<ID=INV,Description="Inversion">')
    header.append('##ALT=<ID=BND,Description="Breakend; Translocation">')
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append(
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
    )
    header.append('##FORMAT=<ID=TC,Number=1,Type=Integer,Description="Total coverage">')
    header.append(
        '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">'
    )
    header.append(
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">'
    )
    header.append(
        '##FORMAT=<ID=GP,Number=1,Type=Float,Description="Genotype probability">'
    )

    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samplenames)
    )
    return header


def svPrimitives_to_vcf(
    svPrimitives: list[SVprimitive_class.SVprimitive],
    reference: Path,
    output: Path,
    min_support: int = 3,
    min_sv_size: int = 30,
    sv_types: list[str] = ["DEL", "INS"],
    tmp_dir_path: Path | None = None,
):
    allowed_svs = {0: "INS", 1: "DEL", 3: "BND", 4: "BND"}
    with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tdir:
        # assert sv_tpes. it can only be DEL or INS or both
        assert all(
            [sv_type in allowed_svs.values() for sv_type in sv_types]
        ), f"sv_types must be {allowed_svs.values()}"
        # generate header
        samplenames = list(
            set(
                [
                    gM.samplename
                    for svp in svPrimitives
                    for gM in svp.genotypeMeasurements.values()
                ]
            )
        )
        header = generate_header(reference=reference, samplenames=samplenames)

        tmp_unsorted_output = tempfile.NamedTemporaryFile(
            dir=tdir, mode="w", delete=True, suffix=".vcf"
        )
        with open(str(tmp_unsorted_output.name), "w") as out:
            for line in header:
                print(line, file=out)
            log.info(f"Writing {len(svPrimitives)} SVprimitives to vcf file...")
            for svPrimitive in svPrimitives:
                if allowed_svs[svPrimitive.sv_type] not in sv_types:
                    continue
                if svPrimitive.size < min_sv_size:
                    continue
                # find start of reference intervals
                # pick minimum start of all intervals
                vcf_line = sv_primitive_to_vcf_line(
                    samplenames, min_support, svPrimitive
                )
                print(*vcf_line, sep="\t", file=out)
        log.info(f"Created vcf file at {str(output)}. zippping and indexing...")
        # sort vcf file
        cmd_sort = f"bcftools sort {str(tmp_unsorted_output.name)} -o {str(output)}"
        subprocess.check_call(split(cmd_sort))
        with open(str(output) + ".gz", "wb") as f:
            cmd_zip = f"bgzip -c {str(output)}"
            subprocess.check_call(split(cmd_zip), stdout=f)
        cmd_index = f"tabix -f -0 -p vcf {str(output)}.gz"
        subprocess.check_call(split(cmd_index))


from scipy.stats import binom


def genotype_likelihood(
    n_alt_reads: int, n_total_reads: int, genotype_probs: dict
) -> tuple[str, dict[str, float]]:
    if n_total_reads == 0:
        return "0/0", {"0/0": 1.0}
    # Compute binomial probabilities
    likelihoods = {
        genotype: binom.pmf(n_alt_reads, n_total_reads, p)
        for genotype, p in genotype_probs.items()
    }
    # Normalize probabilities (optional, converts likelihoods into posterior probabilities)
    total_likelihood = sum(likelihoods.values())
    probabilities = {
        genotype: likelihood / total_likelihood
        for genotype, likelihood in likelihoods.items()
    }
    # Determine the most likely genotype
    best_genotype: str = max(probabilities, key=probabilities.get)
    # replace any nan values with 0.0
    for genotype, probability in probabilities.items():
        if np.isnan(probability):
            probabilities[genotype] = 0.0
    return best_genotype, probabilities


def sv_primitive_to_vcf_line(
    samplenames: list[str], min_support: int, svPrimitive: SVprimitive_class.SVprimitive
) -> list[str]:
    ONE_BASED: int = (
        1  # up to this point, all genomic coordinates are zero-based, but the vcf format is one-based.
    )
    # "assumes that the svPrimitive has the dynamic attribute 'initial_covs', which is a dict samplename:(cov_start, cov_end) of typedict[str,tuple[int,int]]"
    # check input
    # if not hasattr(svPrimitive,'initial_covs'):
    #     raise ValueError(f"svPrimitive {svPrimitive.vcfID} has no initial_covs attribute. They are however necessary to write the vcf line.")
    # svPrimitive needs to have a consensusID (str) longer than 0
    assert len(svPrimitive.consensusID) > 0, "svPrimitive must have a consensusID"
    # svPrimitive needs to have a alignment_to_ref that is not None
    # assert svPrimitive.alignment_to_ref is not None, "svPrimitive must have an alignment_to_ref"
    # assert isinstance(svPrimitive.alignment_to_ref,datatypes.Alignment), "alignment_to_ref must be an instance of Alignment"
    # svPrimitive needs to have a genotypeMeasurements dict, that has ints as keys and GenotypeMeasurement objects as values
    assert (
        len(svPrimitive.genotypeMeasurements) > 0
    ), "svPrimitive must have genotypeMeasurements"
    for samplename, genotype_measurement in svPrimitive.genotypeMeasurements.items():
        assert isinstance(
            genotype_measurement, genotyping.GenotypeMeasurement
        ), "genotype_measurement must be an instance of GenotypeMeasurement"
    # svPrimitive needs to have a list original_cr_signals, with more than 0 elemens and all of type ExtendedSVsignal
    assert (
        len(svPrimitive.original_cr_signals) > 0
    ), "svPrimitive must have original_cr_signals"

    FORMAT_field = "GT:GQ:TC:DR:DV:GP"
    start = svPrimitive.ref_start
    end = (
        (svPrimitive.ref_start + svPrimitive.size)
        if svPrimitive.sv_type == 1
        else (svPrimitive.ref_start + 1)
    )
    chr = svPrimitive.chr
    GTs, GQs, DRs, DVs, TCs, GPs = [], [], [], [], [], []
    dict_genotypes = {samplename: [] for samplename in samplenames}
    # dict_genotypes is filles with all SVprimitives that have the same samplename
    for samplename, genotype_measurement in svPrimitive.genotypeMeasurements.items():
        dict_genotypes[genotype_measurement.samplename].append(genotype_measurement)
        # now iterate over the samplenames in dict_genotypes
        # if the list is empty, add (./.,0,0,0) to the FORMAT field
        # else, merge the list of genotype_measurements to fill the FORMAT field

    for samplename, genotype_measurements in dict_genotypes.items():
        if len(genotype_measurements) == 0:
            if svPrimitive.initial_covs:
                # take measurement from svPrimitive.initial_covs to get the coverage per samplename
                # here is the question what metric is useful. is it max, min or mean? defensive would be min.
                try:
                    initial_cov = min(svPrimitive.initial_covs[samplename])
                except KeyError:
                    import pdb

                    pdb.set_trace()
                GTs.append("0/0" if initial_cov > 0 else "./.")
                # GQs.append(initial_cov / max(svPrimitive.initial_covs[samplename]) if max(svPrimitive.initial_covs[samplename]) > 0 else 0)
                GQs.append(0)
                DRs.append(initial_cov)
                DVs.append(0)
                TCs.append(initial_cov)
                GPs.append(1.0)
                continue
            else:
                GTs.append("./.")
                GQs.append(0)
                DRs.append(0)
                DVs.append(0)
                TCs.append(0)
                GPs.append(0.0)
        else:
            # at first, estimate the support. take the mean numer of supporting_reads from all genotype_measurements
            mean_support = int(
                round(
                    np.mean(
                        [
                            len(genotype_measurement.supporting_reads)
                            for genotype_measurement in genotype_measurements
                        ]
                    )
                )
            )
            mean_depth = int(
                round(
                    np.mean(
                        [
                            genotype_measurement.estimated_total_depth
                            for genotype_measurement in genotype_measurements
                        ]
                    )
                )
            )
            # mean depth estimator might be off. correct it:
            mean_depth = max(mean_depth, mean_support)
            mean_neutral = int(
                round(
                    np.mean(
                        [
                            len(genotype_measurement.neutral_reads)
                            for genotype_measurement in genotype_measurements
                        ]
                    )
                )
            )  # not used right now
            best_genotype, gt_probabilities = genotype_likelihood(
                n_alt_reads=mean_support,
                n_total_reads=mean_depth,
                genotype_probs={"0/0": 0.05, "0/1": 0.5, "1/1": 0.95},
            )
            GTs.append(best_genotype)
            TCs.append(mean_depth)
            # GQs.append(mean_depth-mean_neutral)
            GQs.append(
                int(-10 * np.log10(1.0 - gt_probabilities[best_genotype]))
                if gt_probabilities[best_genotype] < 1.0
                else 60
            )
            DRs.append(mean_depth - mean_neutral - mean_support)
            DVs.append(mean_support)
            GPs.append(gt_probabilities[best_genotype])

    format_content = list(
        map(
            lambda x: f"{x[0]}:{x[1]}:{x[2]}:{x[3]}:{x[4]}:{x[5]}",
            zip(GTs, GQs, TCs, DRs, DVs, GPs),
        )
    )
    # filter section - check if sv passes filter
    # pass_alt_reads = sum(DVs) >= min_support and sum(DVs)/(sum(DRs)+sum(DVs)) >= 0.15 if sum(DRs)+sum(DVs) else 0
    pass_alt_reads = (
        sum(DVs) > min_support
    )  # int(sum(DVs) / sum(TCs) >= 0.15 if sum(TCs) else 0)
    pass_GQ = int(
        any([GQ >= min_support for GQ in GQs])
    )  # TODO: add use of the depth estimator distance here!
    pass_filter = int(pass_alt_reads and pass_GQ)
    alt = svPrimitive.get_alt_sequence() if svPrimitive.sv_type != 1 else "."
    if svPrimitive.sv_type == 1:  # deletion -> add first letter of REF to alt
        alt = svPrimitive.get_ref_sequence()[0]
    # if svPrimitive.alignment_to_ref.to_pysam().is_reverse:
    # might not be necessary, since the alt sequence is already reverse complemented in the signal extraction and merging step in alignments_to_svprimitives.py
    # if svPrimitive.aln_is_reverse:
    #     alt = str(Seq(alt).reverse_complement())
    #
    ref = (
        svPrimitive.get_ref_sequence()
        if svPrimitive.sv_type == 1
        else svPrimitive.get_ref_sequence()[0]
    )  # only add ref sequence if it is a deletion
    if svPrimitive.sv_type == 0:  # insertion -> add first letter of REF to alt
        alt = ref[0] + alt

    if not svPrimitive.adjacent_bnd:
        if svPrimitive.sv_type == 3:
            alt = alt + "."
        elif svPrimitive.sv_type == 4:
            alt = "." + alt
    else:
        # define p and t to replace in adjacent_bnd.tp_str
        # TODO: test and correct. it's probably wrong!
        p = (
            svPrimitive.adjacent_bnd.chr
            + ":"
            + str(svPrimitive.adjacent_bnd.ref_pos + ONE_BASED)
        )
        t = alt[: svPrimitive.adjacent_bnd.gapsize]
        if svPrimitive.sv_type == 3:  # BND left
            t = ref + t
        if svPrimitive.sv_type == 4:  # BND right
            t = t + ref
        # t has the first base same as ref if
        alt = svPrimitive.adjacent_bnd.tp_str.replace("p", p).replace("t", t)

    if svPrimitive.sv_type >= 3 and len(svPrimitive.get_alt_sequence()) == 0:
        raise ValueError(f"svPrimitive {svPrimitive.vcfID} has no alt sequence")

    # TODO: PRECISE should really be attributed if the depth estimator is close AND the distortions of the respective consensus sequences are low.
    imprecise = (
        svPrimitive.original_alt_sequences
        and len(svPrimitive.original_alt_sequences) > 1
    ) or (
        svPrimitive.original_ref_sequences
        and len(svPrimitive.original_ref_sequences) > 1
    )
    imprecise = False  # for now, but it is necessary to look at genotypeMeasurements and maybe original_cr_signals to determine this
    label_imprecise = "IMPRECISE" if imprecise else "PRECISE"

    if imprecise:
        match svPrimitive.sv_type:
            case 0:
                alt = "<INS>"
            case 1:
                alt = "<DEL>"
            case _:
                alt = "<BND>"
    info_fields = {
        "PASS_ALTREADS": pass_alt_reads,
        "pass_GQ": pass_GQ,
        "SVTYPE": datatypes.SV_TYPE_DICT[svPrimitive.sv_type],
        "END": str(
            end + ONE_BASED - 1
        ),  # end is inclusive in vcf, so we need to subtract 1
        "SVLEN": str(svPrimitive.size * (-1 if svPrimitive.sv_type == 1 else 1)),
        "CONSENSUSID": svPrimitive.consensusID,
    }

    if svPrimitive.adjacent_bnd:
        info_fields["MATEID"] = svPrimitive.adjacent_bnd.mate_id
    info = (
        label_imprecise
        + ";"
        + ";".join([f"{key}={value}" for key, value in info_fields.items()])
    )

    vcf_line = [
        chr,
        start + ONE_BASED,
        svPrimitive.vcfID,  # ID
        ref,  # REF
        alt,  # ALT
        60,  # QUAL
        "PASS" if pass_filter else "LowQual",  # FILTER
        info,
        FORMAT_field,
        *format_content,
    ]

    # if svPrimitive.chr == '15':
    #     import pdb; pdb.set_trace()

    return vcf_line


def process_svPrimitives_for_consensus(
    kwargs: dict,
) -> list[SVprimitive_class.SVprimitive]:
    return svPrimitives_from_consensus(**kwargs)


def load_consensus_alignments_from_database(
    path_db: Path, name_reference: str, consensusIDs: list[str] | None = None
) -> dict[str, list[consensus_class.ConsensusAlignment]]:
    # iterate all consensus objects in database and construct a dict consensusID:consensusObject
    # log.info(f"loading consensus alignments from {path_db}...")
    if "." in name_reference:
        name_reference = name_reference.replace(".", "_")
        log.warning(
            f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {name_reference}"
        )
    if consensusIDs is not None:
        query = f"SELECT consensusID,consensusAlignments FROM consensusAlignments__{name_reference} WHERE consensusID IN ({','.join(['?']*len(consensusIDs))})"
        results = {consensusID: [] for consensusID in consensusIDs}
        conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
        c = conn.cursor()
        c.execute(query, consensusIDs)
    else:
        query = f"SELECT consensusID,consensusAlignments FROM consensusAlignments__{name_reference}"
        results = {}
        conn = sqlite3.connect("file:" + str(path_db) + "?mode=ro", uri=True)
        c = conn.cursor()
        c.execute(query)
    # result dict of the form consensusID:list[ConsensusAlignment]
    for row in c:
        results[row[0]] = cattrs.structure(
            pickle.loads(row[1]), list[consensus_class.ConsensusAlignment]
        )
    c.close()
    conn.close()
    if consensusIDs is None:
        consensusIDs = list(results.keys())
    # iterate all consensusIDs and log if there the value lists in results are empty
    for consensusID in consensusIDs:
        if not results[consensusID]:
            log.warning(
                f"consensusID {consensusID} was not found in the database {path_db}"
            )
    # report all lists in results are empty
    if all([not results[consensusID] for consensusID in consensusIDs]):
        log.warning(
            f"all consensusIDs {str(consensusIDs)} of this job were not found in the database {path_db}"
        )
    # check results if all consensusIDs are present
    return results


def generate_svPrimitives(
    samplename: str,
    name_reference: str,
    consensus_alignments: dict[str, list[consensus_class.ConsensusAlignment]],
    path_consensuses: Path,
    threads: int,
    chunksize: int,
    min_signal_size: int = 30,
    DB_timeout: float | None = None,
) -> list[SVprimitive_class.SVprimitive]:
    if DB_timeout is not None:
        SQLITE_TIMEOUT = DB_timeout

    threads = max(1, min(threads, mp.cpu_count()))

    log.info(
        "creating jobs to process consensus objects and their alignments to generate svPrimitives.."
    )
    jobs = [
        {
            "consensusID": consensusID,
            "samplename": samplename,
            "path_consensuses": path_consensuses,
        }
        for consensusID in consensus_alignments.keys()
    ]

    svPrimitives = []
    for i in tqdm(range(0, len(jobs), chunksize)):
        chunk = jobs[i : i + chunksize]  # list of dicts
        chunk_cids = [job["consensusID"] for job in chunk]
        consensus_objects = {
            c.ID: c
            for c in util.yield_consensus_objects(
                consensusIDs=chunk_cids, path_db=path_consensuses
            )
        }
        chunk_jobs = [
            {
                "samplename": samplename,
                "consensus": consensus_objects[cid],
                "consensus_alignments": consensus_alignments[cid],
                "min_signal_size": min_signal_size,
            }
            for cid in chunk_cids
        ]
        # process job_chunk in parallel with process_svPrimitives_for_consensus
        with mp.Pool(threads) as pool:
            # adjust to saving the list outputs of all jobs
            chunksize_mp = int(ceil(len(chunk_jobs) / threads))
            results = [
                svp
                for l in tqdm(
                    pool.imap_unordered(
                        func=process_svPrimitives_for_consensus,
                        iterable=chunk_jobs,
                        chunksize=chunksize_mp,
                    )
                )
                for svp in l
                if svp
            ]

        # svPrimitives = add_vcfIDs_to_svPrimitives(svPrimitives)

        svPrimitives.extend(
            add_adjacencies_to_svPrimitives_breakends(results)
        )  # all adjacencies are present here

        # sanity check: 1) all vcfIDs are unique
        if len(svPrimitives) != len(set([svp.vcfID for svp in svPrimitives])):
            raise ValueError(
                f"vcfIDs are not unique: {len(svPrimitives)} != {len(set([svp.vcfID for svp in svPrimitives]))}, all vcfIDs: {[svp.vcfID for svp in svPrimitives]}"
            )
    log.info("done!")
    return svPrimitives


def merge_genotypeMeasurements(
    this: dict[str, genotyping.GenotypeMeasurement],
    other: dict[str, genotyping.GenotypeMeasurement],
) -> dict[str, genotyping.GenotypeMeasurement]:
    # assert that both input dicts have the same keys
    # assert set(this.keys()) == set(other.keys()), f"Keys of both input dicts must be the same, but are {this.keys()} and {other.keys()}"
    # for each samplename, merge the genotypeMeasurements
    result = dict()
    all_keys = set(this.keys()) | set(other.keys())
    for samplename in all_keys:
        # if sampleID is only present in this or other, but not both, add the one that is present
        if samplename not in this:
            result[samplename] = other[samplename]
            continue
        if samplename not in other:
            result[samplename] = this[samplename]
            continue
        a = this[samplename]
        b = other[samplename]
        # new_position_on_consensus is the one that has the most supporting reads
        new_position_on_consensus = (
            a.position_on_consensus
            if len(a.supporting_reads) >= len(b.supporting_reads)
            else b.position_on_consensus
        )
        # estimated total depth - same case as above
        new_estimated_total_depth = (
            a.estimated_total_depth
            if len(a.supporting_reads) >= len(b.supporting_reads)
            else b.estimated_total_depth
        )
        new_distance_to_closest_depth_estimator = (
            a.distance_to_closest_depth_estimator
            if len(a.supporting_reads) >= len(b.supporting_reads)
            else b.distance_to_closest_depth_estimator
        )
        new_supporting_reads = list(set(a.supporting_reads + b.supporting_reads))
        new_neutral_reads = list(
            set(a.neutral_reads + b.neutral_reads) - set(new_supporting_reads)
        )
        new_contradicting_reads = list(
            set(a.contradicting_reads + b.contradicting_reads)
            - set(new_supporting_reads)
            - set(new_neutral_reads)
        )
        result[samplename] = genotyping.GenotypeMeasurement(
            samplename=samplename,
            position_on_consensus=new_position_on_consensus,
            estimated_total_depth=new_estimated_total_depth,
            distance_to_closest_depth_estimator=new_distance_to_closest_depth_estimator,
            supporting_reads=new_supporting_reads,
            neutral_reads=new_neutral_reads,
            contradicting_reads=new_contradicting_reads,
        )
    return result


def max_support_of_any_samplename(
    genotypeMeasurements: dict[str, genotyping.GenotypeMeasurement],
) -> int:
    try:
        return max(
            [
                len(genotypeMeasurement.supporting_reads)
                for genotypeMeasurement in genotypeMeasurements.values()
            ]
        )
    except ValueError:
        raise ValueError(
            f"genotypeMeasurements must not be empty. Given: {genotypeMeasurements}"
        )


def merge_svPrimitives_pair(
    first: SVprimitive_class.SVprimitive, second: SVprimitive_class.SVprimitive
) -> SVprimitive_class.SVprimitive:
    # check if first and second are fully initialized SVprimitives:

    assert isinstance(
        first, SVprimitive_class.SVprimitive
    ), "first must be an instance of SVprimitive"
    assert isinstance(
        second, SVprimitive_class.SVprimitive
    ), "second must be an instance of SVprimitive"
    # check if both are initialized by checking if their original_signals are not empty lists
    assert len(first.original_signals) > 0, "first must be initialized"
    assert len(second.original_signals) > 0, "second must be initialized"
    assert (
        len(first.original_cr_signals) > 0
    ), "first must be initialized with SV signals from the candidate region"
    assert (
        len(second.original_cr_signals) > 0
    ), "second must be initialized with SV signals from the candidate region"
    assert (
        len(first.genotypeMeasurements) > 0
    ), "first must be initialized with genotypeMeasurements"
    assert (
        len(second.genotypeMeasurements) > 0
    ), "second must be initialized with genotypeMeasurements"

    winner = (
        first
        if max_support_of_any_samplename(first.genotypeMeasurements)
        >= max_support_of_any_samplename(second.genotypeMeasurements)
        else second
    )
    new_genotypeMeasurements = merge_genotypeMeasurements(
        first.genotypeMeasurements, second.genotypeMeasurements
    )
    new_SV_primitive = deepcopy(first) if first == winner else second
    new_SV_primitive.genotypeMeasurements = new_genotypeMeasurements
    # merge original signals
    new_SV_primitive.original_signals = first.original_signals + second.original_signals
    # merge repeatIDs
    new_SV_primitive.repeatIDs = sorted(set(first.repeatIDs + second.repeatIDs))
    new_SV_primitive.reference_name = first.reference_name
    new_SV_primitive.distortion = first.distortion + second.distortion
    assert (
        type(new_SV_primitive.reference_name) == str
    ), "reference_name must be a string"
    return new_SV_primitive


# UNUSED RIGHT NOW!
def kmer_sketch_of_svPrimitive(
    svPrimitive: SVprimitive_class.SVprimitive, k: int, min_abs=2, min_rel=0.02
) -> set[int]:
    if svPrimitive.sv_type > 1:
        log.debug("svPrimitive is a BND. No kmer sketch is computed.")
        return set()
    letter_dict = util.compute_letter_dict("ACGTN")
    seqs = []
    if svPrimitive.sv_type == 0:
        seqs = svPrimitive.original_alt_sequences
    elif svPrimitive.sv_type == 1:
        seqs = svPrimitive.original_ref_sequences
    return util.kmer_sketch_from_strings(
        strings=seqs, letter_dict=letter_dict, k=k, min_abs=min_abs, min_rel=min_rel
    )


def cohen_d(x: list, y: list) -> float:
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(
        ((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof
    )


def can_merge_svPrimitives_pair(
    first: SVprimitive_class.SVprimitive,
    other: SVprimitive_class.SVprimitive,
    first_kmer_sketch: set | None = None,
    other_kmer_sketch: set | None = None,
    min_difference: int = 8,
    very_close: int = 50,
    close: int = 1_000,
    far: int = 10_000,
    max_cohens_d=-2.0,
    verbose: bool = False,
) -> bool:
    assert first.size >= 0 and other.size >= 0, "sizes must be >= 0"
    if verbose:
        print(
            f"trying to merge {first.sv_type} {first.size} {first.ref_start} with {other.sv_type} {other.size} {other.ref_start}"
        )
    if first.sv_type != other.sv_type:
        if verbose:
            print(
                f"can't merge because of different sv_types: {first.sv_type} and {other.sv_type}"
            )
        return False

    ### --- precompute all features --- ###
    first_is_bnd = first.sv_type > 1
    other_is_bnd = other.sv_type > 1
    distance = Interval(first.ref_start, first.ref_end).distance_to(
        Interval(other.ref_start, other.ref_end)
    )
    tolerance_a = int(round(np.mean(first.distortion)))
    tolerance_b = int(round(np.mean(other.distortion)))
    shared_reapeatIDs = set(first.repeatIDs).intersection(set(other.repeatIDs))
    if verbose:
        print(
            f"distance: {distance}, tolerance_a: {tolerance_a}, tolerance_b: {tolerance_b}, shared_reapeatIDs: {shared_reapeatIDs}, first_is_bnd: {first_is_bnd}, other_is_bnd: {other_is_bnd}"
        )

    ### --- RULES --- ###

    def _bnd_merge() -> bool:
        if first_is_bnd and other_is_bnd:
            if distance <= very_close:
                if verbose:
                    print(
                        f"can merge BNDs because of distance: {distance} <= {very_close}"
                    )
                return True
            else:
                if verbose:
                    print(
                        f"can't merge BNDs because of distance: {distance} > {very_close}"
                    )
                return False
        return False

    # ================================ 1) compare distances ================================ #

    def _distance_is_ok() -> bool:
        close_far = far if len(shared_reapeatIDs) > 0 else close
        if distance > tolerance_a + tolerance_b + close_far + min_difference:
            if verbose:
                print(
                    f"can't merge because of distance: {distance} > {tolerance_a} + {tolerance_b} + {close_far} + {min_difference}"
                )
            return False
        else:
            if verbose:
                print(
                    f"can merge because of distance: {distance} <= {tolerance_a} + {tolerance_b} + {close_far} + {min_difference}"
                )
            return True

    # ================================ 2) compare sizes ================================ #

    # adjust size estimates based on cohen's d if both distortion_sum_per_read have more than 2 values
    # if only one has more than 2 values, test if the other is less than 2 stds away
    # if both have less than 3 values, use the mean of the two and merge if their size difference is less than 8 or less than 2% of the smaller size
    def _similar_size() -> bool:
        if (
            len(first.distortion_sum_per_read) > 2
            and len(other.distortion_sum_per_read) > 2
        ):
            # if cohen's d of the sizes + values of distortion_sum_per_read is less than -1, they can't be merged
            dist_a = np.array(list(first.distortion_sum_per_read.values())) + (
                -first.size if first.sv_type == 1 else first.size
            )
            dist_b = np.array(list(other.distortion_sum_per_read.values())) + (
                -other.size if other.sv_type == 1 else other.size
            )
            d = cohen_d(dist_a, dist_b)
            if d < max_cohens_d:
                if verbose:
                    print(f"can't merge because of cohen's d: {d} < {max_cohens_d}")
                    print(f"distortion of A is {dist_a}, distortion of B is {dist_b}")
                return False
            else:
                if verbose:
                    print(f"can merge because of cohen's d: {d} >= {max_cohens_d}")
                    print(f"distortion of A is {dist_a}, distortion of B is {dist_b}")
                return True
        elif (
            len(first.distortion_sum_per_read) > 2
            or len(other.distortion_sum_per_read) > 2
        ):
            if len(first.distortion_sum_per_read) > 2:
                a = first
                b = other
            else:
                a = other
                b = first
            # check if b is in 2 stds of a
            std_a = np.std(list(a.distortion_sum_per_read.values()))
            if abs(b.size - a.size) > 2 * std_a:
                if verbose:
                    print(
                        f"can't merge because of size: {abs(b.size - a.size)} > 2*{std_a}"
                    )
                return False
            else:
                if verbose:
                    print(
                        f"can merge because of size: {abs(b.size - a.size)} <= 2*{std_a}"
                    )
                return True
        else:
            if not (
                abs(first.size - other.size) <= min_difference
                or abs(first.size - other.size) <= 0.02 * min(first.size, other.size)
            ):
                if verbose:
                    print(
                        f"can't merge because of size: {abs(first.size - other.size)} > {min_difference} or {abs(first.size - other.size)} > 0.02 * {min(first.size,other.size)}"
                    )
                return False
            else:
                if verbose:
                    print(
                        f"can merge because of size: {abs(first.size - other.size)} <= {min_difference} or {abs(first.size - other.size)} <= 0.02 * {min(first.size,other.size)}"
                    )
                return True

    # ================================ 3) compare k-mers ================================ #
    # ------ not used right now, because the giab benchmark worsened ----- #
    # -- works not so well and is very slow so far -- #
    # if kmer sketches are provided, first and other can only be merged, if they have similar kmer sketches
    def _kmer_similarity() -> bool:
        if (
            first_kmer_sketch
            and other_kmer_sketch
            and len(first_kmer_sketch) > 0
            and len(other_kmer_sketch) > 0
        ):
            # of one sketch is composed of only a few kmers, it is not reliable
            sketch_ratio: float = 0.0
            if len(first_kmer_sketch) < 10 and len(other_kmer_sketch) < 10:
                # both small
                sketch_ratio = util.intersection_ratio_of_smaller_set(
                    first_kmer_sketch, other_kmer_sketch
                )
                MIN_RATIO = 0.3
                if sketch_ratio > MIN_RATIO:
                    if verbose:
                        print(
                            f"can merge because of kmer sketch: {sketch_ratio} > {MIN_RATIO}"
                        )
                    return True
                else:
                    if verbose:
                        print(
                            f"can't merge because of kmer sketch: {sketch_ratio} <= {MIN_RATIO}"
                        )
                    return False
            if len(first_kmer_sketch) < 10 and len(other_kmer_sketch) > 10:
                sketch_ratio = util.intersection_ratio_of_smaller_set(
                    first_kmer_sketch, other_kmer_sketch
                )
                MIN_RATIO = 0.45
                if sketch_ratio > MIN_RATIO:
                    if verbose:
                        print(
                            f"can merge because of kmer sketch: {sketch_ratio} > {MIN_RATIO}"
                        )
                    return True
                else:
                    if verbose:
                        print(
                            f"can't merge because of kmer sketch: {sketch_ratio} <= {MIN_RATIO}"
                        )
                    return False
            else:  # both large
                MIN_RATIO = 0.6
                sketch_ratio = util.intersection_ratio_of_smaller_set(
                    first_kmer_sketch, other_kmer_sketch
                )
                if sketch_ratio > MIN_RATIO:
                    if verbose:
                        print(
                            f"can merge because of kmer sketch: {sketch_ratio} > {MIN_RATIO}"
                        )
                    return True
                else:
                    if verbose:
                        print(
                            f"can't merge because of kmer sketch: {sketch_ratio} <= {MIN_RATIO}"
                        )
                    return False
        else:
            if verbose:
                print("can merge: no kmer sketches provided")
            return True

    # ================================ CHECK CONDITIONS ================================ #

    similar_size = _similar_size()
    kmer_similarity = _kmer_similarity()
    distance_ok = _distance_is_ok()
    bnd_merge = _bnd_merge()

    if bnd_merge:
        if verbose:
            print("can merge because of BND merge rules")
        return True
    if similar_size and distance_ok and kmer_similarity:
        if verbose:
            print("can merge because of size, distance and kmer similarity")
        return True
    else:
        if verbose:
            print("can't merge because of size, distance or kmer similarity")
        return False


def merge_svPrimitives(
    svPrimitives: list[SVprimitive_class.SVprimitive],
    far: int = 10_000,
    close: int = 1_000,
    very_close: int = 50,
    max_cohens_d: float = -2.0,
    verbose: bool = False,
) -> list[SVprimitive_class.SVprimitive]:
    # --- merge svPrimitives that are most likely the same --- #
    # svps:list[SVprimitive_class.SVprimitive|None] = deepcopy(sorted(svPrimitives,key=lambda x: (x.chr,x.ref_start,x.ref_end)))
    SVPs = deepcopy(svPrimitives)  # don't change the input
    for i in range(len(SVPs) - 1, -1, -1):
        svPrimitive = SVPs[i]
        if svPrimitive is None:
            continue
        j = i - 1
        # if svPrimitive.sv_type <= 1:
        #     first_kmer_sketch = kmer_sketch_of_svPrimitive(svPrimitive,k=5)
        while j >= 0:
            other_svPrimitive = SVPs[j]
            if other_svPrimitive is None:
                j -= 1
                continue
            elif svPrimitive.chr != other_svPrimitive.chr:
                j -= 1
                continue
            elif svPrimitive.sv_type != other_svPrimitive.sv_type:
                j -= 1
                continue

            # TODO: don't merge if the k-mer sketches are very different.
            # count kmers that are more frequent than 5% of the total kmers in both sketches
            # if svPrimitive.sv_type <= 1:
            #     other_kmer_sketch = kmer_sketch_of_svPrimitive(other_svPrimitive,k=5)
            # else:
            #     first_kmer_sketch = None
            #     other_kmer_sketch = None
            can_merge = can_merge_svPrimitives_pair(
                first=svPrimitive,
                other=other_svPrimitive,
                far=far,
                close=close,
                very_close=very_close,
                max_cohens_d=max_cohens_d,
                verbose=verbose,
            )
            # first_kmer_sketch=first_kmer_sketch,
            # other_kmer_sketch=other_kmer_sketch)
            if can_merge:
                SVPs[i] = merge_svPrimitives_pair(svPrimitive, other_svPrimitive)
                SVPs[j] = None
                svPrimitive = SVPs[i]
            j -= 1
    return [svPrimitive for svPrimitive in SVPs if svPrimitive is not None]


def merge_svPrimitives_process_func(
    kwargs: dict,
) -> list[SVprimitive_class.SVprimitive]:
    return merge_svPrimitives(**kwargs)


def merge_svPrimitives_with_mp(
    svPrimitives: list[SVprimitive_class.SVprimitive],
    max_cohens_d: float,
    threads: int = 1,
) -> list[SVprimitive_class.SVprimitive]:
    # split up the list of svPrimitives by chromosome
    svPrimitives_by_chr = dict()
    for svPrimitive in svPrimitives:
        if svPrimitive.chr not in svPrimitives_by_chr:
            svPrimitives_by_chr[svPrimitive.chr] = []
        svPrimitives_by_chr[svPrimitive.chr].append(svPrimitive)
    jobs_args = [
        {"svPrimitives": svPrimitives_by_chr[chr], "max_cohens_d": max_cohens_d}
        for chr in svPrimitives_by_chr.keys()
    ]
    if threads > 1:
        with mp.Pool(min(mp.cpu_count(), threads)) as pool:
            job_results = list(
                pool.map(merge_svPrimitives_process_func, jobs_args, chunksize=1)
            )
    else:
        job_results = [merge_svPrimitives(**job_args) for job_args in jobs_args]
    return [svPrimitive for svPrimitives in job_results for svPrimitive in svPrimitives]


# TODO: needs testing!
def filter_svPrimitives_by_allow_list(
    allow_list: Path, svPrimitives: list[SVprimitive_class.SVprimitive]
) -> tuple[list[SVprimitive_class.SVprimitive], list[SVprimitive_class.SVprimitive]]:
    # build an intervaltree for each chr in the allowlist
    # iterate over all svPrimitives and check if they overlap any interval in the intervaltree (of the same chr)
    # 1) load the intervals of the form chr,start,end per line from allow_list
    with open(allow_list, "rt") as allow_list_handle:
        reader = csv.reader(allow_list_handle, delimiter="\t")
        allow_list_lines = [line for line in reader]
        chrs = {str(chr) for chr, start, end in allow_list_lines}
        allow_list_intervals = {chr: IntervalTree() for chr in chrs}
        for chr, start, end in allow_list_lines:
            allow_list_intervals[str(chr)].addi(int(start), int(end))
        # 2) iterate over all svPrimitives and check if they overlap any interval in the intervaltree (of the same chr)
        # if they overlap, add them to the results
        results: list[SVprimitive_class.SVprimitive] = []
        filtered_svPrimitives = []
        for svPrimitive in svPrimitives:
            if (
                svPrimitive.chr in allow_list_intervals.keys()
                and len(
                    allow_list_intervals[svPrimitive.chr][
                        svPrimitive.ref_start : svPrimitive.ref_end
                    ]
                )
                > 0
            ):
                results.append(svPrimitive)
            else:
                filtered_svPrimitives.append(svPrimitive)
    return results, filtered_svPrimitives


def create_svPrimitves_db(database: Path, name_reference: str):
    if "." in name_reference:
        name_reference = name_reference.replace(".", "_")
    with sqlite3.connect("file:" + str(database) + "?mode=rwc", uri=True) as conn:
        conn.execute(f"DROP TABLE IF EXISTS svPrimitives__{name_reference}")
        conn.execute(
            f"""CREATE TABLE svPrimitives__{name_reference} (
                        svPrimitiveID VARCHAR(30) PRIMARY KEY,
                        crID INTEGER,
                        name_reference VARCHAR(50),
                        consensusID VARCHAR(30),
                        svPrimitive BLOB)"""
        )
        conn.execute(
            f"CREATE INDEX idx_svPrimitives_consensusID ON svPrimitives__{name_reference} (consensusID)"
        )
        conn.execute(
            f"CREATE INDEX idx_svPrimitives_crID ON svPrimitives__{name_reference} (crID)"
        )
        conn.commit()


SQLITE_TIMEOUT = 10.0


# writes batches of fully annotated svPrimitves.
def write_svPrimitives_to_db(
    database: Path, data: list[SVprimitive_class.SVprimitive], name_reference: str
):
    # start_time = datetime.now()
    if "." in name_reference:
        name_reference = name_reference.replace(".", "_")
    cache = []
    for svPrimitive in data:
        if not svPrimitive.vcfID or svPrimitive.vcfID == "":
            raise ValueError(f"svPrimitiveID has no valid vcfID: {svPrimitive.vcfID}")
        consensusID = svPrimitive.consensusID
        crID = int(consensusID.split(".")[0])
        svPrimitiveID = (
            svPrimitive.vcfID
        )  # f"{name_dict[svPrimitive.sv_type]}.{count_dict[unify[svPrimitive.sv_type]]}"
        cache.append(
            (
                svPrimitiveID,
                consensusID,
                crID,
                name_reference,
                pickle.dumps(svPrimitive.unstructure()),
            )
        )
    # preprocess_time = datetime.now() - start_time
    if cache:
        query = f"INSERT INTO svPrimitives__{name_reference} (svPrimitiveID, consensusID, crID, name_reference, svPrimitive) VALUES (?, ?, ?, ?, ?)"
        with sqlite3.connect(str(database), timeout=SQLITE_TIMEOUT) as conn:
            c = conn.cursor()
            c.executemany(query, cache)
            conn.commit()
            c.close()
    # db_writing_time = datetime.now() - preprocess_time


def read_svPrimitives_from_db(
    database: Path,
    reference_name: str,
    consensusIDs: list[str] | None = None,
    crIDs: set[int] | None = None,
) -> dict[str, list[SVprimitive_class.SVprimitive]]:
    if "." in reference_name:
        reference_name = reference_name.replace(".", "_")
        log.warning(
            f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {reference_name}"
        )
    with sqlite3.connect(str(database)) as conn:
        table_name = f"svPrimitives__{reference_name}"
        svPrimitives = []
        c = conn.cursor()
        try:
            if crIDs:
                placeholders = ",".join(["?" for _ in crIDs])
                query = f"SELECT svPrimitive FROM {table_name} WHERE crID IN ({placeholders})"
                log.info(f"reading svPrimitives from database for {len(crIDs)} crIDs")
                c.execute(
                    query,
                    [
                        *crIDs,
                    ],
                )
            elif consensusIDs:
                placeholders = ",".join(["?" for _ in consensusIDs])
                query = f"SELECT svPrimitive FROM {table_name} WHERE consensusID IN ({placeholders})"
                log.info(
                    f"reading svPrimitives from database for {len(consensusIDs)} consensusIDs"
                )
                c.execute(
                    query,
                    [
                        *consensusIDs,
                    ],
                )
            else:
                query = f"SELECT svPrimitive FROM {table_name}"
                log.info("reading all svPrimitives from database")
                c.execute(query)

            for row in c.fetchall():
                # Use pickle.loads to deserialize the binary blob back into an svPrimitive object
                svPrimitives.append(
                    cattrs.structure(
                        pickle.loads(row[0]), SVprimitive_class.SVprimitive
                    )
                )
        except sqlite3.Error as e:
            log.error(f"Error reading svPrimitives from database: {e}")
            # print all tables found in the database
            c.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = c.fetchall()
            log.error(f"Tables in database: {tables}")
            c.close()
            raise e
        c.close()

    result = dict()
    for svPrimitive in svPrimitives:
        if svPrimitive.consensusID not in result:
            result[svPrimitive.consensusID] = []
        result[svPrimitive.consensusID].append(svPrimitive)
    return result


def yield_svPrimitives_from_db(
    database: Path, n_rows: int, reference_name: str
) -> Generator[SVprimitive_class.SVprimitive, None, None]:
    if "." in reference_name:
        reference_name = reference_name.replace(".", "_")
        log.warning(
            f"replacing '.' in name_reference with '_' to avoid sqlite3 error. New reference name is: {reference_name}"
        )
    with sqlite3.connect("file:" + str(database) + "?mode=ro", uri=True) as conn:
        c = conn.cursor()
        try:
            c.execute(f"SELECT svPrimitive FROM svPrimitives__{reference_name}")
            while True:
                rows = c.fetchmany(n_rows)
                if not rows:
                    break
                for row in rows:
                    yield cattrs.structure(
                        pickle.loads(row[0]), SVprimitive_class.SVprimitive
                    )
        except sqlite3.Error as e:
            log.error(f"Error reading svPrimitives from database: {e}")
            # print all tables found in the database
            c.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = c.fetchall()
            log.error(f"Tables in database: {tables}")
            raise e
        c.close()
