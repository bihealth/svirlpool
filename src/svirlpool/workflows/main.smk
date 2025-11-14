import json
from pathlib import Path

# =============================================================================
#  CONFIG
# =============================================================================

wdir = config["workdir"]

workdir: wdir

# =============================================================================
#  HELPER FUNCTIONS
# =============================================================================

def create_regions_from_fai(reference_path, output_path):
    """
    Create a BED file with all chromosomes/contigs from reference FASTA index.
    
    Args:
        reference_path: Path to reference FASTA file
        output_path: Path where regions.bed should be written
    """
    # Look for .fai file (could be .fasta.fai or .fa.fai)
    fai_path = None
    for ext in [".fai", ".fasta.fai", ".fa.fai"]:
        potential_fai = reference_path + ext if not reference_path.endswith(ext) else reference_path
        if Path(potential_fai).exists():
            fai_path = potential_fai
            break
    
    if fai_path is None:
        raise FileNotFoundError(f"Could not find .fai index file for reference: {reference_path}")
    
    # Read .fai file and create BED format (chrom, start, end)
    with open(fai_path, 'r') as fai, open(output_path, 'w') as bed:
        for line in fai:
            fields = line.strip().split('\t')
            chrom = fields[0]
            length = fields[1]
            bed.write(f"{chrom}\t0\t{length}\n")
    
    return output_path


# =============================================================================
#  INPUT FILES
# =============================================================================

samplename      = config["samplename"]
reference       = config["reference"]
alignments      = config["alignments"]
trf             = config["trf"]
# benchmark       = config["benchmark"]

# Handle optional regions parameter
if config.get("regions") is None or config["regions"] == "None":
    # Generate regions.bed from reference .fai file
    regions_file = str(Path(wdir) / "regions.bed")
    create_regions_from_fai(reference, regions_file)
    regions = regions_file
else:
    regions = config["regions"]

N_files_per_dir = config["N_files_per_dir"]
mononucleotides = config["mononucleotides"]
unique_regions  = config["unique_regions"]

# candidate regions
filter_absolute=config["filter_absolute"]
filter_normalized=config["filter_normalized"]
min_cr_size=config["min_cr_size"]
optimal_cr_size=config["optimal_cr_size"]
cr_merge_buffer=config["cr_merge_buffer"]

# consensus
lamassemble_mat = config["lamassemble_mat"]

# min mapq
min_mapq        = config["min_mapq"]
# consensus parameters
min_signal_size = config["min_signal_size"]

# copy number parameters
cn_dispersion   = config.get("cn_dispersion", 0.1)  # Default to 0.1 if not provided

# consensus alignments to proto svs parameters
merge_horizontally = config["cores_per_consensus"]

# svs to vcf parameters
min_alt_reads   = config["min_alt_reads"]
#min_alt_fraction= config["min_alt_fraction"]

# vcf to bed parameters
min_sv_size     = config["min_sv_size"]
cores_per_consensus = config["cores_per_consensus"]
#avg_doc = config["coverage"]
#doc_per_haplotype = avg_doc / 2.0

cores = workflow.cores #config["cores"]


wildcard_constraints:
    sampleID = r'\d+'


ref_name = Path(reference).stem


def get_crIDs(wildcards):
    checkpoint_output = checkpoints.candidate_regions_crs_to_containers.get(**wildcards).output[0]
    path_crIDs = 'crIDs.txt'
    crIDs = [int(line.strip()) for line in open(path_crIDs)]
    return [f"consensus/{str(int(crID / N_files_per_dir))}/consensus.{crID}.txt" for crID in crIDs]

# def get_coverage(wildcards):
#     with open("coverage.txt", "r") as f:
#         return f.read().strip()

# =============================================================================
#  RULES
# =============================================================================
rule all:
    input:
        "QC/crs.bed",
        "QC/rafs_indel_sizes.html",
        "QC/copy_number_tracks.png",
        'crs_containers.db',
        "consensus_containers.txt",
        "svirltile.db"
        
# -----------------------------------------------------------------------------
# SIGNAL PROCESSING
# -----------------------------------------------------------------------------

rule signalprocessing_alignments_to_rafs:
    input:
        alignments=alignments,
        regions=regions,
    output:
        file="rafs.tsv.gz",
        index="rafs.tsv.gz.tbi",
    params:
        samplename=samplename,
        min_signal_size=6,
        min_bnd_size=300,
        min_segment_size=250,
        min_mapq=min_mapq
    resources:
        mem_mb=cores*1024,
        runtime=240
    benchmark:
        "benchmarks/signalprocessing_alignments_to_rafs.txt"
    log:
        "logs/signalprocessing_alignments_to_rafs.log"
    threads:
        cores
    conda:
        "envs/svirlpool.yml"
    shell:
        """set x
        python3 -m svirlpool.signalprocessing.alignments_to_rafs \
        -a {input.alignments} \
        -s {params.samplename} \
        -r {input.regions} \
        -o {output.file} \
        --min-signal-size {params.min_signal_size} \
        --min-bnd-size {params.min_bnd_size} \
        --min-segment-size {params.min_segment_size} \
        --min-mapq {params.min_mapq} \
        -t {threads}"""


rule QC_rafs_indel_sizes:
    input:
        file="rafs.tsv.gz"
    output:
        "QC/rafs_indel_sizes.html"
    threads:
        1
    resources:
        mem_mb=8*1024,
        runtime=10
    shell:
        """python3 -m svirlpool.signalprocessing.rafs_indel_histograms \
        -i {input} \
        -f {output}"""

rule signalprocessing_filter_rafs_sv_signals:
    input:
        file="rafs.tsv.gz",
        index="rafs.tsv.gz.tbi"
    output:
        file="rafs.filtered_sv_signals.tsv.gz",
        index="rafs.filtered_sv_signals.tsv.gz.tbi"
    params:
        reference=reference,
        threshold=1.0
    threads:
        cores
    resources:
        mem_mb=8*1024,
        runtime=15
    shell:
        """python3 -m svirlpool.signalprocessing.filter_rafs_sv_signals \
        -i {input.file} \
        -r {params.reference} \
        -o {output.file} \
        -t {threads} \
        --threshold {params.threshold}"""


# rule signalprocessing_adjust_effective_intervals:
#     input:
#         file="rafs.filtered_sv_signals.tsv.gz",
#         index="rafs.filtered_sv_signals.tsv.gz.tbi"
#     output:
#         file="rafs.effective_intervals.tsv.gz",
#     params:
#         reference=reference,
#         unique_regions=unique_regions,
#         min_raf_size=1000,
#         rafs_dropped="rafs.filtered_sv_signals.dropped.tsv.gz"
#     threads:
#         cores
#     log:
#         "logs/adjust_effective_intervals.log"
#     resources:
#         mem_mb=cores*1024,
#         runtime=60
#     benchmark:
#         "benchmarks/adjust_effective_intervals.txt"
#     shell:
#         """python3 -m svirlpool.scripts.filter_rafs \
#         -i {input.file} \
#         -r {params.reference} \
#         -o {output.file} \
#         -u {params.unique_regions} \
#         -t {threads} \
#         --min-raf-size {params.min_raf_size} \
#         --rafs-dropped {params.rafs_dropped}"""


rule signalprocessing_filter_rafs_indels:
    input:
        file="rafs.filtered_sv_signals.tsv.gz",
        #index="rafs.effective_intervals.tsv.gz.tbi"
    output:
        file="rafs.effective_intervals.filtered_indels.tsv.gz",
        dropped="rafs.effective_intervals.filtered_indels.dropped.tsv.gz"
    params:
        reference=reference,
        rafs_cache_size=200,
        multiplier=8.0,
    threads:
        cores
    resources:
        mem_mb=cores*1024,
        runtime=60
    benchmark:
        "benchmarks/signalprocessing_filter_rafs_indels.txt"
    log:
        "logs/signalprocessing_filter_rafs_indels.log"
    shell:
        """python3 -m svirlpool.signalprocessing.filter_rafs_indels \
        -i {input.file} \
        -o {output.file} \
        -r {params.reference} \
        -t {threads} \
        --multiplier {params.multiplier} \
        --dropped {output.dropped}"""


rule signalprocessing_rafs_to_coverage:
    input:
        file="rafs.effective_intervals.filtered_indels.tsv.gz",
    output:
        file="coverages.db"
    threads:
        1
    resources:
        mem_mb=cores*1024,
        runtime=30
    benchmark:
        "benchmarks/signalprocessing_rafs_to_coverage.txt"
    log:
        "logs/signalprocessing_rafs_to_coverage.log"
    shell:
        """python3 -m svirlpool.signalprocessing.rafs_to_coverage \
        -i {input.file} \
        -o {output.file}"""


rule signalprocessing_generate_copynumber_tracks:
    input:
        coverage_db="coverages.db"
    output:
        bed="copy_number_tracks.bed.gz",
        tbi="copy_number_tracks.bed.gz.tbi",
        db="copy_number_tracks.db",
        plot="QC/copy_number_tracks.png"
    params:
        bed_uncompressed="copy_number_tracks.bed",
        regions=regions,
        reference_fai=reference + ".fai",
        chr_filterlist=[],  # Can be configured via config if needed
        dispersion=cn_dispersion
    threads:
        cores
    resources:
        mem_mb=cores*2048,
        runtime=60
    benchmark:
        "benchmarks/signalprocessing_generate_copynumber_tracks.txt"
    log:
        "logs/signalprocessing_generate_copynumber_tracks.log"
    conda:
        "envs/svirlpool.yml"
    shell:
        """
        set -x
        python3 -m svirlpool.signalprocessing.copynumber_tracks \
        --coverage-db {input.coverage_db} \
        --output-bed {params.bed_uncompressed} \
        --output-db {output.db} \
        --threads {threads} \
        --plot-output {output.plot} \
        --reference-fai {params.reference_fai} \
        --regions {params.regions} \
        --dispersion {params.dispersion} \
        --log-level INFO 2>&1 | tee {log}"""


# rule compute_coverage:
#     input:
#         "coverages.db"
#     output:
#         "coverage.txt"
#     params:
#         regions=regions
#     threads:
#         1
#     resources:
#         mem_mb=2*1024,
#         runtime=10
#     benchmark:
#         "benchmarks/compute_coverage.txt"
#     log:
#         "logs/compute_coverage.log"
#     conda:
#         "envs/svirlpool.yml"
#     shell:
#         """python3 -c "
# from svirlpool.scripts.covtree import median_binned_total_coverage
# coverage = median_binned_total_coverage('{input}', regions='{params.regions}')
# with open('{output}', 'w') as f:
#     f.write(str(coverage))
# " """


rule signalprocessing_rafs_to_signaldepths:
    input:
        file="rafs.effective_intervals.filtered_indels.tsv.gz"
    output:
        file="signaldepths.tsv.gz",
        #index="signaldepths.tsv.gz.tbi"
    threads:
        cores
    resources:
        mem_mb=cores*1024,
        runtime=60
    benchmark:
        "benchmarks/signalprocessing_rafs_to_signaldepths.txt"
    log:
        "logs/signalprocessing_rafs_to_signaldepths.log"
    shell:
        """python3 -m svirlpool.signalprocessing.rafs_to_signaldepths \
        -i {input.file} \
        -o {output.file} \
        -t {threads} \
        --log-file {log}"""


rule signalprocessing_filter_mononucleotides:
    input:
        signaldepths="signaldepths.tsv.gz",
    output:
        outf='filtered_mononucleotides.signaldepths.tsv.gz',
    params:
        margin=5,
        mononucleotides=mononucleotides,
        reference=reference
    resources:
        mem_mb=cores*1024,
        runtime=20
    benchmark:
        "benchmarks/signalprocessing_filter_mononucleotides.txt"
    threads:
        cores
    shell:
        """python3 -m svirlpool.signalprocessing.filter_signaldepths_mononucleotides \
        -m {params.mononucleotides} \
        -s {input.signaldepths} \
        -r {params.reference} \
        -o {output.outf} \
        -g {params.margin} \
        -t {threads}"""


rule signalprocessing_depth_of_coverage_to_signalstrength:
    input:
        signaldepths='filtered_mononucleotides.signaldepths.tsv.gz'
    output:
        out='signalstrength.tsv.gz'
    params:
        ref=reference,
        rep=trf,
        radius=1000,
        flatness=30.0,
    resources:
        mem_mb=cores*1024,
        runtime=90
    benchmark:
        "benchmarks/signalprocessing_depth_of_coverage_to_signalstrength.txt"
    threads:
        cores
    conda:
        "envs/svirlpool.yml"
    shell:
        """python3 -m svirlpool.signalprocessing.signaldepths_to_signalstrength \
        -i {input.signaldepths} \
        -r {params.ref} \
        -o {output.out} \
        -R {params.rep} \
        -t {threads} \
        -s {params.flatness} \
        -k {params.radius}"""

# -----------------------------------------------------------------------------
# CANDIDATE REGIONS
# -----------------------------------------------------------------------------


rule candidate_regions_signalstrength_to_crs:
    input:
        signalstrength='signalstrength.tsv.gz'
    output:
        out='crs.db',
        bedgraph='QC/crs.signal.bedgraph'
    params:
        ref=reference,
        rep=trf,
        filter_absolute=filter_absolute,
        filter_normalized=filter_normalized,
        min_cr_size=min_cr_size,
        bedgraph='QC/crs',
        buffer=cr_merge_buffer,
        dropped='crs.dropped.tsv.gz'
    resources:
        mem_mb=cores*2*1024,
        runtime=60
    benchmark:
        "benchmarks/candidate_regions_signalstrength_to_crs.txt"
    threads:
        cores
    conda:
        "envs/svirlpool.yml"
    shell:
        """set -x
        python3 -m svirlpool.candidateregions.signalstrength_to_crs \
        -i {input.signalstrength} \
        -r {params.ref} \
        -o {output.out} \
        -u {params.rep} \
        -t {threads} \
        -b {params.buffer} \
        --dropped {params.dropped} \
        --bedgraph {params.bedgraph} \
        --filter-absolute {params.filter_absolute} \
        --filter-normalized {params.filter_normalized} \
        --min-cr-size {params.min_cr_size}"""

rule candidate_regions_QC_crs_to_bed:
    input:
        'crs.db'
    output:
        "QC/crs.bed"
    resources:
        mem_mb=cores*2*1024,
        runtime=10
    benchmark:
        "benchmarks/candidate_regions_QC_crs_to_bed.txt"
    threads:
        1
    conda:
        "envs/svirlpool.yml"
    shell:
        """python3 -m svirlpool.candidateregions.crs_to_bed \
        -i {input} \
        -o {output}"""


checkpoint candidate_regions_crs_to_containers:
    input:
        crs='crs.db'
    output:
        containers='crs_containers.db',
        #hist='QC/crs.hist.png',
        crIDs_file='crIDs.txt'
    resources:
        mem_mb=cores*8*1024,
        runtime=40
    benchmark:
        "benchmarks/candidate_regions_crs_to_containers.txt"
    threads:
        1
    conda:
        "envs/svirlpool.yml"
    shell:
        """python3 -m svirlpool.candidateregions.crs_to_containers_db \
        -i {input.crs} \
        -o {output.containers} \
        --crIDs {output.crIDs_file}"""

# # -----------------------------------------------------------------------------
# # CONSENSUS
# # -----------------------------------------------------------------------------

# # assemble the cut reads
rule consensus_consensus:
    input:
        containers='crs_containers.db',
        copynumbertracks='copy_number_tracks.bed.gz',
    output:
        container="consensus/{crIDdir}/consensus.{crID}.txt",
        log="consensus/{crIDdir}/consensus.{crID}.log",
    params:
        alignments=alignments,
        samplename=samplename,
        lamassemble_mat=lamassemble_mat,
    threads:
        cores_per_consensus
    conda:
        "envs/svirlpool.yml"
    resources:
        mem_mb=2*1024,
        runtime=15
    benchmark:
        "benchmarks/consensus/consensus.{crID}.{crIDdir}.txt"
    shell:
        """set -x
        python3 -m svirlpool.localassembly.consensus \
        -s {params.samplename} \
        -cn {input.copynumbertracks} \
        --lamassemble-mat {params.lamassemble_mat} \
        -i {input.containers} \
        -a {params.alignments} \
        -o {output.container} \
        -t {threads} \
        -c {wildcards.crID} \
        --logfile {output.log}"""
        # --verbose"""


rule consensus_aggregate_consensuses:
    input:
        consensus_paths=get_crIDs,
    output:
        "consensus_containers.txt",
    threads:
        1
    resources:
        mem_mb=16*1024,
        runtime=5
    benchmark:
        "benchmarks/consensus/consensus_aggregate.txt"
    shell:
        """for file in {input.consensus_paths}; do cat $file >> {output}; done"""



rule consensus_align:
    input:
        consensus="consensus_containers.txt",
        #coverage="coverages.db",
    output:
        consensus_fasta="consensus.fasta",
        consensus_alignments="consensus.bam",
        svpatterns="svpatterns.db",
    params:
        ref=reference,
        repeats=trf,
        samplename=samplename
    threads:
        cores
    resources:
        mem_mb=16*1024,
        runtime=60
    benchmark:
        "benchmarks/consensus/consensus_align_to_initial_reference.txt"
    log:
        "logs/consensus_align_to_initial_reference.log"
    shell:
        """
        set -x
        python -m svirlpool.localassembly.consensus_align \
            --samplename {params.samplename} \
            --consensus {input.consensus} \
            --reference {params.ref} \
            --trf {params.repeats} \
            --output-consensus-fasta {output.consensus_fasta} \
            --output-consensus-alignments {output.consensus_alignments} \
            --output-svpatterns {output.svpatterns} \
            --threads {threads}"""


rule to_svirltile:
    input:
        svpatterns="svpatterns.db",
        coverage="coverages.db",
        consensus="consensus.fasta",
        copynumber="copy_number_tracks.db"
    params:
        samplename=samplename
    output:
        svirltile_db="svirltile.db",
    threads:
        1
    resources:
        mem_mb=32*1024,
        runtime=10
    benchmark:
        "benchmarks/consensus/to_svirltile.txt"
    log:
        "logs/consensus/to_svirltile.log"
    shell:
        """
        set -x
        python -m svirlpool.localassembly.svirltile \
        --samplename {params.samplename} \
        --consensus {input.consensus} \
        --svpatterns-db {input.svpatterns} \
        --coverage-db {input.coverage} \
        --copynumber-db {input.copynumber} \
        --output-db {output.svirltile_db} \
        --log-file {log} \
        --log-level INFO """
