#!/bin/bash

# Set a name for the job (-J or --job-name).
#SBATCH --job-name=minidel

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=3

# Force allocation of the two cores on ONE node.

# Set the total memory. Units can be given in T|G|M|K.
#SBATCH --mem=16G

# Optionally, set the partition to be used (-p or --partition).
#SBATCH --partition=short

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=30:00

WDIR=$1
PREFIX=$2
REF=$3
BENCHMARK=$4

arg_O0=$5
arg_O1=$6

argID=$7

echo "all inputs are:"
echo $WDIR $PREFIX $REF $BENCHMARK $arg_O0 $arg_O1 $argID
echo "cat and align consensus sequences to reference"
python3 -m scripts.consensus_to_variants \
    -i consensus/*.fasta \
    -r $REF \
    -f fastq/*.fastq \
    -s $PREFIX.signaldepth \
    -d $PREFIX.sampledicts.json \
    -a $PREFIX.$argID.consensus.bam \
    -o $PREFIX.$argID.svs \
    -x \"-O$arg_O0,$arg_O1\"

echo "call svs from consensus-to-reference alignments"
python3 -m scripts.svs_to_vcf \
    -i $PREFIX.$argID.svs \
    -c consensus/*fasta \
    -r $REF \
    -s $PREFIX.sampledicts.json \
    -o $PREFIX.$argID.vcf

echo "benchmark in truvari bench"
rm -r -f truvari.$argID
truvari bench \
    -b $BENCHMARK \
    -c $PREFIX.$argID.vcf.gz \
    -o truvari.$argID \
    -f $REF \
    --includebed minidel.targets.bed \
    --passonly

RESULT=$(python3 -c "import json; json_data = open('truvari.$argID/summary.json'); data = json.load(json_data); json_data.close(); print(data['f1'],data['precision'],data['recall'],sep='\t')")

echo $arg_O0,$arg_O1 $RESULT >$argID.results.txt
