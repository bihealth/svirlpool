# this script runs truvari and creates bed files from the TP, FN, FP vcf files
# it filters the input file for calls in HG002 and then runs truvari

# add input parsing and print how to use the script if the input is not correct
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <input_vcf> <benchmark_vcf> <regions_bed> <reference> <min_size> <samplename> <output>"
    exit 1
fi
# check if min_size is a positive int greater than 1
if ! [[ $5 =~ ^[0-9]+$ ]]; then
    echo "Error: min_size must be a positive integer"
    exit 1
fi
# check if input_vcf and benchmark_vcf end with 'vcf.gz'
if [[ $1 != *vcf.gz ]] || [[ $2 != *vcf.gz ]]; then
    echo "Error: input_vcf and benchmark_vcf must end with 'vcf.gz'"
    exit 1
fi
# check if regions_bed ends with 'bed'
if [[ $3 != *bed ]]; then
    echo "Error: regions_bed must end with 'bed'"
    exit 1
fi
# check if there is a .fai file for the reference (check if a file exists that ends with '.fai')
if ! ls $4.fai 1>/dev/null 2>&1; then
    echo "Error: reference must have a .fai file"
    exit 1
fi

# -----------------------------------------------------------------------------
#  set input variables
input_vcf=$1
benchmark_vcf=$2
regions_bed=$3
reference=$4
min_size=$5
samplename=$6
output=$7

OUTDIR=$output/truvari
TMP_VCF=$output/tmp.vcf.gz
# check if the output directory exists, if so, remove it and create a new one,
# else just create a new one
if [ -d $output ]; then
    rm -r $output
fi
mkdir $output

# -----------------------------------------------------------------------------
# filter input vcf for sample $samplename
set -x
bcftools view -s $samplename $input_vcf -Oz -o $TMP_VCF
# index the tmp vcf file
tabix -f -p vcf $TMP_VCF

# -----------------------------------------------------------------------------
#  run truvari

truvari bench \
    -b $benchmark_vcf \
    -c $TMP_VCF \
    -o $OUTDIR \
    -f $reference \
    --passonly \
    --includebed $regions_bed

# -----------------------------------------------------------------------------
# remove temporary vcf file
#rm $TMP_VCF

# -----------------------------------------------------------------------------
#  run vcf_to_bed

python3 -m scripts.vcf_to_bed \
    -i $OUTDIR/tp-base.vcf.gz -o $OUTDIR/tp-base.bed \
    --min_size $min_size
python3 -m scripts.vcf_to_bed \
    -i $OUTDIR/fp.vcf.gz -o $OUTDIR/fp.bed \
    --min_size $min_size
python3 -m scripts.vcf_to_bed \
    -i $OUTDIR/fn.vcf.gz -o $OUTDIR/fn.bed \
    --min_size $min_size
