# Svirlpool: structural variant detection from long read sequencing by local assembly

Our [preprint](https://www.biorxiv.org/content/10.1101/2025.11.03.686231v1) is now on BioRXiv!

<img src="media/illustration.png" align="right" alt="Svirlpool illustration" width="300" height="300">

- License: MIT
- Contact: [Vinzenz May](mailto:vinzenz.may@bih-charite.de)

**Motivation**

Long-Read Sequencing (LRS) promises great improvements in the detection of structural genome variants (SVs). However, existing methods are lacking in key areas such as the reliable detection of inserted sequence, precise genotyping of variants, and reproducible calling of variants across multiple samples. Here, we present our method Svirlpool, that is aimed at the analysis of Oxford Nanopore Technologies (ONT) sequencing data. Svirlpool uses local assembly of candidate SV regions to obtain high-quality consensus sequences.

**Results**

Svirlpool obtains competitive results to the leading method Sniffles on the widely used Genome in a Bottle (GiaB) benchmark data sets. On trio data, however, Svirlpools shows a clear favorable performance in terms of mendelian consistency. This indicates that Svirlpool shows great promise in clinical applications and beyond benchmark datasets.

## Installing & Running

You can install Svirlpool from sources, see "Developer Notes" below.
Alternatively, we provide a Docker image at `ghcr.io/bihealth/svirlpool` (e.g., `:main` for the latest version).
You can also convert that image to Singularity if needed on a HPC system (e.g., `singularity build svirlpool.sif docker://ghcr.io/bihealth/svirlpool:main`).
The example below explains how to run via Docker.

### Prerequisites

Install Git:

```
# Ubuntu 24.04+
sudo apt install -y git git-lfs
```

Then, install Docker following [these official instructions](https://docs.docker.com/engine/install//).
We assume that you can call `docker` without `sudo` in the following.
Otherwise, you may have to prefix the `docker` commands below with `sudo`.

### Getting Example Data

As the example data is stored in the source repository, we obtain the sources via git:

```
git clone git@github.com:bihealth/svirlpool.git
```

### Get Docker Image

You can explicitely download the docker image now (or it will be downloaded on the fly when used):

```
docker pull ghcr.io/bihealth/svirlpool:main
```

### Run Docker Image

First, run Svirlpool to create the "svirltile" file.

```
# create working directory
mkdir -p /tmp/workdir/result
# now run
docker run \
    --rm \
    -v $(realpath .):/data \
    -v /tmp/workdir/result:/tmp/workdir/result \
    -w /data \
    ghcr.io/bihealth/svirlpool:sha-c8ef958 \
        svirlpool run \
            --threads 1 \
            --samplename muc1test \
            --workdir /tmp/workdir/result \
            --output /tmp/workdir/result/svirltile.db \
            --alignments examples/muc1/data/muc1.bam \
            --reference examples/muc1/data/muc1.fa \
            --trf examples/muc1/data/muc1.trf.bed \
            --mononucleotides examples/muc1/data/muc1.mononucleotides.lt6.bed \
            --lamassemble-mat data/lamassemble-mats/promethion.mat
```

The directory `/tmp/workdir/result` will contain a number of files. The most important one is the file `svirltile.db` that is necessary for the subsequent (potentially joint) SV calling and creation of a VCF file.

Let us now run the SV calling:

```
docker run \
    --rm \
    -v $(realpath .):/data \
    -v /tmp/workdir/result:/tmp/workdir/result \
    -w /data \
    ghcr.io/bihealth/svirlpool:sha-c8ef958 \
        svirlpool sv-calling \
            --threads 1 \
            --reference examples/muc1/data/muc1.fa \
            --input /tmp/workdir/result/svirltile.db \
            --output /tmp/workdir/result/variants.vcf.gz
```

We will end up with the SV file `variants.vcf.gz` which contains the SV calls, see below for an excerpt.

```
# zcat /tmp/workdir/result/variants.vcf.gz | cut -b 1-140
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20251104
##reference=file:///home/mholtgrewe/Development/svirlpool/examples/muc1/data/muc1.fa
##sequences=file:///tmp/workdir/result/variants.variants.fasta
##contig=<ID=1,length=300436>
##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of informative reads.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion, DUP=Duplication, INV=Inversion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=PASS_ALTREADS,Number=1,Type=String,Description="Passed alt reads threshold">
##INFO=<ID=pass_GQ,Number=1,Type=String,Description="Passed Genotype precision threshold">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=CONSENSUSIDs,Number=1,Type=String,Description="ID of the consensus that this SV originiates from. Other consensus sequences can a
##INFO=<ID=SEQ_ID,Number=1,Type=String,Description="ID of sequence in companion FASTA file for symbolic alleles">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend; Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=TC,Number=1,Type=Integer,Description="Total coverage">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=GP,Number=1,Type=Float,Description="Genotype probability">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  muc1test
1       114551  DEL.4   GGGGATTACAGGCGTGAGCCACCGCTCCTAGCGAAAAACATTTTTT  G       60      LowQual PASS_ALTREADS=False;pass_GQ=True;SVTYPE=DEL;END=114596;SVLEN=45;C
1       159296  INS.0   C       Cggctggagcccgagggccggcctggtgtccgggggccgaggtcggacaccgcctggctgggggcggtggagccgccgggccggcctggtgtccgggggccgaggtgacaccgtgggctgggg
1       159351  INS.3   G       Gcggtggagcccgggccggctggtgtccgggggccgaggtgacaccgtgggctggggcggtggagccgggccggcctggtgtccggggccgaggtgacaccgtgggctggggcggtggagccc
1       159704  INS.1   G       Gggctggggcggctgcagcccggggccggcctgctctccggggccgaggtgacaccgccgtgggctgcggcgggcggcgtggagcccggggcccggctgctcctatcttccgggccgaggtgt
1       159707  INS.2   T       Tggggcggtggagccgcgggccggcctggtgtccggggccgaggtgacaccgtgggctgggcggcgtggagcccggggccggcctggtgtccggggccgaggtgacaccgtgggctgggcggc
```

## Developer Notes

### Prerequisites

First of all, you will need to have Git (to get the source code) and Git LFS (to get the example files) installed.

```
# Ubuntu 24.04+
sudo apt install -y git git-lfs
```

Next, we install [Pixi](https://pixi.sh/latest/), for installing Python, the Python libraries we depend on.
Also, Pixi allows us to install non-Python packages from the `conda-forge` and `bioconda` repositories.

```
curl -fsSL https://pixi.sh/install.sh | bash
```

### Steps

Clone repository:

```
git clone git@github.com:bihealth/svirlpool.git
cd svirlpool
```

First, run Svirlpool to create the "svirltile" file.

```
# create working directory
mkdir -p /tmp/workdir/result
# now run
pixi run \
    svirlpool run \
        --threads 1 \
        --samplename muc1test \
        --workdir /tmp/workdir/result \
        --output /tmp/workdir/result/svirltile.db \
        --alignments examples/muc1/data/muc1.bam \
        --reference examples/muc1/data/muc1.fa \
        --trf examples/muc1/data/muc1.trf.bed \
        --mononucleotides examples/muc1/data/muc1.mononucleotides.lt6.bed \
        --lamassemble-mat data/lamassemble-mats/promethion.mat
```

The directory `/tmp/workdir/result` will contain a number of files. The most important one is the file `svirltile.db` that is necessary for the subsequent (potentially joint) SV calling and creation of a VCF file.

Let us now run the SV calling:

```
pixi run \
    svirlpool sv-calling \
        --threads 1 \
        --reference examples/muc1/data/muc1.fa \
        --input /tmp/workdir/result/svirltile.db \
        --output /tmp/workdir/result/variants.vcf.gz
```

We will end up with the SV file `variants.vcf.gz` which contains the SV calls, same as shown in the Docker example above.

## How to call SVs with Svirlpool on my real data

You will need the get the prefab data (I), then create the svirltile(s) of your samples (II), then call SVs from the svirltiles (III) to generate the final vcf file.

Your starting point is the svirlpool directory.

### I prefab data
Some files are needed to run svirlpool. What you need:
1) read alignments - indexed .bam files that were generated with minimap2. They need to have read DNA sequences and sequence quality scores.
2) an indexed reference genome - in .fasta format; indexed with samtools faidx
3) the matrices for lamassemble. They are here in the svirlool-repository. set: SVIRLPOOLDIR=$(realpath .)
4) annotations files - you can get them from: https://github.com/bihealth/svirlpool-data
  - download the data repository and cd into it. Set: DATADIR=$(realpath .)

let's say you have the samples HG002 and HG003 from the GiaB test data set on HG19. You now set the variables and then run svirlpool:

`cd $SVIRLPOOLDIR`

`pixi shell`

```
REFERENCE={hs37d5.fa}
THREADS=16

TRF=$DATADIR/pbsv-annotations/human_hs37d5.trf.bed
MAT=$SVIRLPOOLDIR/data/lamassemble-mats/promethion.mat
MNNTS=$DATADIR/HG19/hs37d5.mononucleotides.lt6.bed.gz
```

Now cd into a directory of your choice, where your data shall be processed, e.g.

`cd dataanalysis`

### II generate svirltiles
#### run the son HG002:
```
SAMPLENAME=HG002
ALIGNMENTS={ALIGNMENTS.HG002.BAM}

svirlpool run \
    --samplename $SAMPLENAME \
    --workdir $SAMPLENAME \
    --alignments $ALIGNMENTS \
    --reference $REFERENCE \
    --trf $TRF \
    --lamassemble-mat $MAT \
    --mononucleotides $ MNNTS \
    --threads $THREADS \
    --min-sv-size 30
```
#### run the father HG003:
```
SAMPLENAME=HG003
ALIGNMENTS={ALIGNMENTS.HG003.BAM}

svirlpool run \
    --samplename $SAMPLENAME \
    --workdir $SAMPLENAME \
    --alignments $ALIGNMENTS \
    --reference $REFERENCE \
    --trf $TRF \
    --lamassemble-mat $MAT \
    --mononucleotides $ MNNTS \
    --threads $THREADS \
    --min-sv-size 30
```
### III vcf

```
svirlpool sv-calling \
    --input $HG002/svirltile.db $HG003/svirltile.db \
    --reference $REFERENCE \
    --output family.vcf.gz \
    --sv-types DEL INS \
    --min-sv-size 50
```

The resulting file should be family.vcf.gz.


### Hints

Open VS Code with pix-installed `dev` environment

```
pixi run -e dev code .
```

Run code formatting

```
make fix
```

Run formatting, lints, other static checks:

```
make check
```

Run tests:

```
make check
```

Or all in one:

```
make fix check test
```
