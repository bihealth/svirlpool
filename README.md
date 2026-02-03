# Svirlpool: structural variant detection from long read sequencing by local assembly

Our [preprint](https://www.biorxiv.org/content/10.1101/2025.11.03.686231v1) is now on BioRXiv!

<img src="media/illustration.png" align="right" alt="Svirlpool illustration" width="300" height="300">

* **License:** MIT
* **Contact:** [Vinzenz May](mailto:vinzenz.may@bih-charite.de)

---

## Overview

**Motivation**

Long-Read Sequencing (LRS) promises great improvements in the detection of structural genome variants (SVs). However, existing methods are lacking in key areas such as the reliable detection of inserted sequence, precise genotyping of variants, and reproducible calling of variants across multiple samples. **Svirlpool** targets Oxford Nanopore Technologies (ONT) sequencing data using local assembly of candidate SV regions to obtain high-quality consensus sequences.

**Results**

Svirlpool obtains competitive results to leading methods like Sniffles on Genome in a Bottle (GiaB) benchmarks. On trio data, Svirlpool shows favorable performance in Mendelian consistency, indicating great promise for clinical applications.

---

## Table of Contents

1. [Installation](https://www.google.com/search?q=%23installation)
    * [Option A: Docker or Singularity](https://www.google.com/search?q=%23option-a-docker-or-singularity)
    * [Option B: From Source (Pixi)](https://www.google.com/search?q=%23option-b-from-source-pixi)
2. [Quick Start (Example Data)](https://www.google.com/search?q=%23quick-start-example-data)
3. [Workflow for Real Data](https://www.google.com/search?q=%23workflow-for-real-data)
    * [I. Required Input Files](https://www.google.com/search?q=%23i-required-input-files)
    * [II. Step-by-Step Execution](https://www.google.com/search?q=%23ii-step-by-step-execution)
4. [Developer & Formatting Hints](https://www.google.com/search?q=%23developer--formatting-hints)
5. [Release Management](https://www.google.com/search?q=%23release-management)

---

## Installation

You can run Svirlpool either via pre-built containers or by installing from the source code. Both methods are fully supported.

### Option A: Docker or Singularity

1. **Install Docker:** Follow the [official instructions](https://docs.docker.com/engine/install/).

2. **Pull Image:**

```bash
docker pull ghcr.io/bihealth/svirlpool:main
```

3. **Singularity (Alternative):** If on an HPC, convert the image:

```bash
singularity build svirlpool.sif docker://ghcr.io/bihealth/svirlpool:main
```

### Option B: From Source (Pixi)

1. **Install Prerequisites:**

```bash
# Ubuntu 24.04+
sudo apt install -y git git-lfs
curl -fsSL https://pixi.sh/install.sh | bash
```

2. **Clone and Setup:**

```bash
git clone git@github.com:bihealth/svirlpool.git
cd svirlpool
# Pixi will manage the environment automatically on the first run
```

---

## Quick Start (Example Data)

This example uses a small MUC1 test dataset to demonstrate the two-step calling process.
We demonstrate how to run with Docker.
You can replace the `docker run ... svirlpool` part with `pixi run svirlpool` to run it using pixi rather than Docker.

### 1. Generate Svirltile

```bash
# Create working directory
mkdir -p /tmp/workdir/result

# Run using Docker (or replace with 'pixi run svirlpool ...')
docker run --rm -v $(realpath .):/data -v /tmp/workdir/result:/tmp/workdir/result -w /data \
    ghcr.io/bihealth/svirlpool:main \
    svirlpool run \
        --threads 1 --samplename muc1test --workdir /tmp/workdir/result \
        --output /tmp/workdir/result/svirltile.db \
        --alignments examples/muc1/data/muc1.bam \
        --reference examples/muc1/data/muc1.fa \
        --trf examples/muc1/data/muc1.trf.bed \
        --mononucleotides examples/muc1/data/muc1.mononucleotides.lt6.bed \
        --lamassemble-mat data/lamassemble-mats/promethion.mat
```

### 2. Generate VCF

```bash
docker run \
    --rm -v $(realpath .):/data -v /tmp/workdir/result:/tmp/workdir/result -w /data \
    ghcr.io/bihealth/svirlpool:main \
    svirlpool sv-calling \
        --threads 1 --reference examples/muc1/data/muc1.fa \
        --input /tmp/workdir/result/svirltile.db \
        --output /tmp/workdir/result/variants.vcf.gz
```

---

## Workflow for Real Data

To call SVs on your own data, follow these three stages: Preparing prefab data, generating tiles, and calling variants.

### I. Required Input Files

| File Type | Description | Requirements |
| --- | --- | --- |
| **Alignments (.bam)** | Indexed long-read alignments | Generated with `minimap2`; must have DNA sequences and quality scores. |
| **Reference (.fa)** | Reference genome | Indexed with `samtools faidx`. |
| **Matrices (.mat)** | Error models for assembly | Included in the repository under `data/lamassemble-mats/`. |
| **Annotations** | TRF and Mononucleotides | Download from [svirlpool-data](https://github.com/bihealth/svirlpool-data). |

### II. Step-by-Step Execution

Below, you will need to add the `docker run ...` or `pixi run` before the commands as explained above.

#### 1. Setup Environment

```bash
export REFERENCE=hs37d5.fa
export THREADS=16
export TRF=$DATADIR/pbsv-annotations/human_hs37d5.trf.bed
export MAT=$SVIRLPOOLDIR/data/lamassemble-mats/promethion.mat
export MNNTS=$DATADIR/HG19/hs37d5.mononucleotides.lt6.bed.gz
```

#### 2. Generate Svirltiles (Per Sample)

Run for each sample in your study (e.g., HG002 and HG003):

```bash
svirlpool run \
    --samplename HG002 \
    --workdir HG002 \
    --alignments HG002.bam \
    --reference $REFERENCE \
    --trf $TRF \
    --lamassemble-mat $MAT \
    --mononucleotides $MNNTS \
    --threads $THREADS \
    --min-sv-size 30
```

#### 3. Joint SV Calling

Combine multiple tiles into a single VCF:

```bash
svirlpool sv-calling \
    --input HG002/svirltile.db HG003/svirltile.db \
    --reference $REFERENCE \
    --output family.vcf.gz \
    --sv-types DEL INS \
    --min-sv-size 50
```

---

## Developer & Formatting Hints

For contributors or users running from source via `pixi`:

* **IDE:** Open VS Code with the pre-configured environment: `pixi run -e dev code .`
* **Format Code:** `make fix`
* **Lint & Tests:** `make check` or `make test`
* **Full Suite:** `make fix check test`

---

## Release Management

For maintainers creating releases:

* **Automated Releases:** The project uses [Release Please](https://github.com/googleapis/release-please) for automated releases. Simply merge commits with conventional commit messages to `main`, and Release Please will create/update a release PR.
* **Manual Releases:** If you need to create a manual release, see [RELEASE_GUIDE.md](RELEASE_GUIDE.md) for detailed instructions.
* **Quick Reference:** See [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for a quick checklist when creating manual releases.
