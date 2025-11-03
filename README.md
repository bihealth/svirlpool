# Svirlpool: structural variant detection from long read sequencing by local assembly

<img src="media/illustration.png" align="right" alt="Svirlpool illustration" width="300" height="300">

- License: MIT
- Contact: [Vinzenz May](mailto:vinzenz.may@bih-charite.de)

**Motivation**

Long-Read Sequencing (LRS) promises great improvements in the detection of structural genome variants (SVs). However, existing methods are lacking in key areas such as the reliable detection of inserted sequence, precise genotyping of variants, and reproducible calling of variants across multiple samples. Here, we present our method Svirlpool, that is aimed at the analysis of Oxford Nanopore Technologies (ONT) sequencing data. Svirlpool uses local assembly of candidate SV regions to obtain high-quality consensus sequences.

**Results**

Svirlpool obtains competitive results to the leading method Sniffles on the widely used Genome in a Bottle (GiaB) benchmark data sets. On trio data, however, Svirlpools shows a clear favorable performance in terms of mendelian consistency. This indicates that Svirlpool shows great promise in clinical applications and beyond benchmark datasets.

## Installing & Running

Currently, the only way to run Svirlpool is in developer mode - so see "Developer Notes" below.

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

TBD: run example

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
