# SEXmer: Fast and Resource Efficient Sex Determination Analysis Using Kmer
![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)
![MIT license](https://img.shields.io/badge/License-MIT-Blue.svg)

<img src="docs/SEXmer-logo.png" alt="SEXmer logo" width="300">

Identifying the sex determination region (SDR) in some plants or animals requires huge effort, especially for XY and ZW sex types. To detect SDR robustly, we generally need population samples for both male and female individuals. This study often produces large whole-genome sequencing (WGS) data. K-mer-based method is a powerful strategy for detecting the sex determination region. However, processing population-scale kmer data requires huge computational resources.

Here we present `SEXmer`, a fast and resource-efficient command-line tool for sex determination region analysis based on kmer. 

> **"`SEXmer` provides a modular workflow from raw reads to sex-specific k-mer discovery, kmer based reads extraction, unknown sex classifier, and genomic localization of candidate SDR signals."**

Currently, SEXmer contain 5 modules:

```bash
dump      Generate filtered canonical k-mer dump files.
scan      Identify sex-specific and sex-biased k-mers.
reads     Extract reads containing sex-specific k-mers.
map       Map sex-specific k-mers or sex-specific reads.
assign    Assign sex using validated sex-specific markers.
```

## Getting Started
SEXmer command line tool currently only available for Linux. SEXmer is implemented as a bash script, embedded with Python, and using some external dependencies.

**Dependencies**
- [Python3](https://www.python.org/) (>3.8, tested on 3.10)
- Python library: 
  - numpy >=1.26 
  - pandas >=2.2
  - matplotlib >=3.8
  - scipy >=1.11
- [KMC](https://github.com/refresh-bio/KMC) (tested on 3.2.4)
- [BBMAP](https://bbmap.org/) (tested on 39.81)

All of these dependencies can easily be installed using conda.
```bash
#create the environment with it depedencies
conda create -n sexmer -c conda-forge -c bioconda python=3.10 kmc bbmap numpy pandas matplotlib scipy

#activate the environment
conda activate sexmer
```

**Instalation**

Currently, `SEXmer` only supports manual installation. Clone the repository or download a specific released package.
```bash
git clone https://github.com/dedee95/SEXmer.git
cd SEXmer
chmod -R 755 bin/
export PATH="$PWD/bin:$PATH"
```

## Quick Usage Guide
After all dependencies are installed, type `SEXmer -h` to verify installation.
```bash
SEXmer: Resource-efficient toolkit for sex determination region analysis based on k-mers

Usage: SEXmer <module> <parameters>

Modules:
dump        Generate filtered canonical k-mer dump files.
scan        Identify sex-specific and sex-biased k-mers.
reads       Extract reads containing sex-specific k-mers.
map         Map sex-specific k-mers or sex-specific reads.
assign      Assign sex using validated sex-specific markers.

Use <module> -h for module usage.
```

## SEXmer Detail Algoritm
