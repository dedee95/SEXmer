# SEXmer: Fast and Resource Efficient Sex Determination Analysis Using Kmer
![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
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

## Table of Contents