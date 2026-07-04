# SEXmer: Fast and Resource Efficient Sex Determination Analysis Using Kmer
![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)
![MIT license](https://img.shields.io/badge/License-MIT-Blue.svg)

<img src="docs/SEXmer-logo.png" alt="SEXmer logo" width="300">

Identify sex determination region (SDR) in some plant and animal require huge effort, especially for XY and ZW sex type. To be able detect SDR robustly, we generally need population samples for both male and female individuals. This study oftem produces large whole genome sequencing (WGS) data. K-mer-based method is a powerful strategy for detecting sex determination region. However, processing population-scale kmer data require huge computational resources. 

Here we present SEXmer, a fast and resource efficient command line tools for sex determination region analysis based on kmer. 