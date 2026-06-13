# SEXmer: Fast and Efficient Sex Determination (ZW & XY) Region Analysis Tool Using Kmer.
![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)
![MIT license](https://img.shields.io/badge/License-MIT-Blue.svg)

<img src="docs/SEXmer-logo.png" alt="findGEVE logo" width="300">

Identify sex determination region (SDR) in some plant and animal require relentless effort, especially for XY and ZW sex type. To be able detect robustly SDR in XY or ZW, we need pooled sample for both male and female. This study always generated tons of WGS data. One of robust method to detect SDR in the pooled samples is kmer based method. However, processing the WGS data will require huge computational resources. Here I present SEXmer, a fast and resource efficient command line tools for sex determination region analysis based on kmer. 