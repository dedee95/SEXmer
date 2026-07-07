# SEXmer detail algorithm
SEXmer was initially inspired by [KMC](https://github.com/refresh-bio/KMC) for its disk-oriented architecture and [quarTeT](https://github.com/aaranyue/quarTeT) for its dispatcher design. This tool is developed as a memory-efficient, disk-based, and [streaming computational framework](https://en.wikipedia.org/wiki/Streaming_algorithm) to minimize RAM consumption when processing large-scale kmer. Several components of SEXmer utilize [2-bit DNA sequence encoding](https://medium.com/analytics-vidhya/bioinformatics-2-bit-encoding-for-dna-sequences-9b93636e90e2) to accelerate operation and reduce memory consumption. SEXmer did not rely on WGS read alignment, as in the GWAS/SNP method, resulting in greater resource efficiency. Moreover, SEXmer adopts a modular architecture in which each module can running independenly. Therefore, users can execute only the required function without installing all dependencies. Here are the details of the computational algorithm and implementation of each SEXmer module.

![SEXmer architecture](./SEXmer_architecture.png)

## SEXmer dump
## SEXmer scan
## SEXmer reads
## SEXmer map
## SEXmer assign