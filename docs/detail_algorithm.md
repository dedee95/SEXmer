# SEXmer detail algorithm
SEXmer was initially inspired by [KMC](https://github.com/refresh-bio/KMC) for its disk-oriented architecture and [quarTeT](https://github.com/aaranyue/quarTeT) for its dispatcher design. This tool is developed as a memory-efficient, disk-based, and [streaming computational framework](https://en.wikipedia.org/wiki/Streaming_algorithm) to minimize RAM consumption when processing large-scale kmer. Several components of SEXmer utilize [2-bit DNA sequence encoding](https://medium.com/analytics-vidhya/bioinformatics-2-bit-encoding-for-dna-sequences-9b93636e90e2) to accelerate operation and reduce memory consumption. SEXmer did not rely on WGS read alignment, as in the GWAS/SNP method, resulting in greater resource efficiency. Moreover, SEXmer adopts a modular architecture in which each module can running independenly. Therefore, users can execute only the required function without installing all dependencies. Here are the details of the computational algorithm and implementation of each SEXmer module.

![SEXmer architecture](./SEXmer_architecture.png)

## SEXmer dump
This module basically is very simple; it generates kmer frequency (dump) from each raw WGS read. The main backbone of this module is KMC. We chose KMC over other kmer counting tools like Jellyfish or Meryl because KMC uses a disk-oriented architecture to perform kmer counting. In short, KMC uses less RAM, and it operates very fast. 
To reduce output data, we implemented a minimum kmer count and trigger sequence. Kmer with a count less than 3 is most likely an artifact or sequencing error, so better remove it. On the other hand, the trigger sequence makes the SEXmer dump only keep kmer with the start nucleotide specified in the trigger sequence. For example, if we specify the trigger sequence as "AG", only retained kmer with start "AG" will be kept. This can significantly reduce kmer output (1/16) while still representing the important information. 

## SEXmer scan
Lorem ipsum sir dolor sit amet.

## SEXmer reads
SEXmer reads use bbduk.sh from the BBMAP package to extract reads using kmer. We chose to implement bbduk over other tools because bbduk runs very fast and uses less RAM. Not just for WGS short paired reads, bbduk can also be used to extract specific reads from long reads sequencing. This module works by matching the given kmer sequence to reads. Whoever reads hits into the given kmer, bbduk will keep it. For WGS paired-end, bbduk writes both mates if either mate has at least the minimum marker hits. The minimum hits for WGS short reads is 1, while short read sequencing should have a minimum of more than 1. Based on our experience, at least 3 hits should be sufficient for long-read extraction.

## SEXmer map
The core function of this module is to identify genomic regions enriched with a specific kmer. In simple terms, the SEXmer map asks, "Where are sex-specific kmer sequences located and enriched in the genome?" SEXmer takes a sex-specific kmer sequence (MSK or FSK), found where they are enriched in the genome by using a rapid kmer-based genome scanning approach. This method is much faster than traditional alignment like BWA-MEM2. Moreover, there is no need to split the chromosome into several chunks if the size is too big, like in BWA-MEM2. Here are the details on how this algorithm works


```
1. Converted a sex-specific kmer sequence into a 2-bit encoded integer and stored it in a hash-based marker index.
2. Scan the reference genome sequentially using a rolling integer kmer scanning algorithm. 
3. Generate each genomic kmer at the current position using efficient 2-bit operations (bit shifting, masking, and integer encoding). 
4. Query each genomic kmer against the marker index. 5. If a match is detected, record the chromosome and genomic coordinate of the kmer hit. 
6. Assign each detected kmer hit to overlapping sliding windows based on its genomic coordinate. 
7. Count the number of marker hits within each window and normalize the hit density based on the number of valid genomic kmer sites. 
8. Generate a genome-wide kmer enrichment profile for visualization and identification of candidate sex-determining regions.
```

SEXmer map also provides sex-specific reads mapping to give more evidence of the sex determination region. This module asks: "Do actual sequencing reads support the genomic location identified by sex-specific kmers?" The main computational engine for this step is BBMap, which can map both short-read and long-read sequencing data. The main actor of this part is bbmap. It can map both WGS short reads and long reads. But bbmap has a limit for mapping large genomes, so SEXmer map splits the large genome into several chunks. Also, bbmap split the long reads into 6000 bp pieces. So, the reads mapping here is just to facilitate double validation of SDR. If comprehensive mapping is needed, BWA-MEM2 or minimap2 is a better option.


## SEXmer assign