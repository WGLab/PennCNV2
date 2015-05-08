# PennCNV2

PennCNV2 is the second major release of the popular program PennCNV, which was originally designed for inferring copy number variation in germline DNA using data from genotyping microarrays. PennCNV2 now supports inference of copy number abberrations in tumor cells (PennCNV-tumor). New functionality for inference of copy number variation using next gen sequencing data is also developed under the same framework (PennCNV-Seq).

In PennCNV2, we estimate stromal contamination by applying a maximum likelihood approach over multiple discrete genomic intervals. By conditioning on signal intensity across the genome, our method accounts for both aneuploidy and genomic waves. Finally, our method uses a hidden Markov model to integrate multiple sources of information, including total and allele-specific signal intensity at each SNP, as well as physical maps to make posterior inferences of CNAs. 

The software is written in C++ and targeted for unix/linux platforms.

## Reference:

- Chen GK, Chang X, Curtis C, Wang K. [Precise inference of copy number alterations in tumor samples from SNP arrays](http://bioinformatics.oxfordjournals.org/content/29/23/2964.long). _**Bioinformatics**_ 29(23):2964-70, 2013

