# NTMtools
This repository contains scripts and manuals for whole genome sequencing based analysis of non-tuberculous mycobacteria. Questions? Ask mdiricks@fz-borstel.de

The wrapper script starter_NTMseq.sh runs on linux and accepts FastQ (illumina) and FastA files as input.
It includes:

1.	Quality control of raw sequence reads
[Input: FastQ files; Required tools: FastQC and multiQC]
2.	Contamination detection 
[Input: FastQ files or FastA files; Required tools: kraken2]
3.	Preprocessing of raw sequence reads (e.g. adapter removal) 
[Required tools: fastp]
4.	Mycobacterium (sub)species and resistance prediction 
[Input: FastQ files; Required tools: NTMprofiler]
5.	Multi-locus sequence typing (MLST)
[Input: FastQ files; Required tools: srst2; Required database: pubMLST (only available for M. abscessus)]
6.	Assembly of sequence reads
[Input: FastQ files; Required tools: Shovill; Output: FastA files]
7.	Fast phylogenetic analysis using Mashtree
[Input: FastQ files or FastA files; Required tools: Mashtree]
8.	Detection of known plasmids
[Input: FastQ files; Required tools: srst2 and seqkit; Required database: PLSDB or custom]
9.	De novo prediction of plasmid contigs
[Input: FastQ files; Required tools: platon and/or plasmidspades]
10.	Resistance and virulence gene prediction
[Input: FastA files; Required tools: AMRfinder+]

# Getting Started
For complete installation instructions of NTMseq, description and usage examples please send a mail to mdiricks@fz-borstel.de.

# Citation
https://github.com/ngs-fzb/NTMtools/tree/main/scripts/NTMseq or https://doi.org/10.1186/s12866-025-04563-7





