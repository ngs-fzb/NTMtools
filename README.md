# NTMtools
This is a repository with scripts and manuals for analysing NGS data from nontuberculous mycobacteria

## Scripts ##

#### Usage ####
_Install required packages (see script), e.g. using conda_ <br />
_Activate conda environment where packages are installed_ <br /> 
$ conda activate EnvName <br /> 
_run script_ <br /> 
$ bash PathToScript/scriptname.sh <br /> 

### starter_Assemble_usingShovill_Mab.sh ###
#### Function ####
Wrapper script for creating assemblies from illumina fastQ files for Mycobacterium abscessus. <br /> 
This bash script calls the shovill pipeline (https://github.com/tseemann/shovill) and summarises the output (assembly statistics) into a text file.

### starter_SRA-Download_usingFasterq-dump.sh ###
#### Function ####
Download FastQ files from SRA, very fast.

### starter_PhyloTree_usingMashtree_Mab.sh ###
#### Function ####
This script calculates mash distances and creates a phylogenetic tree out of it".

## Manuals ##
### cgMLST_Mabs.docx ###
This manual includes detailed information to perform cgMLST analysis for M. abscessus, including steps to identify subspecies or DCC status of newly sequenced isolates.



