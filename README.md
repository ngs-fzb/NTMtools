# NTMtools
This is a repository with scripts and manuals for analysing NGS data from nontuberculous mycobacteria

## Scripts ##
### Assemble_usingShovill_Mabs.sh ###
#### Function ####
Wrapper script for creating assemblies from illumina fastQ files for Mycobacterium abscessus. <br /> 
This bash script calls the shovill pipeline (https://github.com/tseemann/shovill) and summarises the output (assembly statistics) into a text file.

### sra_download_parallel_MD.sh ###
#### Function ####
Download FastQ files in parallel (to reduce total time for downloading the files).

### starter_Mashtree_Mab_MD.sh ###
#### Function ####
This script calculates mash distances and creates a phylogenetic tree out of it".

#### Usage ####
_Install required packages (see script)_ <br />, e.g. using conda
_activate conda environment where shovill is installed_ <br /> 
$ conda activate EnvName <br /> 
_run script_ <br /> 
$ bash $PathToScript/scriptname.sh <br /> 
 
## Manuals ##
### cgMLST_Mabs.docx ###
This manual includes detailed information to perform cgMLST analysis for M. abscessus, including steps to identify subspecies or DCC status of newly sequenced isolates.



