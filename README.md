# NTMtools
This is a repository with scripts and manuals for analysing NGS data from nontuberculous mycobacteria

## Scripts ##
### Assemble_usingShovill_Mabs.sh ###
#### Function ####
Wrapper script for creating assemblies from illumina fastQ files for Mycobacterium abscessus. <br /> 
This bash script calls the shovill pipeline (https://github.com/tseemann/shovill) and summarises the output (assembly statistics) into a text file.
#### Packages to be installed (e.g. with conda) ####
1) shovill
#### Usage ####
_activate conda environment where shovill is installed_ <br /> 
$ conda activate EnvName <br /> 
_run script_ <br /> 
$ bash Assemble_usingShovill_Mabs.sh <br /> 

Input: FastQ files (illumina paired end reads) <br /> 
Output: FastA Assemblies (in folder "FinalAssemblies") and a text file containing assembly statistics <br /> 

Parameters that need to be adjusted: <br /> 
_PATH_input_= path where your fastQ files are stored <br /> 
_PATH_output_= path where you want to store assemblies <br /> 

Parameters that might need to be adjusted: <br /> 
_Fw_= forward read notation (default R1, according to naming convention SampleName_R1.fastq.gz ) <br /> 
_Rv_= forward read notation (default R2, according to naming convention SampleName_R2.fastq.gz ) <br /> 
_ass_= assembler (default skesa) <br /> 
_cov_= coverage to which you want to downsample (default 100) <br /> 
_cpu_= resources you want to use (default 0 = all cpus) <br /> 
Note: Additional parameters for shovill can be added or removed (see https://github.com/tseemann/shovill) <br /> 

#### Estimated time #### 
Between 4 and 20 min per sample depending on coverage, assembler and resources used <br /> 

### sra_download_parallel_MD.sh ###
#### Function ####
Download FastQ files in parallel (to reduce total time for downloading the files).
#### Packages to be installed (e.g. with conda) ####
1) fastq-dump
2) GNU parallel
#### Usage ####
_activate conda environment where these two packages are installed_ <br /> 
$ conda activate EnvName <br /> 
_run script_ <br /> 
$ bash sra_download_parallel_MD.sh <br /> 

Input: / <br /> 
Output: FastQ files (in folder "Outdir") <br /> 

Parameters that need to be adjusted: <br /> 
_Outdir_= path where your fastQ files will be stored <br /> 
_group_sra_download_= accession numbers that you want to download <br /> 
_Jobs_= Amount of fastQ files that need to be downloaded in parallel. This number is limited by your hardware configuration. Note that you also might run into full disk error as there is a ncbi folder created in /home folder which fills up space quickly if you download too much in parallel. A possible solution for this is to work with symlinks (see line 53-58;74) <br /> 

Parameters that might need to be adjusted: <br /> 
Line 73: fastq-dump.2.8.2 --> change e.g. to fastq-dump if you have downloaded this package via conda <br /> 


#### Estimated time #### 
Depends on file size and resources used <br /> 

## Manuals ##
### cgMLST_Mabs.docx ###
This manual includes detailed information to perform cgMLST analysis, including steps to identify subspecies or DCC status of newly sequenced isolates.



