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

## Manuals ##
### cgMLST_Mabs.docx ###
This manual includes detailed information to perform cgMLST analysis in SeqSphere+ (https://www.ridom.de/seqsphere/), including steps to identify subspecies or DCC status of newly sequenced isolates.



