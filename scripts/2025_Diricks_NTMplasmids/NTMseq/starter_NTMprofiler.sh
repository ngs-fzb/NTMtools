#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#"This script starts NTMprofiler (Author: Jody Phelan) to predict NTM species and resistance"

###Required packages###
# NTMprofiler (https://github.com/jodyphelan/NTM-Profiler)

###Required parameters - via command line###
#-i PATH_input=""
#-o PATH_output=""
#-n filetype=""

###Optional parameters that can be changed###
Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz
Rv="_R2" #Change if your fastQ files are not SampleName_R2.fastq.gz
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') #Uses all available threads
set="SampleSet"
conda_env="NTMprofiler"
platform="illumina"

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "-i     Full path to folder where input files (fastQ or fastA) are stored"
   echo "-o     Full path to folder where result files need to be stored"
   echo "-n     Input file type: .fastq.gz or .fasta"
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: All"
   echo "-e     Name of conda environment; Default: NTMprofiler"
   echo "-f     Forward read notation; Default: _R1"
   echo "-r     Reverse read notation; Default: _R2"
   echo "-s     Name of sample set - used for file naming; Default: Sampleset"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:n:c:f:r:s:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_input=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
         PATH_output=$PATH_output_tmp/NTMprofiler ;;
      n) # 
         filetype=$OPTARG;;
      c) # 
         cpu=$OPTARG;;
      f) # 
         Fw=$OPTARG;;
      r) # 
         Rv=$OPTARG;;
      s) # 
         set=$OPTARG;;
      e) # 
         conda_env=$OPTARG;;
      v) # display Version
         echo $version
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

###Check if required parameters are provided###
if [[ -z "$PATH_input" ]] || [[ -z "$PATH_output" ]] || [[ -z "$filetype" ]]
then
	echo "Please provide all required arguments (-i PATH_input, -o PATH_output and -d db_subspecies_Mab)! use starter_NTMprofiler.sh -h for help on syntax"
	exit
fi

#No changes required below this point
###############################################################################CODE#################################################################################
###Remove previous files
#rm $PATH_output/Failed.txt

###Create folders###
mkdir -p $PATH_output

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "NTMprofiler version: "$(ntm-profiler profile --version) >> $PATH_output/info.txt
echo "Input files: "$PATH_input >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt

cd $PATH_output

if [[ "$filetype" == ".fastq" ]] || [[ "$filetype" == ".fastq.gz" ]]
then
	echo "Amount of samples in input folder:" $(ls $PATH_input/*1$filetype | wc -l) >> $PATH_output/info.txt
	for fastq in $PATH_input/*$Fw$filetype
	do
		SampleName=$(basename $fastq| cut -d '_' -f 1)
		if [[ ! -s $SampleName".results.txt" ]]
		then
			ntm-profiler profile --read1 $fastq --read2 $(echo $fastq | sed "s/$Fw$filetype$/$Rv$filetype/1") --platform $platform --dir $PATH_output --threads $cpu --csv --txt -p $SampleName
			#Clean up big files (GB file sizes)
		else
			echo "Sample was already analyzed"
		fi
	done
elif [[ "$filetype" == ".fasta" ]]
then
	echo "Amount of samples in input folder:" $(ls $PATH_input/*$filetype | wc -l) >> $PATH_output/info.txt
	for fasta in $PATH_input/*$filetype 
	do
		SampleName=$(basename $fastq | cut -d '_' -f 1 | cut -d '.' -f 1 )
		if [[ ! -s $SampleName".results.txt" ]] 
		then 
			ntm-profiler profile -a $fasta --dir $PATH_output --threads $cpu --csv --txt -p $SampleName 
		else
			echo "Sample was already analyzed"
		fi
	done
fi

###Summarize results
ntm-profiler collate



date >> $PATH_output/info.txt

#Clean up big files (GB file sizes) + small files (do it only at this point because I don´t know what collate uses)
rm $PATH_output/*.bam
rm $PATH_output/*-*-*-*.kmers.txt #BIG!
rm $PATH_output/*.fq.gz #BIG!
rm $PATH_output/*.mash_dist.txt
#rm $PATH_output/*.json Don´t do this, otherwise you cannot use the collate function again (important if you process samples in different batches)
rm $PATH_output/*.results.csv

conda deactivate
echo "Script finished!"
exit

####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################




#####################################################################HELP##########################################################################################


# usage: ntm-profiler profile [-h] [--read1 READ1] [--read2 READ2] [--bam BAM]
                            # [--fasta FASTA] [--vcf VCF]
                            # [--platform {illumina,nanopore}]
                            # [--resistance_db RESISTANCE_DB]
                            # [--external_resistance_db EXTERNAL_RESISTANCE_DB]
                            # [--species_db SPECIES_DB]
                            # [--external_species_db EXTERNAL_SPECIES_DB]
                            # [--prefix PREFIX] [--dir DIR] [--csv] [--txt]
                            # [--add_columns ADD_COLUMNS] [--call_whole_genome]
                            # [--mapper {bwa,minimap2,bowtie2,bwa-mem2}]
                            # [--caller {bcftools,gatk,freebayes}]
                            # [--calling_params CALLING_PARAMS]
                            # [--min_depth MIN_DEPTH] [--af AF]
                            # [--reporting_af REPORTING_AF]
                            # [--coverage_fraction_threshold COVERAGE_FRACTION_THRESHOLD]
                            # [--missing_cov_threshold MISSING_COV_THRESHOLD]
                            # [--species_only] [--no_trim] [--no_flagstat]
                            # [--no_clip] [--no_delly] [--no_mash]
                            # [--threads THREADS] [--verbose {0,1,2}]
                            # [--version] [--no_cleanup]

# optional arguments:
  # -h, --help            show this help message and exit

# Input options:
  # --read1 READ1, -1 READ1
                        # First read file (default: None)
  # --read2 READ2, -2 READ2
                        # Second read file (default: None)
  # --bam BAM, -a BAM     BAM file. Make sure it has been generated using the
                        # H37Rv genome (GCA_000195955.2) (default: None)
  # --fasta FASTA, -f FASTA
                        # Fasta file (default: None)
  # --vcf VCF             VCF file (default: None)
  # --platform {illumina,nanopore}, -m {illumina,nanopore}
                        # NGS Platform used to generate data (default: illumina)
  # --resistance_db RESISTANCE_DB
                        # Mutation panel name (default: None)
  # --external_resistance_db EXTERNAL_RESISTANCE_DB
                        # Path to db files prefix (overrides "--db" parameter)
                        # (default: None)
  # --species_db SPECIES_DB
                        # Mutation panel name (default: ntmdb)
  # --external_species_db EXTERNAL_SPECIES_DB
                        # Path to db files prefix (overrides "--db" parameter)
                        # (default: None)

# Output options:
  # --prefix PREFIX, -p PREFIX
                        # Sample prefix for all results generated (default:
                        # ntmprofiler)
  # --dir DIR, -d DIR     Storage directory (default: .)
  # --csv                 Add CSV output (default: False)
  # --txt                 Add text output (default: False)
  # --add_columns ADD_COLUMNS
                        # Add additional columns found in the mutation database
                        # to the text and csv results (default: None)
  # --call_whole_genome   Call whole genome (default: False)


# Algorithm options:
  # --mapper {bwa,minimap2,bowtie2,bwa-mem2}
                        # Mapping tool to use. If you are using nanopore data it
                        # will default to minimap2 (default: bwa)
  # --caller {bcftools,gatk,freebayes}
                        # Variant calling tool to use. (default: freebayes)
  # --calling_params CALLING_PARAMS
                        # Override default parameters for variant calling
                        # (default: None)
  # --min_depth MIN_DEPTH
                        # Minimum depth required to call variant. Bases with
                        # depth below this cutoff will be marked as missing
                        # (default: 10)
  # --af AF               Minimum allele frequency to call variants (default:
                        # 0.1)
  # --reporting_af REPORTING_AF
                        # Minimum allele frequency to use variants for
                        # prediction (default: 0.1)
  # --coverage_fraction_threshold COVERAGE_FRACTION_THRESHOLD
                        # Cutoff used to calculate fraction of region covered by
                        # <= this value (default: 0)
  # --missing_cov_threshold MISSING_COV_THRESHOLD
                        # Cutoff used to positions/codons in genes which are
                        # missing (this argument has now been merged with
                        # --min_depth argument and will be deprecated in future
                        # releases) (default: None)
  # --species_only        Predict species and quit (default: False)
  # --no_trim             Don't trim files using trimmomatic (default: False)
  # --no_flagstat         Don't collect flagstats (default: False)
  # --no_clip             Don't clip reads (default: True)
  # --no_delly            Don't run delly (default: False)
  # --no_mash             Don't run mash if kmers speciation fails (default:
                        # False)
  # --threads THREADS, -t THREADS
                        # Threads to use (default: 1)

# Other options:
  # --verbose {0,1,2}, -v {0,1,2}
                        # Verbosity increases from 0 to 2 (default: 0)
  # --version             show program's version number and exit
  # --no_cleanup          Don't remove temporary files on error (default: False)

echo "Script Finished!"
exit
