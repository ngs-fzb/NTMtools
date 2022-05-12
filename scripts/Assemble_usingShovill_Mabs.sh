#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo Diricks!"

###Function###
#Wrapper script for creating assemblies from illumina fastQ files for M. abscessus
#This bash script calls the shovill pipeline (https://github.com/tseemann/shovill) and summarises the output (assembly statistics) into a text file.

###Required packages###
#https://github.com/tseemann/shovill
	
###Usage###
# bash file_name.sh

###Input and output files###
#Input Type: Paired-end raw fastQ files (with sufficient coverage!)
PATH_input=""
#Change if your fastQ files are not SampleName_R1.fastq.gz and SampleName_R2.fastq.gz
Fw="R1"
Rv="R2"
#Output Type: FastA assemblies
PATH_output=""

###Shovill Parameters###
#Expected genome size
egs="5.1M"
#Downsampling to
cov=100
#Assembler (choose skesa velvet megahit or spades)
ass="skesa"
#Resources
cpu=0

###Estimated time### 
#Between 4 and 20 min per sample depending on coverage and assembler used

################################################################################################################################################################
mkdir -p $PATH_output
mkdir $PATH_output/FinalAssemblies
rm $PATH_output/FinalAssemblies/Failed.txt
#Loop through all the raw fastQ files
for fastq in $PATH_input/*_$Fw.fastq.gz 
do
	SampleName=$(basename $fastq| cut -d '_' -f 1)
	if [ -f "$PATH_output/FinalAssemblies/$SampleName.fasta" ]
	then
		 echo "Sample was already analysed" >> $PATH_output/FinalAssemblies/AlreadyAnalyzed.txt
	else
		shovill --R1 $fastq --R2 $(echo $fastq | sed "s/$Fw.fastq.gz$/$Rv.fastq.gz/1") --outdir $PATH_output/$SampleName --gsize $egs --depth $cov --trim --assembler $ass 
#Move the final assembly to a seperate folder and rename file
		if [ -f "$PATH_output/$SampleName/contigs.fa" ]
		then
			mv $PATH_output/$SampleName/contigs.fa $PATH_output/FinalAssemblies/$SampleName"_Assembly_Shovill_"$ass.fasta
		else
			echo $SampleName >> $PATH_output/FinalAssemblies/Failed.txt
		fi
	fi
done

# Create summary file
rm $PATH_output/"shovillOutput_summary.txt"
#Make header of the summary file 
echo -e "Sample\tEstimated_depth(x)\tRead_max_len\tRead_avg_len\tRead_min_len\tSurviving_Read_Pairs_trimmomatic_perc\tFw_surviving_trimmomatic_perc\tRv_surviving_trimmomatic_perc\tDropped_trimmomatic_perc\tFlash_CombinedPairs_perc\tWalltime\tContigAmount\tMinContigLength\tAssemblyLength\tAssemblyRel%("$egs")\tAssembler\tREC\tShovillVersion\tFullCode" >> $PATH_output/shovillOutput_summary.txt

#Loop through all the generated log files and extract necessary info
for log in $PATH_output/*/shovill.log
do
	#ReadStats
	SampleName=$(basename $(dirname $log))
	max_len=$(grep 'max_len' $log| cut -d '=' -f 2)
	avg_len=$(grep 'avg_len' $log| cut -d '=' -f 2)
	min_len=$(grep 'min_len' $log| cut -d '=' -f 2)
	EstSeqDepth=$(grep 'Estimated.sequencing.depth' $log| cut -d ':' -f 2 | sed 's/\s//' | sed 's/x//')
	Flash_CombinedPairs_perc=$(grep -o 'Percent.combined.*%' $log| cut -d ':' -f 2 | sed 's/\s//' | sed 's/%//')
	Surviving_Read_Pairs_trimmomatic_perc=$(grep -Po 'Both.Surviving.*?%' $log| grep -Po '.....%' | sed 's/%//') #-P means that you want to use non-greedy expression .*? which means any character multiple times until the first character you put after this expression; -o means output only matching part, not complete line
	Fw_surviving_trimmomatic_perc=$(grep -P -o 'Forward.Only.Surviving.*?%' $log| grep -Po '....%' | sed 's/%//')
	Rv_surviving_trimmomatic_perc=$(grep -P -o 'Reverse.Only.Surviving.*?%' $log| grep -Po '....%' | sed 's/%//')
	Dropped_trimmomatic_perc=$(grep -Po 'Dropped.*?%' $log| grep -Po '....%' | sed 's/%//')
	Walltime=$(grep 'Walltime' $log| cut -d ':' -f 2 | sed 's/\s//')
	Contig_amount=$(grep 'It.contains' $log| cut -d ' ' -f 4)
	min_contig_length=$(grep -P -o '\(min=.*?\)' $log| cut -d '=' -f 2 | sed -e 's/)//')
	Assembly_length=$(grep -P -o 'Assembly.is.*?,' $log| cut -d ' ' -f 3 | sed 's/,//')
	Assembly_rel=$(grep 'Assembly.is' $log| grep -Po '\(.*?\)'| sed -e 's/)//' | sed -e 's/(//' | sed -e 's/%//')
	GC=$(awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.1f", (gc*100)/(gc+at) }' $PATH_output/FinalAssemblies/$SampleName"_Assembly_Shovill_"$ass.fasta)
	Assembler=$(grep -oP -m 1 'assembler\s.*' $log | cut -d ' ' -f 2) #-P to use Perl regular expression, -o to output only match, \s to match space, \K to omit the match, .* to match rest of the string(s)
	REC="Unknown" #Read error correction option from shovill
	if grep -Fxq 'noreadcorr' $log; then REC="no"; else REC="Yes"; fi
	Shovill_v=$(grep 'This.is.shovill' $log | cut -d ' ' -f 5)
	FullCode=$(grep 'You.ran' $log)
	paste --delimiter='\t'  <(echo $SampleName)	<(echo $EstSeqDepth) <(echo $max_len) <(echo $avg_len) <(echo $min_len) <(echo $Surviving_Read_Pairs_trimmomatic_perc) <(echo $Fw_surviving_trimmomatic_perc) <(echo $Rv_surviving_trimmomatic_perc) <(echo $Dropped_trimmomatic_perc) <(echo $Flash_CombinedPairs_perc) <(echo $Walltime) <(echo $Contig_amount) <(echo $min_contig_length) <(echo $Assembly_length) <(echo $Assembly_rel) <(echo $Assembler) <(echo $REC) <(echo "v"$Shovill_v) <(echo $FullCode)>> $PATH_output/shovillOutput_summary.txt
done

exit 0
################################################################################################################################################################
#QUICK START: shovill --outdir out --R1 test/R1.fq.gz --R2 test/R2.fq.gz

#MAIN STEPS
# Estimate genome size and read length from reads (unless --gsize provided)
# Reduce FASTQ files to a sensible depth (default --depth 100)
# Trim adapters from reads (with --trim only)
# Conservatively correct sequencing errors in reads
# Pre-overlap ("stitch") paired-end reads
# Assemble with SPAdes/SKESA/Megahit with modified kmer range and PE + long SE reads
# Correct minor assembly errors by mapping reads back to contigs
# Remove contigs that are too short, too low coverage, or pure homopolymers
# Produce final FASTA with nicer names and parseable annotations

# SYNOPSIS
  # De novo assembly pipeline for Illumina paired reads
# USAGE
  # shovill [options] --outdir DIR --R1 R1.fq.gz --R2 R2.fq.gz
# GENERAL
  # --help          This help
  # --version       Print version and exit
  # --check         Check dependencies are installed
# INPUT
  # --R1 XXX        Read 1 FASTQ (default: '')
  # --R2 XXX        Read 2 FASTQ (default: '')
  # --depth N       Sub-sample --R1/--R2 to this depth. Disable with --depth 0 (default: 150)
  # --gsize XXX     Estimated genome size eg. 3.2M <blank=AUTODETECT> (default: '')
# OUTPUT
  # --outdir XXX    Output folder (default: '')
  # --force         Force overwite of existing output folder (default: OFF)
  # --minlen N      Minimum contig length <0=AUTO> (default: 0)
  # --mincov n.nn   Minimum contig coverage <0=AUTO> (default: 2)
  # --namefmt XXX   Format of contig FASTA IDs in 'printf' style (default: 'contig%05d')
  # --keepfiles     Keep intermediate files (default: OFF)
# RESOURCES
  # --tmpdir XXX    Fast temporary directory (default: '/tmp/tseemann')
  # --cpus N        Number of CPUs to use (0=ALL) (default: 8)
  # --ram n.nn      Try to keep RAM usage below this many GB (default: 16)
# ASSEMBLER
  # --assembler XXX Assembler: skesa velvet megahit spades (default: 'spades')
  # --opts XXX      Extra assembler options in quotes eg. spades: "--untrusted-contigs locus.fna" ... (default: '')
  # --kmers XXX     K-mers to use <blank=AUTO> (default: '')
# MODULES
  # --trim          Enable adaptor trimming (default: OFF)
  # --noreadcorr    Disable read error correction (default: OFF)
  # --nostitch      Disable read stitching (default: OFF)
  # --nocorr        Disable post-assembly correction (default: OFF)
  
#--trim option:
#[trimmomatic]  -threads 8 -phred33 R1.sub.fq.gz R2.sub.fq.gz R1.fq.gz /dev/null R2.fq.gz /dev/null ILLUMINACLIP:/home/mdiricks/bin/miniconda2/db/trimmomatic.fa:1:30:11 LEADING:3 TRAILING:3 MINLEN:30 TOPHRED33

# This will perform the following in this order
# Remove Illumina adapters provided in the trimmomatic.fa file (provided). Initially
# Trimmomatic will look for seed matches (16 bases) allowing maximally 1
# mismatch. These seeds will be extended and clipped if in the case of paired end
# reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a
# score of 11, (about 17 bases).
# Remove leading low quality or N bases (below quality 3)
# Remove trailing low quality or N bases (below quality 3)
# Drop reads which are less than 30 bases long after these steps