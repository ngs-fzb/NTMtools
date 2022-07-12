#!/bin/bash

###Author###
echo "hello World, this script was written by Margo Diricks!"

###Function###
echo "This script calculates mash distances and creates a phylogenetic tree out of it"

###Required packages###
echo "You need to install mashtree: https://github.com/lskatz/mashtree, e.g. using conda"
	#
	
###Usage###
# sbatch file_name.sh

###!!!PARAMATERS THAT NEED TO BE CHANGED!!!###

###Input and output files###
FileType=".fasta" # Following filetypes can be used:_R1.fastq(.gz) .fasta .gbk .msh
PATH_input="" #Path to input files
PATH_output="" #Path where output files will be stored

###OPTIONAL PARAMATERS THAT YOU CAN CHANGE###

TreeName="AwesomeTree_S"$SketchSize"_k"$kmerLength"_D"$MinDepth
cpu=1 #Resources to be used (max. amount depends on your hardware configuration)
TempDir=$HOME/tmp #Path to directory where temporary files are stored (change if you have space disk issues)
MinDepth=5 #Default 5, D; Use 0 if you want to do bootstrap or jacknife
SketchSize=10000 #Default, 10000 S 
kmerLength=21 #Default, 21 k
reps=100 #Default, 100; For confidence values

###PARAMATERS THAT SHOULD NOT BE CHANGED!###

GenomeSize=5100000

###EstimatedTime###
echo "You should have your results before tomorrow!"
#Without confidence values: 1800 samples in a couple of hours


###############################################################################CODE#################################################################################

mkdir $TempDir
mkdir $PATH_output
mkdir -p $PATH_output/$TreeName
SampleCount=$(ls -l $PATH_input/*_R1.$FileType | wc -l)
echo $(mashtree --version) > $PATH_output/Info_$TreeName.txt
echo "SampleCount:"$SampleCount >> $PATH_output/Info_$TreeName.txt
echo "Input:"$PATH_input >> $PATH_output/Info_$TreeName.txt
echo "Treename:"$TreeName >> $PATH_output/Info_$TreeName.txt
echo "GenomeSize:"$GenomeSize >> $PATH_output/Info_$TreeName.txt
echo "CPU:"$cpu >> $PATH_output/Info_$TreeName.txt
echo "kmer length:"$kmerLength >> $PATH_output/Info_$TreeName.txt
echo "Sketch Size:"$SketchSize >> $PATH_output/Info_$TreeName.txt
echo "reps:"$reps >> $PATH_output/Info_$TreeName.txt

mashtree --genomesize $GenomeSize --mindepth $MinDepth --kmerlength $kmerLength --sketch-size $SketchSize --tempdir $TempDir --outmatrix $PATH_output/$TreeName/distance --numcpus $cpu $PATH_input/*$FileType > $PATH_output/$TreeName"_mash.dnd"
mv $TempDir/* $PATH_output/$TreeName

echo "Script is finished!"
echo "Tip: open the tree (.dnd file) in FigTree and export from there as nexus or Newick"

####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################

#Adding confidence values
#Jackknife
#mashtree_jackknife.pl --reps $reps --tempdir $TempDir --outmatrix $PATH_output/$TreeName/distance --numcpus $cpu $PATH_input/*$FileType -- --mindepth $MinDepth --kmerlength $kmerLength --sketch-size $SketchSize > $PATH_output/$TreeName"_mash_jackknife.dnd"

#Adding confidence values
#bootstrap
#mashtree_bootstrap.pl --reps $reps --tempdir $TempDir --outmatrix $PATH_output/$TreeName/distance --numcpus $cpu $PATH_input/*$FileType -- --mindepth $MinDepth --kmerlength $kmerLength --sketch-size $SketchSize > $PATH_output/$TreeName"_mash_bootstrap.dnd"


#####################################################################HELP#########################################################
# Usage: mashtree [options] *.fastq *.fasta *.gbk *.msh > tree.dnd
# NOTE: fastq files are read as raw reads;
      # fasta, gbk, and embl files are read as assemblies;
      # Input files can be gzipped.
# --tempdir            ''   If specified, this directory will not be
                          # removed at the end of the script and can
                          # be used to cache results for future
                          # analyses.
                          # If not specified, a dir will be made for you
                          # and then deleted at the end of this script.
# --numcpus            1    This script uses Perl threads.
# --outmatrix          ''   If specified, will write a distance matrix
                          # in tab-delimited format
# --file-of-files           If specified, mashtree will try to read
                          # filenames from each input file. The file of
                          # files format is one filename per line. This
                          # file of files cannot be compressed.
# --outtree                 If specified, the tree will be written to
                          # this file and not to stdout. Log messages
                          # will still go to stderr.
# --version                 Display the version and exit

# TREE OPTIONS
# --truncLength        250  How many characters to keep in a filename
# --sort-order         ABC  For neighbor-joining, the sort order can
                          # make a difference. Options include:
                          # ABC (alphabetical), random, input-order

# MASH SKETCH OPTIONS
# --genomesize         5000000
# --mindepth           5    If mindepth is zero, then it will be
                          # chosen in a smart but slower method,
                          # to discard lower-abundance kmers.
# --kmerlength         21
# --sketch-size        10000