#! /bin/bash

###Usage###
# bash $PathToScript/NameScript.sh

###Author###
echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Function###
echo "Running this script will download fastQ files from SRA, very fast!"

###Required packages###
echo "You need to install fasterq-dump (included in SRA-tools) and pigz, e.g. using conda" #https://anaconda.org/bioconda/sra-tools
echo "Don´t forget to activate the conda environment where you have installed these tools before running this script!"
echo "Don´t forget to change the required parameters in this script!"

###!!!PARAMATERS THAT NEED TO BE CHANGED!!!###

###Output folder###
PATH_output="$HOME/auto/Thecus_SeqData/Hinfluenzae_mdiricks/SRA_samples/Test8" #Directory where fastQ files will be stored
TempDir="$HOME/auto/Thecus_Homes/Thecus_mdiricks/Tmp" #It is important that this temporary directory has enough space, otherwise the tool wil fail!! This temporary directory will use approximately up to 10 times the size of the final output-file. this directory and its content will be deleted.

###Required FastQ files###
#Enter the accession numbers of the datasets you want to download (SRR or ERR numbers). Each accession number should be on a seperate line (you can e.g. easily copy paste a column from excel here)
group_toDownload=(
SRR12805684
SRR12805671
)

###OPTIONAL PARAMATERS THAT YOU CAN CHANGE###

cpu=16 #Resources to be used (max. amount depends on your hardware configuration)

###EstimatedTime###
#This script takes less than 2 min for 2 Haemophilus samples with compressed file sizes of about 100 MB (SRR12805684 and SRR12805671) using 16 threads

#No changes required below this point
###############################################################################CODE#################################################################################
###Pre-processing###
mkdir $PATH_output
rm $PATH_output/Info.txt

group_toDownload2=()
#Check if sample was already downloaded
for id in "${group_toDownload[@]}"
do 
	if [[ -s $PATH_output/$id"_1.fastq.gz" ]] || [[ -s $PATH_output/$id"_1.fastq" ]] || [[ -s $PATH_output/$id".fastq" ]]
	then
		echo "Sample was already downloaded"
	else 
		group_toDownload2+=( $id )
	fi
done

###Create info file###
i=0
SampleCount=${#group_toDownload2[@]}
echo $SampleCount
echo "SampleCount="$SampleCount > $PATH_output/Info.txt
echo "FirstSample="${group_toDownload2[0]} >> $PATH_output/Info.txt
echo "LastSample="${group_toDownload2[-1]} >> $PATH_output/Info.txt

###Actual command###
for id in "${group_toDownload2[@]}";
do echo "Processing $id"; 
	fasterq-dump --outdir $PATH_output --temp $TempDir --threads $cpu --split-files $id 
done

###Post-processing###
cache-mgr.2.8.2 -c #This is not required if you configure SRA tools in such a way that no caching occurs!

rm $PATH_output/Failed.txt
rm $PATH_output/downloaded.txt
for id in "${group_toDownload[@]}"
do
	if [[ -s $PATH_output/$id"_4.fastq" ]]
	then
		echo $id "	was downloaded as more than 2 files" >> $PATH_output/Failed.txt
	elif [[ -s $PATH_output/$id"_1.fastq" ]] && [[ -s $PATH_output/$id"_2.fastq" ]] && [[ -s $PATH_output/$id".fastq" ]]
	then
		echo $id "	was downloaded as paired-end files and unpaired file" >> $PATH_output/downloaded.txt
	elif [[ -s $PATH_output/$id"_1.fastq" ]] && [[ -s $PATH_output/$id"_2.fastq" ]]
	then
		echo $id "	was downloaded as 2 paired-end files" >> $PATH_output/downloaded.txt
	elif [[ -s $PATH_output/$id".fastq" ]]
	then
		echo $id "	was downloaded as 1 file" >> $PATH_output/Failed.txt
	else 
		echo $id >> $PATH_output/Failed.txt
	fi
done

####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################

cd $PATH_output
#mmv "*_*.fastq" "#1_SRA_Download_Date_xbp_R#2.fastq"
#gzip *.fastq # much slower!
pigz *.fastq

echo "Script Finished!"


#####################################################################HELP#########################################################

# Usage: fasterq-dump [ options ] [ accessions(s)... ]

# Parameters:

  # accessions(s)                    list of accessions to process


# Options:

  # -o|--outfile <path>              full path of outputfile (overrides usage
                                     # of current directory and given accession)
  # -O|--outdir <path>               path for outputfile (overrides usage of
                                     # current directory, but uses given
                                     # accession)
  # -b|--bufsize <size>              size of file-buffer (dflt=1MB, takes
                                     # number or number and unit where unit is
                                     # one of (K|M|G) case-insensitive)
  # -c|--curcache <size>             size of cursor-cache (dflt=10MB, takes
                                     # number or number and unit where unit is
                                     # one of (K|M|G) case-insensitive)
  # -m|--mem <size>                  memory limit for sorting (dflt=100MB,
                                     # takes number or number and unit where
                                     # unit is one of (K|M|G) case-insensitive)
  # -t|--temp <path>                 path to directory for temp. files
                                     # (dflt=current dir.)
  # -e|--threads <count>             how many threads to use (dflt=6)
  # -p|--progress                    show progress (not possible if stdout used)
  # -x|--details                     print details of all options selected
  # -s|--split-spot                  split spots into reads
  # -S|--split-files                 write reads into different files
  # -3|--split-3                     writes single reads into special file
     # --concatenate-reads           writes whole spots into one file
  # -Z|--stdout                      print output to stdout
  # -f|--force                       force overwrite of existing file(s)
  # -N|--rowid-as-name               use rowid as name (avoids using the name
                                     # column)
     # --skip-technical              skip technical reads
     # --include-technical           explicitly include technical reads
  # -P|--print-read-nr               include read-number in defline
  # -M|--min-read-len <count>        filter by sequence-lenght
     # --table <name>                which seq-table to use in case of pacbio
     # --strict                      terminate on invalid read
  # -B|--bases <bases>               filter output by matching against given
                                     # bases
  # -A|--append                      append to output-file, instead of
                                     # overwriting it
     # --ngc <path>                  <path> to ngc file
     # --perm <path>                 <path> to permission file
     # --location <location>         location in cloud
     # --cart <path>                 <path> to cart file
  # -V|--version                     Display the version of the program
  # -v|--verbose                     Increase the verbosity of the program
                                     # status messages. Use multiple times for
                                     # more verbosity.
  # -L|--log-level <level>           Logging level as number or enum string.
                                     # One of
                                     # (fatal|sys|int|err|warn|info|debug) or
                                     # (0-6) Current/default is warn
     # --option-file file            Read more options and parameters from the
                                     # file.
  # -h|--help                        print this message

# fasterq-dump version 2.11.0

exit
