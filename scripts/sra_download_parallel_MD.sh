#! /bin/bash
###Author###
echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Function###
echo "Running this script will start downloading fastQ files in parallel using a combination of fastq-dump and GNU parallel (much faster than using fastq-dump alone)!"

echo "Activate a conda environment where parallel (and fastq-dump) is installed"
echo "Note that in this script fastq-dump.2.8.2 is called (which was not installed via conda in our case). If you have another version installed or installed fastq-dump via conda, change the call name in this script (line 73)"
echo "Change Outdir in this script"

###Required packages###
#  fastq-dump and GNU parallel --> these can be installed e.g. using conda.

###Usage###
# bash $PathToScript/Name_Script.sh"

#Change parameters below!!
######################################################################################################################################################

###Output folder###
#Directory where fastQ files will be stored
Outdir=""

###Required FastQ files###
#Enter the accession numbers of the datasets you want to download (SRR or ERR numbers) starting from line 28. Each accession number should be on a seperate line (you can e.g. easily copy paste a column from excel here)
group_sra_download=(

)

###Jobs###
#Amount of fastQ files that need to be downloaded in parallel. This number is limited by your hardware configuration. Note that you also might run into full disk error as there is a ncbi folder created in /home folder which fills up space quickly if you download too much in parallel. A possible solution for this is to work with symlinks (see line 53-58;74) 
Jobs=8

#No changes required below this point, unless you need to change call name for fastq-dump (line 73)
######################################################################################################################################################
cache-mgr.2.8.2 -c
group_toDownload=()
###Check if sample was already downloaded###
for id in "${group_sra_download[@]}"
do 
	if test -f $Outdir/$id"_1.fastq.gz"
	then
		echo "Sample was already downloaded"
	else 
		group_toDownload+=( $id )
	fi
done

###Download in parallel###

TmpDi="$HOME/tmp"
# if test -L ~/ncbi
# then 
	# echo "symlink already exists"
# else 
	# ln -s $HOME/auto/Fastq/ncbi ~/ncbi #otherwise all temporary files are stored in the home folder and you will run into full disk error
# fi

mkdir $Outdir
rm $Outdir/Info.txt

i=0
SampleCount=${#group_toDownload[@]}
echo $SampleCount
echo "SampleCount="$SampleCount > $Outdir/Info.txt
echo "FirstSample="${group_toDownload[0]} >> $Outdir/Info.txt
echo "LastSample="${group_toDownload[-1]} >> $Outdir/Info.txt

source `which env_parallel.bash`
until [ $i -gt $SampleCount ]
do
	env_parallel -j $Jobs --compress --tmpdir $TmpDi 'fastq-dump.2.8.2 --outdir $Outdir --split-files --gzip {}' ::: "${group_toDownload[@]:$i:120}"
	#rm $HOME/auto/Fastq/ncbi/public/sra/*
	cache-mgr.2.8.2 -c
	i=$(($i+119))
done

rm $Outdir/Failed.txt
for id in "${group_sra_download[@]}";
do
if test -f $Outdir/$id"_1.fastq.gz"
	then
		echo "Sample was downloaded"
	else
		echo $id >> $Outdir/Failed.txt
fi
done
cd $Outdir
###If you want, change name of your fastQ file###
#mmv "*_*.fastq.gz" "#1_SRA_Download_Date_xbp_R#2.fastq.gz"