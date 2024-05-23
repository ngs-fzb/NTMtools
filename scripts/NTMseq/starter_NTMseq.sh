#! /bin/bash

###Usage###
# bash $PathToScript/Scriptname.sh

###Author###
echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Function###
echo "This is a starter script for NTMseq, a pipeline for analyis of WGS data from non-tuberculous mycobacteria"

###Version###
version="1.0.2"

#######################
###REQUIRED PACKAGES###
#######################
#Create conda environments with the names in between "".

conda_env_SRST2="SRST2" #Required packages: SRST2 (tested with v0.2.0) and GNU parallel (tested with v. 20210222)
conda_env_spades="spades" #Required packages: spades (tested with v3.15.5)
conda_env_shovill="Shovill" #Required packages: shovill (tested with v1.1.0)
conda_env_platon="platon" #Required packages: platon (tested with v1.6)
conda_env_mashtree="mashtree" #Required packages: mashtree (tested with v.1.2.0)
conda_env_NTMprofiler="NTMprofiler" #Reguired packages: ntm-profiler tested with (v.0.2.1)
conda_env_kraken2="kraken2" #Required packages: Kraken (tested with v.2.1.2) and ktImportText. Note that you also need to install an appropriate kraken2 database and krona database!
conda_env_seqkit="seqkit" #Required packages: seqkit (tested with v2.3.1) and pigz tested with (v2.6)
conda_env_multiqc="multiqc" #Required packages: fastqc (tested with v0.11.4) and multiqc (tested with v1.13.dev0)
conda_env_fastp="fastp" #Required packages: fastp (tested with v.0.23.2)
conda_env_AMRfinder="AMRFinder" #Required packages: AMRfinder (tested with v.3.11.2 and database version 2022-12-19.1)
conda_env_fastANI="fastani" #Required packages: fastani (tested with v.1.33)

#Note: The conda enviroment is activated from within the script!"
#Note: Don´t forget to download the latest databases (e.g. PLSDB, MLST, platon) !" # See below: CODE THAT MIGHT BE USED IN ADDITION

##############################################
###!!!PARAMATERS THAT NEED TO BE CHANGED!!!###
##############################################

##Input/Output##
PATH_scripts="" #Path to folder where scripts are stored
PATH_fastQ="" # Path to folder where fastQ files are stored
PATH_output="" # Path to folder where result files will be stored

##Databases##
db_platon="" # Syntax: $Path_to_folder; Don´t forget to download the latest databases: https://zenodo.org/record/4066768/files/db.tar.gz
db_PLSDB_SRST2="" # Syntax: $Path_to_folder/plsdb.fna; Don´t forget to download the latest database: https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/
db_custom_SRST2="" # Syntax: $Path_to_folder/customDB.fasta
db_subspecies_Mab="" # Syntax: $Path_to/Minias2020-Steindor2019_SubspeciesMarkers.fasta #Available @ https://github.com/ngs-fzb/NTMtools/Gene_databases ?
db_MLST_SRST2="" # Syntax: $Path_to_folder; The most recent version of the pubMLST database is automatically downloaded, no need to download manually
db_kraken2="" # Syntax: $Path_to_folder; Don´t forget to download the latest databases: https://benlangmead.github.io/aws-indexes/k2
db_krona="" # Syntax: $Path_to_folder; Don´t forget to download the latest databases, if errors, see https://github.com/bioconda/bioconda-recipes/issues/10959
db_AMRfinder="" # Syntax: $Path_to_folder; Don´t forget to download the latest databases
fastANI_ref="" # Syntax: $Path_to_folder/Ref.fasta

##Species##
species="" #Options: e.g. Mycobacteroides abscessus, Mycobacterium avium
genome_size="" #Expected genome size in Mbp; For example "5.1M" for M. abscessus. Other species: check median genome size on NCBI
genome_size_full="" #Expected genome size in bp; For example "5100000" for M. abscessus. Other species: check median genome size on NCBI

###########################################
###!!!CHOOSE ANALYSES YOU WANT TO DO!!! ###
###########################################

#[options: Yes or No]
#Estimated times for 30 Mab samples (technical validation set Diricks et al. 2022) using Intel® Xeon® processor (E5-2650 v4 @ 2.2 GHz) with 48 Gb RAM and 16 CPUs
##Read quality control##
Do_multiqc="No"
###FastQ preprocessing###
Do_fastp="No"
##Subsampling
Do_subsampling="No" 
cov=100 #reads will be downsampled to reach this theoretic coverage; Default: 100"
#Estimate time: ~8 min 
##Assembly##
Do_assembly="No" #Using shovill
#Estimate time: ~3.5 hours 
##Taxonomy and phylogeny##
Do_subspecies="No" #Based on subspecies specific markers for Mab (Steindor et. al 2019 and Minias et. al 2020)
#Estimated time: ~10 min
Do_kraken2="No"
#Estimated time: ~250 min (4 hours)
Do_fastANI="No"
#Estimated time:"
Do_mash_fastQ="No"
#Estimated time: 2 min
Do_mash_fastA="No"
#Estimated time: 10 sec
Do_MLST_fastQ="No" # compare with mlst_tseemann?
#Estimated time: 18 min
##Ristance prediction##
Do_NTMprofiler="No" #Taxonomy and resistance prediction starting from raw reads
#Estimated time: 98 min
##Plasmid prediction##
Do_SRST2_PLSDB="No" #Against PLSDB database
#Estimated time: ~8 hours for creation of database files. This should only be done once (except if there is a new updated database available: check https://ccb-microbe.cs.uni-saarland.de/plsdb/ - typically updated once every year or every two years)
#Estimated time: ~2.5 hours for actual analysis
Do_SRST2_customDB="No" # Against custom database
#Estimated time: depends on size of plasmid database!
Do_extract_plasmids="No" #Extract fasta sequences for which a hit was found (in SRST2_PLSDB analysis); Default: No
#Estimated time: ~0.5 min per plasmid
Do_plasmidspades="No"
#Estimated time: ~200 min
Do_platon="No" # Requires assemblies!
#Estimated time: ~67 min
Do_amrfinder="No" #Requires assemblies!
#Estimated time: 

#############################################
###OPTIONAL PARAMATERS THAT YOU CAN CHANGE###
#############################################

set="Test" #Name of sample set - used for file naming"
#cpu: Resources to be used (max. amount depends on your hardware configuration)
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') # $(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') Uses all available threads

Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz
Rv="_R2" #Change if your fastQ files are not SampleName_R2.fastq.gz
PATH_tmp="$HOME/tmp" #Directory were temporary files are stored; Default: $HOME/tmp
#Note: It is important that PATH_tmp has enough space, otherwise some tools might fail!! 
PATH_fastA="$PATH_output/Assemblies/FinalAssemblies" # Path to folder where fastA files are stored; $PATH_output/Assemblies/FinalAssemblies = Default output if Do_assembly = Yes
ass="spades" #assembler that needs to be used (choose skesa velvet megahit or spades); Default: spades"

##########################################
###NO CHANGES REQUIRED BELOW THIS POINT###
##########################################

###############################################################################CODE#################################################################################
#mkdir -p $PATH_output

###QUALITY CONTROL READS###
if [[ "$Do_multiqc" == "Yes" ]]
then
	if [[ ! -z "$PATH_fastQ" ]] || [[ ! -z "$PATH_output" ]]
	then
		echo "Starting to perform quality control on reads"
		bash $PATH_scripts/starter_multiqc.sh -i $PATH_fastQ -o $PATH_output -c $cpu -e $conda_env_multiqc
	else
		echo "Please provide all required arguments (-i PATH_fastQ and -o PATH_output )! use starter_multiqc.sh -h for help on syntax"
	fi
elif [[ "$Do_multiqc" == "No" ]]
then
	echo "quality control on reads skipped."
else 
	echo "Please decide on whether you want to do quality control on reads or not!"
fi

###Fastq preprocessing###
if [[ "$Do_fastp" == "Yes" ]]
then
	if [[ ! -z "$PATH_fastQ" ]] || [[ ! -z "$PATH_output" ]]
	then
		echo "Starting to preprocess reads (e.g. trimming)"
		bash $PATH_scripts/starter_fastp.sh -i $PATH_fastQ -o $PATH_output -c $cpu -e $conda_env_fastp
	else
		echo "Please provide all required arguments (-i PATH_fastQ and -o PATH_output )! use starter_fastp.sh -h for help on syntax"
	fi
elif [[ "$Do_fastp" == "No" ]]
then
	echo "Preprocessing of reads skipped."
else 
	echo "Please decide on whether you want to do preprocessing of reads or not!"
fi

###SUBSAMPLING###
if [[ "$Do_subsampling" == "Yes" ]]
then
	if [[ ! -z "$PATH_fastQ" ]] || [[ ! -z "$PATH_output" ]] || [[ ! -z "$genome_size_full" ]] || [[ ! -z "$cov" ]]
	then
		echo "Starting to subsample"
		bash $PATH_scripts/starter_subsampling.sh -i $PATH_fastQ -o $PATH_output -g $genome_size_full -d $cov -c $cpu -e $conda_env_seqkit
	else
		echo "Please provide all required arguments (-i PATH_fastQ, -o PATH_output, -g genome_size, -d cov )! use starter_subsampling.sh -h for help on syntax"
	fi
elif [[ "$Do_subsampling" == "No" ]]
then
	echo "Subsampling skipped."
else 
	echo "Please decide on whether you want to do subsampling or not!"
fi


###ASSEMBLY###
if [[ "$Do_assembly" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$genome_size" ]]
	then
		#genome_size_shovill=$()
		echo "Starting to assemble"
		bash $PATH_scripts/starter_Assemble_usingShovill.sh -i $PATH_fastQ -o $PATH_output -g $genome_size -d $cov -c $cpu -e $conda_env_shovill -f $Fw -r $Rv -s $set -a $ass
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output and -d db_subspecies_Mab!"
	fi
elif [[ "$Do_assembly" == "No" ]]
then
	echo "Assembly skipped."
else 
	echo "Please decide on whether you want to make an assembly or not!"
fi

###TAXONOMY AND PHYLOGENY###
if [[ "$Do_subspecies" == "Yes" ]] && [[ "$species" == "Mycobacteroides abscessus" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_subspecies_Mab" ]]
	then 
		echo "Starting subspecies classification"
		bash $PATH_scripts/starter_Subspecies_Mab_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $db_subspecies_Mab -c $cpu -e $conda_env_SRST2 -f $Fw -r $Rv -s $set -t $PATH_tmp
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output and -d db_subspecies_Mab!"
	fi
elif [[ "$Do_subspecies" == "Yes" ]] && [[ "$species" != "Mycobacteroides abscessus" ]]
then
	echo "We don´t have subspecies classification available for "$species". Please contact mdiricks@fz-borstel.de if this should be implemented."
elif [[ "$Do_subspecies" == "No" ]]
then
	echo "Subspecies analysis skipped."
else 
	echo "Please decide on whether you want to do subspecies classification or not (Choose Yes or No)!"
fi

if [[ "$Do_mash_fastQ" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]]
	then 
		echo "Starting mash on FastQ files"
		bash $PATH_scripts/starter_Mashtree.sh -i $PATH_fastQ -o $PATH_output -n .fastq.gz -e $conda_env_mashtree -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastQ and -o PATH_output!"
	fi
elif [[ "$Do_mash_fastQ" == "No" ]]
then
	echo "mash on FastQ files skipped."
else 
	echo "Please decide on whether you want to do mash on FastQ files or not (Choose Yes or No)!"
fi

if [[ "$Do_mash_fastA" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastA" ]]
	then 
		echo "Starting mash on FastA files"
		bash $PATH_scripts/starter_Mashtree.sh -i $PATH_fastA -o $PATH_output -n .fasta -e $conda_env_mashtree -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastA and -o PATH_output!"
	fi
elif [[ "$Do_mash_fastA" == "No" ]]
then
	echo "mash on FastA files skipped."
else 
	echo "Please decide on whether you want to do mash on FastA files or not (Choose Yes or No)!"
fi

if [[ "$Do_kraken2" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_kraken2" ]] && [[ ! -z "$db_krona" ]]
	then 
		echo "Starting kraken2 on FastQ files"
		bash $PATH_scripts/starter_kraken2.sh -i $PATH_fastQ -o $PATH_output -d $db_kraken2 -k $db_krona -e $conda_env_kraken2
		perl $PATH_scripts/kraken_parse_results.v2.0.0.pl -s "$species" $PATH_output/kraken2/*.report
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output, db_krona and db_kraken2!"
	fi
elif [[ "$Do_kraken2" == "No" ]]
then
	echo "kraken2 on FastQ files skipped."
else 
	echo "Please decide on whether you want to do kraken2 on FastQ files or not (Choose Yes or No)!"
fi

if [[ "$Do_fastANI" == "Yes" ]]
then
	if [[ ! -z "$PATH_fastA" ]] && [[ ! -z "$PATH_output" ]] && [[ ! -z "$fastANI_ref" ]]
	then 
		echo "Starting fastANi"
		bash $PATH_scripts/starter_fastANI.sh -i $PATH_fastA -o $PATH_output -r $fastANI_ref -e $conda_env_fastANI
	else
		echo "Please provide all required arguments: -i PATH_fastA, -o PATH_output and -r reference!"
	fi
elif [[ "$Do_fastANI" == "No" ]]
then
	echo "FastANI skipped."
else 
	echo "Please decide on whether you want to do fastANI or not (Choose Yes or No)!"
fi

if [[ "$Do_MLST_fastQ" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_MLST_SRST2" ]] && [[ ! -z "$species" ]]
	then 
		echo "Starting MLST on FastQ files"
		bash $PATH_scripts/starter_MLST_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $db_MLST_SRST2 -e $conda_env_SRST2 -p "$species"
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output, db_MLST_SRST2 and species!"
	fi
elif [[ "$Do_MLST_fastQ" == "No" ]]
then
	echo "MLST on FastQ files skipped."
else 
	echo "Please decide on whether you want to do MLST on FastQ files or not (Choose Yes or No)!"
fi

###TAXONOMY AND RESISTANCE PREDICTION###

if [[ "$Do_NTMprofiler" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]]
	then 
		echo "Starting running NTMprofiler"
		bash $PATH_scripts/starter_NTMprofiler.sh -i $PATH_fastQ -o $PATH_output -n .fastq.gz -c $cpu -e $conda_env_NTMprofiler -f $Fw -r $Rv -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastQ and -o PATH_output!"
	fi
elif [[ "$Do_NTMprofiler" == "No" ]]
then
	echo "running NTMprofiler skipped."
else 
	echo "Please decide on whether you want to run NTMprofiler or not (Choose Yes or No)!"
fi

####PLASMID DETECTION###
###PART1
if [[ "$Do_plasmidspades" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]]
	then
		echo "Starting plasmid assembly using plasmidspades"
		bash $PATH_scripts/starter_plasmidspades.sh -i $PATH_fastQ -o $PATH_output -c $cpu -e $conda_env_spades -f $Fw -r $Rv -s $set -t $PATH_tmp
		bash $PATH_scripts/starter_Platon.sh -i $PATH_output/Plasmidspades/FinalAssemblies -o $PATH_output/Plasmidspades/FinalAssemblies -d $db_platon -c $cpu -e $conda_env_platon -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastQ -o PATH_output and -d db_PLSDB_mash!"
	fi
elif [[ "$Do_plasmidspades" == "No" ]]
then
	echo "Plasmid assembly using plasmidspades skipped."
else 
	echo "Please decide on whether you want to do assembly using plasmidspades or not (Choose Yes or No)!"
fi

###PART2a
if [[ "$Do_SRST2_PLSDB" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_PLSDB_SRST2" ]]
	then
		if [ ! -f "$(dirname $db_PLSDB_SRST2)/plsdb_SRST2.fasta" ]
		then
			echo "Starting converting plsdb database into SRST2 compatible database"
			bash $PATH_scripts/PLSDB_to_SRST2_db.sh -d $db_PLSDB_SRST2
		else
			echo "SRST2 compatible plsdb database found"
		fi
		echo "Starting plasmid detection with SRST2 using PLSDB database"
		bash $PATH_scripts/starter_PLSDB_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $(dirname $db_PLSDB_SRST2)/plsdb_SRST2.fasta -c $cpu -e $conda_env_SRST2 -f $Fw -r $Rv -s $set -t $PATH_tmp
	else
		echo "Please provide all required arguments: -i PATH_fastQ -o PATH_output and -d db_PLSDB_SRST2!"
	fi
elif [[ "$Do_SRST2_PLSDB" == "No" ]]
then
	echo "Plasmid detection with SRST2 using PLSDB database skipped."
else 
	echo "Please decide on whether you want to do plasmid detection with SRST2 using PLSDB database or not (Choose Yes or No)!"
fi

if [[ "$Do_extract_plasmids" == "Yes" ]]
then
	bash $PATH_scripts/Extract_fasta_from_multifasta.sh -i $PATH_output/PLSDB_SRST2/Plasmid_hits.txt -o $PATH_output/PLSDB_SRST2/Plasmid_hits -e $conda_env_seqkit -m $(dirname $db_PLSDB_SRST2)/plsdb_SRST2.fasta
	bash $PATH_scripts/starter_Mashtree.sh -i $PATH_output/PLSDB_SRST2/Plasmid_hits -o $PATH_output/PLSDB_SRST2/Plasmid_hits -e $conda_env_mashtree -n .fasta
fi

###PART2b
if [[ "$Do_SRST2_customDB" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_custom_SRST2" ]]
	then
		customDB=$(basename $db_custom_SRST2| cut -d '.' -f 1)
		customDB_SRST2=$customDB"_SRST2.fasta"
		if [ ! -f "$(dirname $db_custom_SRST2)/$customDB_SRST2" ]
		then
			echo "Starting converting custom database into SRST2 compatible database"
			bash $PATH_scripts/customDB_to_SRST2_db.sh -d $db_custom_SRST2
		else
			echo "SRST2 compatible custom database found"
		fi
		echo "Starting plasmid detection with SRST2 using custom database"
		bash $PATH_scripts/starter_customDB_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $(dirname $db_custom_SRST2)/$customDB_SRST2 -c $cpu -e $conda_env_SRST2 -f $Fw -r $Rv -s $set -t $PATH_tmp
	else
		echo "Please provide all required arguments: -i PATH_fastQ -o PATH_output and -d db_custom_SRST2!"
	fi
elif [[ "$Do_SRST2_customDB" == "No" ]]
then
	echo "Plasmid detection with SRST2 using custom database skipped."
else 
	echo "Please decide on whether you want to do plasmid detection with SRST2 using custom database or not (Choose Yes or No)!"
fi

###PART3
if [[ "$Do_platon" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastA" ]] && [[ ! -z "$db_platon" ]]
	then
		echo "Plasmid detection with platon"
		bash $PATH_scripts/starter_platon.sh -i $PATH_fastA -o $PATH_output -d $db_platon -c $cpu -e $conda_env_platon -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastA, -o PATH_output and -d db_platon!"
	fi
elif [[ "$Do_platon" == "No" ]]
then
	echo "Plasmid detection with platon skipped."
else 
	echo "Please decide on whether you want to do plasmid detection with platon or not (Choose Yes or No)!"
fi

echo "Script Finished!"

####RESISTANCE AND VIRULENCE DETECTION###
if [[ "$Do_amrfinder" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastA" ]] && [[ ! -z "$db_AMRfinder" ]]
	then
		echo "Resistance gene and virulence detection with AMRFinderplus"
		bash $PATH_scripts/starter_amrfinder.sh -i $PATH_fastA -o $PATH_output -c $cpu -e $conda_env_AMRfinder -d $db_AMRfinder
	else
		echo "Please provide all required arguments: -i PATH_fastA, -d db_AMRfinder and -o PATH_output!"
	fi
elif [[ "$Do_amrfinder" == "No" ]]
then
	echo "Resistance gene and virulence detection with AMRFinderplus skipped."
else 
	echo "Please decide on whether you want to do resistance gene and virulence detection with AMRFinderplus or not (Choose Yes or No)!"
fi

echo "Script Finished!"
exit


####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################
###Download and extract PLSDB database### 
#https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/
#db_PLSDB_SRST2: Download FASTA archive. Use bzip2 -d db_name to extract
#Convert fna file to SRST2 compatible database: automated in NTMseq pipeline

###Download and extract latest platon database### 
#Info: Platon depends on a custom database based on MPS, RDS, RefSeq Plasmid database, PlasmidFinder db as well as manually curated MOB HMM models from MOBscan, custom conjugation and replication HMM models and oriT sequences from MOB-suite.
#Download page: https://github.com/oschwengers/platon#database or https://zenodo.org/search?page=1&size=20&q=conceptrecid:3349651&all_versions&sort=-version

