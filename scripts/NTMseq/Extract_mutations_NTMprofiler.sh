#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#To extract mutations from a specific gene, mentioned in the sample-specific result file of NTMprofiler into excel file

###Required packages###
#None

###Required parameters - via command line###
# None


PATH_input="$HOME/auto/Thecus_mdiricks/Mabscessus/2020_Cronoclone/NTMseq/WithFastP/NTMprofiler"
output="mutation_results.txt"
GOI="MAB_r5051" #Gene of interest

# Initialize the Excel file with header
echo -e "Isolate ID\tFraction_GOI_missing\tMAB_r5051 Mutations\tMutation_fraction" > $PATH_input/$output

# Loop through each result file in the directory
for resultfile in $PATH_input/*.results.txt; do
    isolate_id=$(basename "$resultfile" .results.txt)
    mutations=$(grep $GOI "$resultfile" | awk '{ print $4 }' | head -n -1 | paste -s -d";")
	fraction=$(grep $GOI "$resultfile" | awk '{ print $5 }' | paste -s -d";")
	coverage=$(grep "MAB_r5051" -A1 "$resultfile" | tail -n 1 | awk '{ print $4 }')
	#echo $mutations

    echo -e "$isolate_id\t$coverage\t$mutations\t$fraction" >> $PATH_input/$output
done

echo "Data has been saved to $output"