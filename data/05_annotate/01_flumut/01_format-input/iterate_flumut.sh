#!/bin/bash

#USAGE
#copy the script into the folder
#./iterate_flumut.sh

#read all files matching pattern .fasta in the current folder
alignfiles=(*.fasta)
echo "identified ${#alignfiles[@]} alignments"

#go through each alignment, identify tree file corresponding and run all through iqtree model selection
for i in "${alignfiles[@]}"; do
    echo "initiating flumut for $i"
    flumut -m "${i}_markers_output.tsv" -M "${i}_mutations_output.tsv" -l "${i}_lit_output.tsv" "$i" 
done

