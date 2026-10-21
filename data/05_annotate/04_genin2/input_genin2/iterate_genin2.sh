#!/bin/bash

#USAGE
#copy the script into the folder
#./iterate_genin2.sh

#read all files matching pattern .fasta in the current folder
alignfiles=(*.fasta)
echo "identified ${#alignfiles[@]} alignments"

#run genin on all the alignfiles
for i in "${alignfiles[@]}"; do
    echo "initiating genin for $i"
    genin2 -o "${i}_output.tsv" "$i" 
done