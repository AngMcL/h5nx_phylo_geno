#!/bin/bash

#USAGE
#./iterate_genoflu.sh

#read all files matching pattern .fasta in the current folder
files=(*.fasta)
echo "identified ${#files[@]} fasta"

#go through each fasta, call genoflu, export results to file by samplename
for i in "${files[@]}"; do
    echo "calling genoflu for $i"
    genoflu.py -f "$i"
done
