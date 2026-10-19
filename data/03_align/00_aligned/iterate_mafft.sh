#!/bin/bash

#USAGE
# cd to folder with unaligned fasta
#./ Iterate_mafft.sh

#read all files matching pattern .fasta in the current folder
files=(*.fasta)
echo "identified ${#files[@]} un-aligned seq files"

#go through each alignment to align with mafft
for i in "${files[@]}"; do
    echo "aligning $i"
    outfile=$(echo "$i" | sed 's,.fasta,_aligned.fasta,g')
    mafft "$i" > "$outfile"
done

#make a folder and mv output files into it
mkdir "aligned_fasta"
mv *aligned.fasta aligned_fasta/.