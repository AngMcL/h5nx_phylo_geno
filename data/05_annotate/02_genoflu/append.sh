#!/bin/bash

#to append all tsvs with matching colnames into one file
#Usage: ./append.sh
output="combined_geno_samples.tsv"
tsv_files=(*.tsv)

# Remove output file if it already exists - overwrites
rm -f "$output"

# All tsv files 
tsv_files=(*.tsv)

# Check
if [ ${#tsv_files[@]} -eq 0 ]; then
  echo "No TSV files found."
  exit 1
fi

# Write header from the first file
head -n 1 "${tsv_files[0]}" > "$output"

# Loop through all files and append content (excluding headers)
for file in "${tsv_files[@]}"; do
  tail -n +2 "$file" >> "$output"
done

echo "Combined ${#tsv_files[@]} TSV files into '$output'"