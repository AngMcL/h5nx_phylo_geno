#!/bin/bash

# Define the genome segments (edit as needed)
segments=("PB2" "PB1" "PA" "HA" "NP" "NA" "M" "NS")

# Define the clustering thresholds to evaluate
thresholds=(85 86 87 88 89 90 91 92 93 94 95 95.5 96 96.5 97 97.5 98 98.5 99 99.5)

# Output directory (make sure it exists)
outdir="cdhit_out_g1"
mkdir -p "$outdir"

# Loop over each segment and threshold
for seg in "${segments[@]}"; do
    infile="fasta_in/trim_nodup_nuc_${seg}_rooted_rem.no_gaps.fa"
    if [[ ! -f "$infile" ]]; then
        echo "Warning: input file not found for $seg — skipping."
        continue
    fi

    for thr in "${thresholds[@]}"; do
        outfile="${outdir}/cdhit_${seg}_${thr}_g1.fasta"

        echo "Running CD-HIT for $seg at ${thr}% identity..."
        cd-hit-est -i "$infile" -o "$outfile" -c $(echo "$thr / 100" | bc -l) -n 10 -d 0 -g 1 -M 0 -T 8
    done
done

echo "All CD-HIT runs complete."
