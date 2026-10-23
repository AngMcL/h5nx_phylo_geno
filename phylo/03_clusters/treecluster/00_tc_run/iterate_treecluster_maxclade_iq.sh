#!/bin/bash

# Get current date (format: YYYY-MM-DD)
date_str=$(date +%F)

# Define output directory with date
outdir="${date_str}_ft_treecluster_maxclade_out"
mkdir -p "$outdir"

# Define the 8 influenza segments
segments=("HA" "NA" "M" "NS" "NP" "PA" "PB1" "PB2")

# Loop over each segment and threshold (0.01 to 0.25 in steps of 0.01)
for seg in "${segments[@]}"; do
    for t in $(seq 0.01 0.01 0.25); do
        echo "Running TreeCluster for segment ${seg} at threshold ${t}"

        # Input tree file
        infile="../../../02_iqtree/iq_gtrfr10_rooted/trim_nodup_nuc_${seg}_rooted_remiqtree_rooted.treefile"
        # Output file in dated subfolder
        outfile="${outdir}/${seg}_cluster_maxclade_t${t}.csv"

        # Run TreeCluster
        TreeCluster.Py -i "$infile" -o "$outfile" -t "$t" -m max_clade
    done
done
