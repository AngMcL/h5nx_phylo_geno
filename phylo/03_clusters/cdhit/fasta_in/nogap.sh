#!/bin/bash

for f in *.fasta; do
  [ -e "$f" ] || continue
  out="${f%.*}.no_gaps.fa"
  seqkit seq -g "$f" > "$out"
  echo "Processed $f -> $out"
done