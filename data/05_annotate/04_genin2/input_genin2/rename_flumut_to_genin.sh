#!/bin/bash

# Rename files in the current directory:
# replace all instances of "flumut" with "genin"

for file in *flumut*; do
    # skip if no matching files
    [ -e "$file" ] || continue

    newfile="${file//flumut/genin}"

    # avoid overwriting an existing file
    if [ -e "$newfile" ]; then
        echo "Skipping: $file -> $newfile already exists"
    else
        echo "Renaming: $file -> $newfile"
        mv -- "$file" "$newfile"
    fi
done