# 06_host

## annotate_host_v2/3.Rmd
merge host_old into common by correcting spelling, redundancies, latin v common
requires ncbi api access for taxa lookup

## mallard_duck_fix.Rmd
identified issue on 2026-06-18 where too many mallards in final dataset, and found mistaken duck species grouped as mallards
fixed and this had to be fixed in multi_taxa

## treemaps/
move the treemap chunk to own subfolder
treemaps_all: full cleaned dataset H5Nx (with mallard fix)
treemaps_2344b: filtered to 2.3.3.4b since 2021 for malcolm

## multi_taxa/
repeat this with mallard fix 
also add treemaps for taxa2 final groupings
outputs lookup of sample to mallard-corrected taxa2
TO DO add larger grouping threshold for fewer classes: taxa3