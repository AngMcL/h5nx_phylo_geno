# compare_phylo_clusters

Compare cluster and genotype labels on shared ML trees using `ggtree` + `gheatmap`.

## Goal

Visualize treecluster (threshold-based) cluster assignments and GenoFlu genotype labels side-by-side on the same phylogeny for each of the 8 influenza segments. Each tree is annotated with two `gheatmap` columns — one per labelling scheme — to allow direct visual comparison of concordance and discordance.

## Inputs

| Source | Description |
|---|---|
| `02_iqtree/iq_gtrfr10_rooted/` | Rooted ML trees (one `.treefile` per segment) |
| `03_clusters/01_tc_evaluate/` | Treecluster cluster assignments per segment |
| `data/05_annotate/02_genoflu/` | GenoFlu per-segment genotype assignments |

## Outputs

- `YYYY-MM-DD_compare_clusters_out/` — per-segment ggtree + gheatmap figures (`.tiff`)
