# composite_subtree_geno

Manuscript composite figure combining pre-rendered subtree-cluster-matrix and genotype-tree panels

## Goal

Assemble a 2x3 composite for the main manuscript: concordance matrices comparing treecluster genotypes (`tc_geno`) against GenoFlu and Genin2 genotypes, each paired with example genotype trees for segments that disagree within a top genotype for that scheme

## Inputs (pre-rendered TIFFs, not regenerated here)

| Panel | Source | File |
|---|---|---|
| A | `../2026-08-28_subtree_cluster_matrix_out/` | `tc_geno_vs_genoflu_geno_matrix_500_50_highlight.tiff` |
| B | `../../compare_phylo_clusters/2026-07-01_compare_clusters_out/genotype_trees/genoflu_D1.1/` | `PB2_genoflu_D1.1_gheatmap.tiff` |
| C | `../../compare_phylo_clusters/2026-07-01_compare_clusters_out/genotype_trees/genoflu_B3.13/` | `PA_genoflu_B3.13_gheatmap.tiff` |
| D | `../2026-08-28_subtree_cluster_matrix_out/` | `tc_geno_vs_genin_geno_matrix_100_50_highlight.tiff` |
| E | `../../compare_phylo_clusters/2026-07-01_compare_clusters_out/genotype_trees/genin_EA-2020-C/` | `NP_genin_EA-2020-C_gheatmap.tiff` |
| F | `../../compare_phylo_clusters/2026-07-01_compare_clusters_out/genotype_trees/genin_EA-2020-C/` | `PA_genin_EA-2020-C_gheatmap.tiff` |

## Layout

- Row 1 (A, B, C): "Genoflu vs. TreeCluster" matrix + two GenoFlu genotype trees (D1.1, B3.13)
- Row 2 (D, E, F): "Genin2 vs. TreeCluster" matrix + two Genin2 EA-2020-C genotype trees (NP, PA segments)
- A and D are the `_highlight` matrix variants — dark-red outline around the B3.13/D1.1 (A) and EA-2020-C (D) columns, matching the genotypes shown in B/C and E/F
- Column widths shared across both rows (matrix column sized off the average native aspect ratio of A/D, tree columns off the trees' shared aspect ratio) so rows align
- Row heights `1.4:1` (matrix row taller than tree row)
- Each panel gets its own caption strip above the image rather than an overlaid label, to avoid covering the image's own axis/column-name text

## Outputs

- `YYYY-MM-DD_composite_subtree_geno_out/composite_subtree_geno.tiff`
