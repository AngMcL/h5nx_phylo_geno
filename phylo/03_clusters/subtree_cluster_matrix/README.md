# subtree_cluster_matrix

Compute and visualize concordance matrices comparing whole-genome treecluster genotypes against GenoFlu and Genin2 genotype labels.

## Goal

Build a pairwise sequence-count matrix between treecluster genotypes (`tc_geno`) and each of the GenoFlu (`genoflu_geno`) and Genin2 (`genin_geno`) genotype labels. Each cell is a heatmap tile coloured and annotated by sequence count, quantifying where the labelling schemes agree or diverge.

## Inputs

| Source | Description |
|---|---|
| `data/08_clean-share/meta/metadata.csv` | Compiled metadata with `tc_geno`, `genoflu_geno`, `genin_geno`, and per-segment cluster columns |

## Threshold sweep

Instead of single fixed minimum-count cutoffs, every matrix is regenerated for all unique combinations of two threshold vectors (set in the `ins/outs` chunk):

- `thresh.geno` — minimum total sequences per GenoFlu / Genin2 genotype: `100, 500, 1000, 10000`
- `thresh.tc` — minimum total sequences per `tc_geno`: `10, 50, 100, 500`

This produces `length(thresh.geno) * length(thresh.tc)` = 16 matrices per genotype scheme (32 files total). Combinations where no group passes a threshold are skipped and logged.

## Outputs

- `YYYY-MM-DD_subtree_cluster_matrix_out/` — concordance matrix heatmaps (`.tiff`), one per threshold combination. File names carry both thresholds as `tc_geno_vs_genoflu_geno_matrix_{thresh.geno}_{thresh.tc}.tiff` and `tc_geno_vs_genin_geno_matrix_{thresh.geno}_{thresh.tc}.tiff`.
