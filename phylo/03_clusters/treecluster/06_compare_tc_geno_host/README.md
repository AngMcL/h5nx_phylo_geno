# 06_compare_tc_geno_host

Compare host composition (taxa2, order) within treecluster genotypes.

## Goal

For each treecluster genotype (`tc_geno`), evaluate how many unique host types are represented and how sequences are distributed across hosts. Genotypes with low sample counts or only one host type are excluded.

## Inputs

| File | Description |
|------|-------------|
| `data/08_clean-share/meta/metadata.csv` | Sample metadata with `tc_geno`, `taxa2`, `order` columns |
| `data/06_host/multi-taxa/taxa_palette/taxa_colors_refined_contrast.csv` | Hex colour palette for taxa2 |

## Parameters (top of Rmd)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_n` | 100 | Minimum sequences per genotype to include |

Genotypes must also have > 1 unique taxa2 to be included.

## Outputs

Saved to `<YYYY-MM-DD>_out/`:

| File | Description |
|------|-------------|
| `*_tc_geno_host_taxa2_composition.pdf/png` | Stacked bar: n sequences per genotype, filled by taxa2, ranked by n unique taxa2 |
| `*_tc_geno_host_order_composition.pdf/png` | Stacked bar: n sequences per genotype, filled by order, ranked by n unique orders |
| `*_tc_geno_host_combined.pdf/png` | Patchwork panel combining both plots |
| `*_tc_geno_host_summary.csv` | Per-genotype counts of n_seq, n_taxa2, n_order |

## Analysis

`compare_tc_geno_host.Rmd`
