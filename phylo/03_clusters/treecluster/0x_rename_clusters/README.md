# 0x_rename_clusters

Renumber per-segment TreeCluster clusters so the number reflects **phylogenetic position**
instead of TreeCluster's traversal order.

## Why
TreeCluster's original cluster numbers correlate only weakly/inconsistently with display
position (Spearman 0.32–0.75 across segments), so the integer carries no reliable
phylogenetic meaning. `max_clade` clusters are perfectly contiguous clades (contiguity = 1),
so ordering by tip position is unambiguous.

## What `rename_clusters_v1.Rmd` does
For each segment, using the rooted ML IQ-TREE (`../../../02_iqtree/iq_gtrfr10_rooted/`, via the
ggtree ladderized layout — same as the figures):
- orders the numbered clusters by **lowermost tip (min y) ascending = bottom → top**
- renumbers `{SEG}-1 … {SEG}-N` (bottom of the tree = 1)
- leaves the singleton bucket `{SEG}-s` unchanged

The `tc_geno` / `TC#` alias rebuild is **validation only** here (proves the relabel is a clean
bijection): distinct `tc_geno` unchanged at 5042; old→new 1:1; 0 unmapped IDs. The
authoritative regeneration of `tc_geno` / alias happens in the downstream 08_all-together build.

## Outputs (`<date>_rename_clusters_out/`)
**Authoritative:**
| File | Contents |
|---|---|
| `tc_cluster_rename_lookup_{SEG}.csv` | old → new per segment (+ min_y, mean_y, size) |
| `tc_cluster_rename_lookup_all.csv` | all segments combined |
| `newnumber_vs_position.tiff` | sanity plot: new number vs tip position (monotone) |

**`validation/` (reference only, NOT authoritative):**
| File | Contents |
|---|---|
| `VALIDATION_tc_cluster_renamed_samples.csv` | per-sample old + `_new` columns + `tc_geno_new` |
| `VALIDATION_tc_geno_rename_lookup.csv` | old `tc_geno` → new `tc_geno` |
| `VALIDATION_tc_alias_rename.csv` | `TC#` alias → old/new `tc_geno` |

## Downstream (must be re-run with the new IDs)
The **08_all-together** metadata build (which feeds the `data/08_clean-share/` copy) should
apply the per-segment lookups to `tc_cluster{SEG}` and rebuild `tc_geno` /
`treecluster_genotype_alias` from the renamed columns. The `08_clean-share` copy
(`rename_samp_seq` step) then re-runs, and
`phylo/03_clusters/compare_phylo_clusters/compare_phylo_clusters_all.Rmd` +
`phylo/03_clusters/subtree_cluster_matrix/` pick up the new numbering from the metadata.

## Note
Ordering is tied to ggtree's default ladderization (deterministic for a fixed tree). If the
IQ-TREE trees are regenerated, re-run this step.
