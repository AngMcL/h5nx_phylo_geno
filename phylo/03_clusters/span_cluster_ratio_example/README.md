# span_cluster_ratio_example

A toy, simulated phylogeny illustrating the span / cluster ratio reported for genotypes in `subtree_cluster_matrix` (`geno_tc_cluster_summary.csv`, `geno_tc_cluster_spread.csv`, Tables S5/S6). No real data are used.

## Definition

For one genotype and one segment, using numeric cluster IDs only (singleton `{SEG}-s` buckets excluded):

- `n_clusters` = number of distinct clusters the genotype's sequences fall in
- `span` = max cluster ID - min cluster ID + 1
- `ratio` = span / n_clusters (1 decimal)

Cluster IDs are numbered bottom to top of the ladderized segment tree, so a ratio of 1 means the genotype occupies a contiguous block of clusters, and a ratio above 1 means it is scattered across the tree with other clusters in between.

## Scenarios (HA only, one tree with 7 clusters)

| Scenario | Genotype | Genotype clusters | n_clusters | span | ratio |
|---|---|---|---|---|---|
| `01_single_cluster` | X | HA_2 | 1 | 1 | 1.0 |
| `02_contiguous_ratio1` | Y | HA_4, HA_5, HA_6 | 3 | 3 | 1.0 |
| `03_noncontiguous` | Z | HA_1, HA_3, HA_7 | 3 | 7 | 2.3 |

Every cluster on the tree belongs to exactly one genotype.

## Simulation

- Caterpillar backbone (`ape::stree(type = "left")`) with one tip per cluster and long random branch lengths, so the ladder is visible
- Each backbone tip is replaced (`ape::bind.tree`) by a random subtree of 2-4 tips with short branches, so every cluster is a monophyletic group of closest sequences
- Tree is ladderized with `ape::ladderize(right = FALSE)` and drawn with `ggtree(ladderize = FALSE)` so the plotted order is exactly the ladderized order; clusters are numbered bottom to top by mean tip position, matching `02_rename_clusters`
- Seed, number of clusters, tips per cluster, branch-length ranges and scenarios are set in the `ins/outs` chunk

## Plot

One tree with all three genotypes drawn together (`span_ratio_tree.png`):

- Drawn with `ggtree`; tips coloured by cluster with a discrete viridis scale
- Each cluster in a genotype gets a triangle (apex at its MRCA, base just past its tips) filled by genotype: X red, Y purple, Z orange
- Cluster labels (`HA_1`, ...) sit immediately right of each cluster's furthest tip, vertically centred
- Summary table above the tree (Genotype, clusters, span, n_clusters, ratio); the Genotype column is the genotype legend (fill key + name, centred)
- "HA cluster" legend directly below the table, highest cluster on top to match the tree; table and legend both left-aligned at the plot edge
- No title, no x axis (branch lengths are arbitrary), white background, enlarged text

## Outputs

`YYYY-MM-DD_summary/`:

- `span_ratio_tree.png` — the combined tree
- `span_ratio_summary.csv` — one row per genotype with member clusters, min/max ID, n_clusters, span, ratio
- `span_ratio_tip_clusters.csv` — tip-to-cluster assignments with the genotype each cluster belongs to (NA = none)
- `sim_tree.nwk` — the simulated tree
