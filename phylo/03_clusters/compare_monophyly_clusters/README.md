# compare_monophyly_clusters

Compare cluster quality (monophyly, purity, size) across four clustering methods
— **genoflu**, **genin**, **cdhit**, **treecluster** — over the 8 influenza segments
(PB2, PB1, PA, HA, NP, NA, M, NS).

## Script

`compare_monophyly_clusters.R` (run from this folder). Outputs go to a dated
`YYYY-MM-DD_cluster_metrics_compare/` folder (runs before 2026-09-22 are in
`cluster_metrics_compare/`).

## Inputs

Per-cluster metrics files (one row per cluster) are the **single source of truth**;
the segment-level aggregate files are no longer used:

| method | path |
|---|---|
| genoflu | `../genoflu_clades/2026-08-20_out/genoflu_cluster-metrics_all-segments.csv` |
| genin | `../genin_clades/2026-08-20_out/genin_cluster-metrics_all-segments.csv` |
| cdhit | `../cdhit/cluster_eval/cdhit_cluster-metrics_all-segments.csv` |
| treecluster | `../treecluster/02_tc_trees/cluster_eval/treecluster_cluster-metrics_all-segments.csv` |

Note: the cdhit and treecluster paths are bare `cluster_eval/` folders (not nested in a
dated `..._trees_out/` folder) because both runs picked up a stale `f.out` variable left
over from an earlier session (a different upstream script sets `f.out`, and
`if(exists("f.out"))out_dir<-f.out` in these two scripts can't tell it came from
somewhere else) — falls back to the literal `"cluster_eval"` default. Doesn't affect
correctness, just location; re-point these if a future re-run lands somewhere else.

Key columns: `segment, method, cluster, n_tips, type` (monophyletic / paraphyletic /
singleton), `purity`, `faith_pd`.

`faith_pd` is Faith's phylogenetic diversity per cluster: the total branch length
(subs/site) of the minimal subtree connecting a cluster's tips to their MRCA. For a
monophyletic cluster this equals the clade's total branch length; for a
paraphyletic/polyphyletic cluster it sums only the tip-to-MRCA paths of member tips,
so it captures how much evolutionary distance a cluster spans even when that
distance isn't reducible to a single clade. Computed upstream in the same
tree-traversal loop that assigns `purity`/`type`
(`genoflu_clades/plot_fill_clades_v3.Rmd`, `genin_clades/plot_fill_clades_genin_v3.Rmd`,
`cdhit/cdhit_plot_trees_g1.Rmd`, `treecluster/02_tc_trees/treecluster_plot_trees_v2.Rmd`).

## Conventions

- **Singletons excluded everywhere.** A single `filter(n_tips > 1)` is applied once,
  then *all* metrics are recomputed from the non-singleton clusters. Singletons carry
  `type == "singleton"` (not "monophyletic"), so including them would deflate
  `prop_clusters_monophyletic` and inflate purity; cdhit/treecluster have ~30–40%
  singletons, so this matters.
- **Median, not mean, for purity and Faith PD**, in every plot and summary statistic.
  Checked via the ridgeline plots (`faith-pd_ridgeline_by-method.png`,
  `purity_ridgeline_by-method.png`): both distributions are visibly right-skewed/
  non-normal, so mean isn't representative — a dashed mean line was tried on both
  ridgelines during development and dropped once this was confirmed. `N clusters` and
  cluster size use their natural count/median/max, not an average.
  `prop_clusters_monophyletic` (a proportion, not a central-tendency stat) is
  unaffected by this.
- **One size-bin scheme, used everywhere** (`size_breaks_seg`/`size_labs_seg`: 2-5,
  6-10, 11-50, 51-100, 100-200, 200-500, 500-1000, 1000-10000, 10000+). Earlier had
  separate coarse (pooled) vs fine (per-segment) bins; unified to fine bins throughout
  after finding real GenoFlu/Genin2 clusters well above 200 tips, which a single
  "200+" catch-all bin was hiding.
- **Colorblind-safe method palette**: Paul Tol's "vibrant" qualitative set, one fixed
  color per method (GenoFlu `#0077BB`, Genin2 `#EE7733`, CD-HIT `#009988`, TreeCluster
  `#CC3311`). Segment colors elsewhere in this pipeline (e.g. `geno_reassort.Rmd`)
  already use Tol's sister "bright" palette, so these are deliberately drawn from the
  other Tol set to avoid the two encodings looking alike.
- **Faith PD >11 tips**: any Faith-PD panel/plot restricted to clusters with >11 tips
  (not the standard >1 singleton filter) is doing so because at smaller sizes — especially
  the sparse 6-10 bin — a single cluster (literally n=1 for some genoflu/genin
  segment×bin combinations) can swing a median and produce an artifact (e.g. Genin2
  briefly reading lower than TreeCluster) that isn't a real pattern. This is a
  data-driven threshold, not a pre-specified one — justify it in the manuscript by
  pointing at `faith-pd_by-cluster-size.csv`/`.png` (the evidence for where the small-n
  instability actually stops), not by asserting 11 as inherently meaningful.
- **genin HA and M** each have only 1 non-singleton cluster (spans the segment), so
  purity/monophyly/PD metrics are trivially defined and are dropped for those two
  segment×metric cells (size/N-clusters panels keep them).
- **Bar charts on a log y-axis, and why `coord_trans()` not `scale_*_log10()`:**
  `geom_col()` always draws bars from 0, which a log scale can't represent. Using
  `scale_y_continuous(trans="log10")` (confirmed empirically, not assumed) makes ggplot
  substitute log10(1)=0 as the bar baseline instead of the panel's real floor — bars for
  values <1 render *downward* from y=1 and bars >1 *upward* from y=1, so bars don't
  share one visual baseline. `coord_trans(y="log10", ylim=c(floor, max))` instead
  transforms only the rendering, leaving `geom_col`'s baseline math in real
  (untransformed) space, so every bar correctly runs from the coordinate system's floor
  up to its value. `coord_trans()` is soft-deprecated for `coord_transform()` as of
  ggplot2 4.0, but the replacement wasn't exported yet as of ggplot2 4.0.2 — swap it in
  once available. A `scale_y_continuous(trans="log1p", ...)` alternative (avoids both
  the baseline-snap and log10's negative-looking low-end tick labels) was tried and
  works, but log1p isn't an intuitive axis for a general readership, hence log10 +
  `coord_trans` + explicit minor gridlines instead.
- **Zero is a real value, not missing data, and is outlined to show it.** E.g. GenoFlu
  genuinely has 0% monophyletic clusters in several segments (matches the manuscript's
  own finding). A zero-height bar is otherwise invisible against the axis and easy to
  misread as an absent segment×method combination (which also happens, deliberately,
  for genin HA/M — see above). The main bar-chart figure outlines only true-zero bars
  (`aes(colour = is_zero)`, mapped per-row, not a fixed `colour=` on every bar) so a
  real zero is visually distinguishable from genuinely absent data.

## Outputs (`YYYY-MM-DD_cluster_metrics_compare/`)

- `cluster-metrics_by-segment_bars_facets_filt.png` — **main comparison figure.** 6
  metrics (N clusters, median/max cluster size, median purity, proportion clusters
  monophyletic, median Faith PD per tip [size >11]), each built as an independent
  ggplot (own y-scale/transform) and combined via `cowplot::plot_grid` into a 3×2 grid
  with one shared legend — not a single `facet_wrap()`. Tried `ggh4x::facetted_pos_scales()`
  first for per-facet log/linear scales within one `facet_wrap()`; its positional
  matching to panels didn't behave as documented (confirmed by disabling it — the
  warnings/dropped bars traced to non-log panels being silently log-transformed too)
  and wasn't worth chasing further once the independent-panels approach worked
  immediately and generalizes trivially to any future per-metric scale change.
- `..._nomonotip.png`, `..._nomonotip_sq.png` — older facet_wrap-based variants,
  pre-dating the independent-panel rebuild above; not kept in sync with it.
- `cluster-metrics_long_filtered.csv` — tidy metric table backing the bar plots, plus
  an `N singletons` metric (not plotted). `N clusters` counts non-singleton clusters only
- `cluster-counts_by-segment.csv` — `N_cluster` (excl. singletons) and `N_single` per
  segment × method, counted from the per-cluster files before the singleton filter.
  Use this for cluster counts in manuscript text. The upstream
  `*_segment-metrics_all-segments.csv` `n_clusters` column counts singletons too
  (e.g. CD-HIT HA 674 = 460 + 214), so don't quote it
- `cluster-metrics_long_informative-only.csv`, `cluster-metrics_method-summary_by-metric.csv`,
  `*_coarse-bins.csv` — **stale**, from the Feb 2026 `compare_clusters.R` and earlier
  coarse-bin runs, and no longer written. Their `n_clusters` and proportions include
  singletons, and genin is missing
- `monophyly_by-cluster-size.png` / `.csv` — **prop. clusters monophyletic vs cluster
  size bin, pooled across segments** (lead figure for the size-controlled comparison).
- `monophyly_by-cluster-size_by-segment.png` / `.csv` — same, faceted by segment (supplement).
- `faith-pd_by-cluster-size.png` / `.csv` — **median Faith's PD (log scale) vs cluster
  size bin, pooled across segments**, the size-controlled PD analogue of the monophyly
  plot above.
- `faith-pd_by-cluster-size_by-segment.png` / `.csv` — same, faceted by segment (supplement).
- `purity_by-cluster-size.png` / `.csv`, `purity_by-cluster-size_by-segment.png` / `.csv`
  — same structure, for median purity.
- `faith-pd_ridgeline_by-method.png`, `purity_ridgeline_by-method.png` — density
  ridgeline per method, all non-singleton clusters pooled across segments (one row per
  method, median marked via ggridges' own quantile line). Density, not raw counts,
  since non-singleton cluster totals differ hugely by method (genoflu ~140, genin ~70,
  cdhit ~3700, treecluster ~4800). Built specifically to check mean vs median for
  purity's summary statistic (see Conventions above) — a dashed mean line was tried and
  removed once both distributions were confirmed skewed.
- `composite_faith-pd_overall-and-by-segment.png`, `composite_monophyly_...`,
  `composite_purity_...` — all three identical in structure: panel A (pooled
  metric-vs-cluster-size plot) on top, panel B (by-segment breakdown) below at
  2.1x A's height, combined via `cowplot::plot_grid` (`align="hv", axis="tblr"`
  so axes line up despite differing legend widths). Faith PD's composite first
  used its ridgeline as panel A in a row layout (ridgeline left/narrow,
  by-segment right/wide), since a ridgeline has no size categories to align
  against B — switched to match monophyly/purity's structure exactly instead;
  the ridgeline remains available as its own standalone file
  (`faith-pd_ridgeline_by-method.png`).

## Key caveats / interpretation

- **treecluster monophyly = 1.0 is an algorithmic artifact, not a result.** This run
  is `tc_maxclade_0.04` (max-clade mode), which defines every cluster as a monophyletic
  clade by construction, so its monophyly metrics are guaranteed ~1.0 and are
  uninformative. The meaningful monophyly contrasts are among genoflu / genin / cdhit.
- **`prop_clusters_monophyletic` is size-confounded.** It averages over cluster *count*,
  which over-weights the many tiny clusters (monophyly is easy at small size). This makes
  cdhit's aggregate bar (~0.44) look deceptively high. The **monophyly-vs-size** plot is
  the honest, size-controlled view: cdhit degrades from ~0.75 (pairs) to ~0.05 (200+ tips),
  far below treecluster at matched sizes. `prop_tips_monophyletic` (tip-weighted) tells the
  same story from the other direction (cdhit ~0.05, since big non-monophyletic clusters hold
  most tips).
- **Per-segment facets are reliable only for cdhit and treecluster** (thousands of clusters).
  genoflu (~140) and genin (~67 total non-singleton clusters) are sparse per segment, so their
  faceted curves are jagged/broken and should be read as suggestive only; use the pooled figure
  for them.
- **genoflu and genin metrics are computed on a pruned tree; cdhit and treecluster are
  not — the four methods aren't evaluated on strictly the same tree per segment.**
  GenoFlu and Genin2 leave some sequences unassigned (NA, or "Unassigned" for genin);
  their upstream scripts drop those tips from `cluster_df` and then `keep.tip()`-prune
  the tree to assigned tips only *before* computing MRCA, `desc_tip_n`, purity,
  monophyly type, and `faith_pd` — unassigned tips are excluded from the phylogeny for
  that evaluation, not imputed and not counted as impurity. CD-HIT and TreeCluster
  assign every tip to some cluster (down to singletons), so there's no unassigned
  category and no pruning — their scripts `stop()` if any cluster assignment is
  missing. Practically: a genoflu/genin cluster's purity, `mrca_desc_n`, and
  `faith_pd` are relative to a smaller, assignment-dependent tree, while cdhit/
  treecluster's are relative to the full segment tree. This doesn't bias any single
  method's own metrics, but it means segment-level cross-method comparisons aren't
  on an identical tree when genoflu/genin non-assignment is substantial for that
  segment (e.g. PB2, PB1, NP — see Fig. 3 assignment-rate results).
- **Faith PD is not a construction artifact for treecluster the way monophyly is.**
  max-clade guarantees monophyly==1.0, but a monophyletic clade can still span a lot of
  branch length if its members are deeply diverged, so treecluster's PD-vs-size curve is a
  genuine (not trivial) result and the most direct test of the manuscript's "minimum span to
  include all labels" comparison: at matched cluster size, genoflu/genin/cdhit clusters are
  expected to need more tree distance to connect their members than treecluster's monophyletic
  clusters do.
