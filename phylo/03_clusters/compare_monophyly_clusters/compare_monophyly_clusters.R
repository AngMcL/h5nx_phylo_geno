# Compare cluster metrics across methods + segments
# - single source of truth: per-cluster metrics files
# - ALL metrics recomputed after excluding singletons (n_tips > 1) for consistency
# - robust CSV ingest (keeps "NA" segment as a real string)
# - side-by-side bars + facet by metric
# - N clusters and cluster-size metrics shown on log scale

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(fs)
  library(purrr)
  library(ggridges)
  library(cowplot)
})

# standing project ggplot theme (matches every other script in this pipeline):
# grey95 panel background, no minor gridlines (major kept, visible against the
# grey background), black axis lines, bold text, angled x-axis by default.
plot.theme <- theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect("grey95"),
                    axis.line = element_line(colour = "black"),
                    text = element_text(size = 10, face = "bold"),
                    legend.text = element_text(size = 8),
                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# INPUTS (per-cluster metrics: one row per cluster)
in_files_cluster <- tibble::tribble(
  ~source,        ~path,
  "genoflu",      "../genoflu_clades/2026-08-20_out/genoflu_cluster-metrics_all-segments.csv",
  "genin",        "../genin_clades/2026-08-20_out/genin_cluster-metrics_all-segments.csv",
  "cdhit",        "../cdhit/cluster_eval/cdhit_cluster-metrics_all-segments.csv",
  "treecluster",  "../treecluster/02_tc_trees/cluster_eval/treecluster_cluster-metrics_all-segments.csv"
)

out_dir <- paste0(Sys.Date(), "_cluster_metrics_compare")  # dated; runs before 2026-09-22 are in cluster_metrics_compare/
dir_create(out_dir)

segments <- c("PB2","PB1","PA","HA","NP","NA","M","NS")

# READ PER-CLUSTER METRICS
read_clusters <- function(path, source) {
  if (!file.exists(path)) stop(path, " does not exist.")
  readr::read_csv(
    path,
    show_col_types = FALSE,
    na = c("", "NaN"),                 # IMPORTANT: do NOT treat "NA" as missing
    col_types = cols(.default = col_character())) %>%
    mutate(
      source = source,
      method = recode(method, "genin_cluster" = "genin"),
      n_tips = as.integer(n_tips),
      purity = suppressWarnings(as.numeric(purity)),
      faith_pd = suppressWarnings(as.numeric(faith_pd)))
}

clusters_raw <- purrr::pmap_dfr(in_files_cluster, read_clusters) %>%
  mutate(
    segment = factor(ifelse(is.na(segment), "NA", segment), levels = segments),
    method  = factor(method, levels = c("genoflu", "genin", "cdhit", "treecluster")))

# COUNT CLUSTERS vs SINGLETONS (before the singleton filter below)
# N_cluster excludes singletons; N_single reported separately so the two are
# never conflated (the old segment-metrics n_clusters summed both, e.g. CD-HIT
# HA 674 = 460 clusters + 214 singletons)
cluster_counts <- clusters_raw %>%
  group_by(segment, method) %>%
  summarise(N_cluster = sum(n_tips > 1),
            N_single  = sum(n_tips == 1),
            .groups = "drop") %>%
  arrange(method, segment)
write_csv(cluster_counts, file.path(out_dir, "cluster-counts_by-segment.csv"))

# EXCLUDE SINGLETONS everywhere (n_tips > 1)
clusters <- clusters_raw %>% filter(n_tips > 1)

# RECOMPUTE ALL METRICS from non-singleton clusters
# Values are RAW here (not pre-logged) — the log-scale metrics get an actual
# log10 y-axis via scale_y_log10() at plot time instead (see the split
# log/linear plot construction below), so the facet label just reads
# "N clusters" etc., not "log10 (N clusters)".
# Distributions of both purity and Faith PD are visibly non-normal (checked via
# the ridgeline plots), so both use MEDIAN, not mean, throughout — mean cluster
# size and prop. tips monophyletic (redundant with prop. clusters monophyletic
# once things are median-based) are dropped rather than kept alongside.
metric_levels <- c("N clusters", "Median cluster size", "Maximum cluster size",
                   "Median purity", "Proportion clusters monophyletic",
                   "Median Faith's PD [size >11] (s/s)")

metrics_summary <- clusters %>%
  group_by(segment, method) %>%
  summarise(
    `N clusters`                    = n(),
    `Median cluster size`           = median(n_tips),
    `Maximum cluster size`          = max(n_tips),
    `Median purity`                 = median(purity, na.rm = TRUE),
    `Proportion clusters monophyletic` = mean(type == "monophyletic"),
    # >11 tips only: at smaller sizes (esp. the sparse 6-10 bin) a single
    # cluster can swing the median and make a method look lower than
    # TreeCluster when it isn't a real pattern — see faith-pd_by-cluster-size.csv
    # Raw (not per-tip): per-tip obscures how PD is actually distributed and
    # how much large clusters contribute; kept raw and plotted with a log1p
    # y-axis instead (see metric_specs below) so the wide dynamic range is
    # still readable without per-tip normalization.
    `Median Faith's PD [size >11] (s/s)` = median(faith_pd[n_tips > 11], na.rm = TRUE),
    .groups = "drop")

metrics_long_f <- metrics_summary %>%
  pivot_longer(cols = all_of(metric_levels), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_levels)) %>%
  # genin HA and M have only 1 (non-singleton) cluster spanning the segment —
  # purity/monophyly/PD metrics are trivially defined (a single genome-spanning
  # cluster), so drop them for those panels
  filter(!(method == "genin" & segment %in% c("HA", "M") &
             metric %in% c("Median purity", "Proportion clusters monophyletic",
                            "Median Faith's PD [size >11] (s/s)"))) %>%
  # zero is a real value here (e.g. GenoFlu genuinely has 0% monophyletic
  # clusters in several segments), not missing data — outline just the zero
  # bars so they're visible against the axis instead of reading as absent
  mutate(is_zero = value == 0)

# -----------------------------
# COLORS: Paul Tol "vibrant" qualitative palette (colorblind-safe), fixed one
# color per method. Segment colors elsewhere in this pipeline (geno_reassort.Rmd
# etc.) already use Tol's "bright" palette, so these are deliberately drawn from
# the sister "vibrant" set instead, avoiding overlap between the two encodings.
# -----------------------------
methods_all <- levels(metrics_long_f$method)

tol_vibrant_method <- c(genoflu = "#0077BB", genin = "#EE7733",
                        cdhit = "#009988", treecluster = "#CC3311")
method_cols <- tol_vibrant_method[methods_all]

# Each metric is its OWN independent ggplot (own y scale/transform), combined
# via plot_grid into a 3x2 grid — NOT one facet_wrap with a shared/mixed scale.
# Tried two facet_wrap-based approaches first (ggh4x::facetted_pos_scales() for
# a single unified facet_wrap: its positional matching to panels didn't line up
# as documented, confirmed by disabling it — warnings/dropped bars traced to
# non-log panels being silently log-transformed too; then a split log-row/
# linear-row pair of facet_wraps), but different metrics keep needing different
# transforms (log10 for counts/sizes/Faith PD, identity for proportions).
# Building one independent panel per metric makes each metric's transform a
# one-line choice instead of fighting facet_wrap/ggh4x to keep several
# transforms straight in one call.
# Faith PD tried log1p first (avoids log10's negative-looking tick labels for
# values <1) but that scale isn't intuitive to a general audience either — back
# to log10, with faint minor gridlines at the standard 2-9 sub-decade positions
# (scales::minor_breaks_log()) so the log nature of the axis reads visually
# instead of relying on the reader to parse tick labels alone.
metric_specs <- tibble::tribble(
  ~metric,                                ~trans,    ~log_minor,
  "N clusters",                           "log10",   FALSE,
  "Median cluster size",                  "log10",   FALSE,
  "Maximum cluster size",                 "log10",   FALSE,
  "Median purity",                        "identity", FALSE,
  "Proportion clusters monophyletic",     "identity", FALSE,
  "Median Faith's PD [size >11] (s/s)",   "log10",   TRUE
)

zero_outline_scale <- scale_colour_manual(values = c(`FALSE` = NA, `TRUE` = "black"), guide = "none")

make_metric_panel <- function(m, trans, log_minor = FALSE) {
  p <- ggplot(metrics_long_f %>% filter(metric == m),
             aes(x = segment, y = value, fill = method, colour = is_zero)) +
    geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.75,
            linewidth = 0.4)

  # geom_col draws bars from 0. scale_y_continuous(trans="log10") can't
  # represent 0, and (confirmed empirically) substitutes log10(1)=0 as the
  # bar baseline instead of the panel's actual floor — bars for values <1 then
  # render downward from y=1 and bars >1 upward from y=1, instead of every bar
  # sharing one common bottom. coord_trans() transforms the rendering only,
  # leaving geom_col's own baseline math in real (untransformed) space, so
  # every bar correctly runs from the coordinate system's floor up to its
  # value. (coord_trans() is soft-deprecated for coord_transform() as of
  # ggplot2 4.0, but the replacement isn't exported yet in this version.)
  p <- if (log_minor) {
    max_val <- max(metrics_long_f$value[metrics_long_f$metric == m], na.rm = TRUE)
    p + scale_y_continuous(breaks = c(0.01, 0.1, 1, 10),
                           minor_breaks = scales::minor_breaks_log()) +
      coord_trans(y = "log10", ylim = c(0.01, max_val * 1.1))
  } else {
    p + scale_y_continuous(trans = trans)
  }

  p +
    zero_outline_scale +
    scale_fill_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
    labs(x = NULL, y = m, fill = "Cluster method") +
    theme_bw(base_size = 12) +
    plot.theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = if (log_minor) element_line(colour = "grey85", linewidth = 0.15) else element_blank(),
      axis.title.y = element_text(size = 10, face = "bold"),
      legend.position = "top",
      legend.key.size = unit(1, "lines")
    ) +
    guides(fill = guide_legend(nrow = 1, title.position = "left"))
}

metric_panels <- purrr::pmap(metric_specs, function(metric, trans, log_minor) make_metric_panel(metric, trans, log_minor))

# shared legend: extract once, then strip legends from every panel before
# combining, so the composite has a single legend instead of one per panel
metrics_shared_legend <- get_legend(metric_panels[[1]])
metric_panels_nolegend <- lapply(metric_panels, function(p) p + theme(legend.position = "none"))
p_metrics_grid <- plot_grid(plotlist = metric_panels_nolegend, ncol = 2, nrow = 3,
                            align = "hv", axis = "tblr")
p_metrics <- plot_grid(metrics_shared_legend, p_metrics_grid, ncol = 1,
                       rel_heights = c(0.08, 1))

ggsave(
  filename = file.path(out_dir, "cluster-metrics_by-segment_bars_facets_filt.png"),
  plot = p_metrics,
  width = 8.6, height = 8, dpi = 300
)


## restrict to fewer metrics (drop prop. tips monophyletic)
p_metrics2 <- metrics_long_f %>%
  filter(metric != "Prop. tips monophyletic") %>%
  ggplot(aes(x = segment, y = value, fill = method)) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.75) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_fill_manual(values = method_cols,
                    labels = c(genoflu = "GenoFlu", genin = "Genin2",
                               cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  labs(x = NULL, y = NULL, fill = "Method") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "top",
    legend.key.size = unit(1, "lines"),
    text = element_text(size = 13)
  ) +
  guides(fill = guide_legend(nrow = 1, title.position = "left"))

ggsave(
  filename = file.path(out_dir, "cluster-metrics_by-segment_bars_facets_filt_nomonotip.png"),
  plot = p_metrics2,
  width = 7.5, height = 14, dpi = 300
)

## square-ish version (drop prop. tips monophyletic) — 2 rows
p_metrics2h <- metrics_long_f %>%
  filter(metric != "Prop. tips monophyletic") %>%
  ggplot(aes(x = segment, y = value, fill = method)) +
  geom_col(position = position_dodge2(width = 0.25, preserve = "single"), width = 0.9) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2, strip.position = "left") +
  scale_fill_manual(values = method_cols,
                    labels = c(genoflu = "GenoFlu", genin = "Genin2",
                               cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  labs(x = NULL, y = NULL, fill = "Method") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10),
    legend.position = "top",
    legend.key.size = unit(1, "lines"),
    text = element_text(size = 13)) +
  guides(fill = guide_legend(nrow = 1, title.position = "left"))

ggsave(filename = file.path(out_dir, "cluster-metrics_by-segment_bars_facets_filt_nomonotip_sq.png"),
  plot = p_metrics2h,
  width = 8, height = 5, dpi = 300)


#writedata object
# N singletons appended as its own metric (not plotted) alongside N clusters,
# which counts non-singleton clusters only
stopifnot(all.equal(metrics_summary$`N clusters`,
                    cluster_counts$N_cluster[match(paste(metrics_summary$segment, metrics_summary$method),
                                                   paste(cluster_counts$segment, cluster_counts$method))]))
metrics_long_out <- bind_rows(
  metrics_long_f %>% mutate(metric = as.character(metric)),
  cluster_counts %>% transmute(segment, method, metric = "N singletons",
                               value = N_single, is_zero = N_single == 0)) %>%
  arrange(method, segment)
write_csv(metrics_long_out, file.path(out_dir, "cluster-metrics_long_filtered.csv"))

# -----------------------------
# MONOPHYLY vs CLUSTER SIZE (size-controlled comparison)
# prop_clusters_monophyletic averaged over cluster COUNT is dominated by the many
# tiny clusters, so it over-states cdhit. Stratifying by size shows cdhit degrading
# with size while TreeCluster (max-clade) is monophyletic by construction (flat ~1).
# -----------------------------
# Same fine high-end bins used for both the pooled plot and the per-segment facets
# (kept as one scheme throughout, not coarse-vs-fine, so panel A/B x-axes always
# match in composites): genoflu/genin clusters get large enough that collapsing
# everything above 200 into one bin loses real resolution.
size_breaks_seg <- c(1, 5, 10, 50, 100, 200, 500, 1000, 10000, Inf)
size_labs_seg   <- c("2-5", "6-10", "11-50", "51-100", "100-200",
                     "200-500", "500-1000", "1000-10000", "10000+")

# pooled across segments
mono_by_size <- clusters %>%
  mutate(size_bin = cut(n_tips, breaks = size_breaks_seg, labels = size_labs_seg)) %>%
  group_by(method, size_bin) %>%
  summarise(n_clusters = n(),
            prop_monophyletic = mean(type == "monophyletic"),
            .groups = "drop")

# per segment
mono_by_size_seg <- clusters %>%
  mutate(size_bin = cut(n_tips, breaks = size_breaks_seg, labels = size_labs_seg)) %>%
  group_by(segment, method, size_bin) %>%
  summarise(n_clusters = n(),
            prop_monophyletic = mean(type == "monophyletic"),
            .groups = "drop")

write_csv(mono_by_size, file.path(out_dir, "monophyly_by-cluster-size.csv"))
write_csv(mono_by_size_seg, file.path(out_dir, "monophyly_by-cluster-size_by-segment.csv"))

# overall per method (all segments and sizes pooled; non-singleton clusters only).
# Genin2 pools 6 segments (no HA/M)
mono_overall <- clusters %>%
  group_by(method) %>%
  summarise(n_clusters = n(),
            n_monophyletic = sum(type == "monophyletic"),
            prop_monophyletic = n_monophyletic / n_clusters,
            .groups = "drop") %>%
  arrange(desc(prop_monophyletic))

write_csv(mono_overall, file.path(out_dir, "monophyly_overall_by-method.csv"))

p_mono_size <- ggplot(mono_by_size,
                      aes(x = size_bin, y = prop_monophyletic,
                          colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_clusters)) +
  scale_colour_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_size_area(max_size = 6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Cluster size (n tips, singletons excluded)",
       y = "Prop. clusters monophyletic",
       colour = "Cluster method", size = "N clusters") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(filename = file.path(out_dir, "monophyly_by-cluster-size.png"),
       plot = p_mono_size, width = 7, height = 4.5, dpi = 300)

# per-segment facets (supplement) — sparse for genin/genoflu given few clusters
p_mono_size_seg <- ggplot(mono_by_size_seg,
                          aes(x = size_bin, y = prop_monophyletic,
                              colour = method, group = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_clusters)) +
  facet_wrap(~ segment, ncol = 4) +
  scale_colour_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_size_area(max_size = 4, breaks = c(5, 10, 50, 100, 150, 200)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Cluster size (n tips, singletons excluded)",
       y = "Prop. clusters monophyletic",
       colour = "Cluster method", size = "N clusters") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

ggsave(filename = file.path(out_dir, "monophyly_by-cluster-size_by-segment.png"),
       plot = p_mono_size_seg, width = 14, height = 6.5, dpi = 300)

# -----------------------------
# FAITH'S PD vs CLUSTER SIZE (size-controlled comparison)
# Faith PD (total branch length of the minimal subtree spanning a cluster's tips)
# scales with cluster size even for a well-behaved (monophyletic) cluster, so raw
# PD is only comparable within a size bin, not across bins or methods pooled.
# A monophyletic TreeCluster cluster's PD reflects only its clade depth; a
# non-monophyletic GenoFlu/Genin2/CD-HIT cluster of the same size additionally
# pays for the branch length connecting its scattered members, so at matched
# size its PD is expected to exceed TreeCluster's.
# -----------------------------
pd_by_size <- clusters %>%
  mutate(size_bin = cut(n_tips, breaks = size_breaks_seg, labels = size_labs_seg)) %>%
  group_by(method, size_bin) %>%
  summarise(n_clusters = n(),
            median_faith_pd = median(faith_pd, na.rm = TRUE),
            .groups = "drop")

pd_by_size_seg <- clusters %>%
  mutate(size_bin = cut(n_tips, breaks = size_breaks_seg, labels = size_labs_seg)) %>%
  group_by(segment, method, size_bin) %>%
  summarise(n_clusters = n(),
            median_faith_pd = median(faith_pd, na.rm = TRUE),
            .groups = "drop")

write_csv(pd_by_size, file.path(out_dir, "faith-pd_by-cluster-size.csv"))
write_csv(pd_by_size_seg, file.path(out_dir, "faith-pd_by-cluster-size_by-segment.csv"))

p_pd_size <- ggplot(pd_by_size,
                    aes(x = size_bin, y = median_faith_pd,
                        colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_clusters)) +
  scale_colour_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_size_area(max_size = 6) +
  scale_y_log10() +
  labs(x = "Cluster size (n tips, singletons excluded)",
       y = "Median Faith's PD (subs/site)",
       colour = "Cluster method", size = "N clusters") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(filename = file.path(out_dir, "faith-pd_by-cluster-size.png"),
       plot = p_pd_size, width = 7, height = 4.5, dpi = 300)

# per-segment facets (supplement) — sparse for genin/genoflu given few clusters
p_pd_size_seg <- ggplot(pd_by_size_seg,
                        aes(x = size_bin, y = median_faith_pd,
                            colour = method, group = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_clusters)) +
  facet_wrap(~ segment, ncol = 4) +
  scale_colour_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_size_area(max_size = 4, breaks = c(5, 10, 50, 100, 150, 200)) +
  scale_y_log10() +
  labs(x = "Cluster size (n tips, singletons excluded)",
       y = "Median Faith's PD (subs/site)",
       colour = "Cluster method", size = "N clusters") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

ggsave(filename = file.path(out_dir, "faith-pd_by-cluster-size_by-segment.png"),
       plot = p_pd_size_seg, width = 14, height = 6.5, dpi = 300)

# -----------------------------
# FAITH'S PD: OVERALL DISTRIBUTION PER METHOD (ridgeline)
# Density (not raw counts) so different cluster totals per method don't distort
# shape comparison. All non-singleton clusters (no additional size filter).
# Median marked on each ridge via a quantile line (ggridges' own quantile_lines,
# quantiles=2 -> the 50th percentile) — median only, not mean: both this and the
# purity ridgeline below show visibly non-normal, skewed distributions, so mean
# isn't a representative summary here and was dropped after inspecting both.
# Single panel, all 4 methods on one shared x-axis, for direct side-by-side
# comparability. Linear x-axis, not log: cdhit has real non-singleton clusters
# with faith_pd==0, which a log scale can't display.
# -----------------------------
pd_ridge_data <- clusters

p_pd_ridge <- ggplot(pd_ridge_data, aes(x = faith_pd, y = method, fill = method)) +
  geom_density_ridges(alpha = 0.8, scale = 2, rel_min_height = 0.005,
                      colour = "white", linewidth = 0.3,
                      quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = method_cols,
                    labels = c(genoflu = "GenoFlu", genin = "Genin2",
                               cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_y_discrete(labels = c(genoflu = "GenoFlu", genin = "Genin2",
                              cdhit = "CD-HIT", treecluster = "TreeCluster"),
                   expand = expansion(add = c(0.3, 1.0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Faith's PD (s/s)", y = NULL) +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "none")

ggsave(filename = file.path(out_dir, "faith-pd_ridgeline_by-method.png"),
       plot = p_pd_ridge, width = 7, height = 2.25, dpi = 300)

# -----------------------------
# PURITY: OVERALL DISTRIBUTION PER METHOD (ridgeline)
# Same structure as the Faith PD ridgeline above. Built to check whether mean or
# median better represents purity's distribution shape; visibly skewed/non-normal
# on inspection, so median only (mean dropped), matching the Faith PD decision.
# -----------------------------
purity_ridge_data <- clusters

p_purity_ridge <- ggplot(purity_ridge_data, aes(x = purity, y = method, fill = method)) +
  geom_density_ridges(alpha = 0.8, scale = 2, rel_min_height = 0.005,
                      colour = "white", linewidth = 0.3,
                      quantile_lines = TRUE, quantiles = 2) +
  scale_fill_manual(values = method_cols,
                    labels = c(genoflu = "GenoFlu", genin = "Genin2",
                               cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_y_discrete(labels = c(genoflu = "GenoFlu", genin = "Genin2",
                              cdhit = "CD-HIT", treecluster = "TreeCluster"),
                   expand = expansion(add = c(0.3, 1.0))) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Purity", y = NULL) +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "none")

ggsave(filename = file.path(out_dir, "purity_ridgeline_by-method.png"),
       plot = p_purity_ridge, width = 7, height = 2.25, dpi = 300)

# -----------------------------
# CLUSTER PURITY vs CLUSTER SIZE (size-controlled comparison)
# Same structure as the monophyly-vs-size and Faith-PD-vs-size sections above: mean
# purity per size bin per method, pooled and per-segment. Mean (not median) to match
# "Mean purity" as already named in the main by-segment bar chart panel.
# -----------------------------
purity_by_size <- clusters %>%
  mutate(size_bin = cut(n_tips, breaks = size_breaks_seg, labels = size_labs_seg)) %>%
  group_by(method, size_bin) %>%
  summarise(n_clusters = n(),
            median_purity = median(purity, na.rm = TRUE),
            .groups = "drop")

purity_by_size_seg <- clusters %>%
  mutate(size_bin = cut(n_tips, breaks = size_breaks_seg, labels = size_labs_seg)) %>%
  group_by(segment, method, size_bin) %>%
  summarise(n_clusters = n(),
            median_purity = median(purity, na.rm = TRUE),
            .groups = "drop")

write_csv(purity_by_size, file.path(out_dir, "purity_by-cluster-size.csv"))
write_csv(purity_by_size_seg, file.path(out_dir, "purity_by-cluster-size_by-segment.csv"))

p_purity_size <- ggplot(purity_by_size,
                        aes(x = size_bin, y = median_purity,
                            colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_clusters)) +
  scale_colour_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_size_area(max_size = 6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Cluster size (n tips, singletons excluded)",
       y = "Median purity",
       colour = "Cluster method", size = "N clusters") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(filename = file.path(out_dir, "purity_by-cluster-size.png"),
       plot = p_purity_size, width = 7, height = 4.5, dpi = 300)

# per-segment facets (supplement) — sparse for genin/genoflu given few clusters
p_purity_size_seg <- ggplot(purity_by_size_seg,
                            aes(x = size_bin, y = median_purity,
                                colour = method, group = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = n_clusters)) +
  facet_wrap(~ segment, ncol = 4) +
  scale_colour_manual(values = method_cols,
                      labels = c(genoflu = "GenoFlu", genin = "Genin2",
                                 cdhit = "CD-HIT", treecluster = "TreeCluster")) +
  scale_size_area(max_size = 4, breaks = c(5, 10, 50, 100, 150, 200)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Cluster size (n tips, singletons excluded)",
       y = "Median purity",
       colour = "Cluster method", size = "N clusters") +
  theme_bw(base_size = 12) +
  plot.theme +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

ggsave(filename = file.path(out_dir, "purity_by-cluster-size_by-segment.png"),
       plot = p_purity_size_seg, width = 14, height = 6.5, dpi = 300)

# -----------------------------
# COMPOSITES
# All three: same layout/composition — A (pooled metric-vs-cluster-size, top) +
# B (by-segment breakdown, bottom) at 2.1x A's height (more headroom for A's
# legend than the first pass, which was clipping the legend title). Panel A uses
# the same fine size bins as B (see MONOPHYLY/PURITY/FAITH-PD-vs-CLUSTER-SIZE
# sections above), so they line up directly without any rebinning. Faith PD's
# panel A was a ridgeline in a row layout at first (paired with panel B side by
# side, since a ridgeline has no size categories to match against B) — switched
# to the same pooled-vs-size plot type and top/bottom structure as monophyly/
# purity for consistency across all three composites; the ridgeline still exists
# standalone (faith-pd_ridgeline_by-method.png) for anyone who wants it.
# align="hv" (cowplot) aligns the actual panel areas, not just the outer image
# bounds, so axes line up cleanly between A and B despite differing legend/strip
# widths.
# -----------------------------
composite_faithpd <- plot_grid(p_pd_size, p_pd_size_seg, ncol = 1,
                               rel_heights = c(1.4, 2.1), align = "hv", axis = "tblr",
                               labels = c("A", "B"), label_size = 16)
ggsave(filename = file.path(out_dir, "composite_faith-pd_overall-and-by-segment.png"),
       plot = composite_faithpd, width = 10, height = 10, dpi = 300)

composite_monophyly <- plot_grid(p_mono_size, p_mono_size_seg, ncol = 1,
                                 rel_heights = c(1.4, 2.1), align = "hv", axis = "tblr",
                                 labels = c("A", "B"), label_size = 16)
ggsave(filename = file.path(out_dir, "composite_monophyly_overall-and-by-segment.png"),
       plot = composite_monophyly, width = 10, height = 10, dpi = 300)

composite_purity <- plot_grid(p_purity_size, p_purity_size_seg, ncol = 1,
                              rel_heights = c(1.4, 2.1), align = "hv", axis = "tblr",
                              labels = c("A", "B"), label_size = 16)
ggsave(filename = file.path(out_dir, "composite_purity_overall-and-by-segment.png"),
       plot = composite_purity, width = 10, height = 10, dpi = 300)
