# Controls for the cluster-ID span vs number-of-clusters screen.
#
# The screen (subtree_cluster_matrix.Rmd, geno_tc_cluster_spread.csv) flags a
# genotype x segment when the genotype's member tc_cluster IDs are scattered
# across a wide swath of that segment's tree: ratio = span / n_clusters, where
# span = max - min + 1 of the NUMERIC cluster IDs and singleton "-s" buckets are
# excluded from both. Cluster IDs are numbered bottom -> top of each segment's
# ladderized tree, so ratio ~1 means a contiguous block and ratio >> 1 means the
# genotype's sequences sit in unrelated parts of the tree.
#
# This script reports the full per-segment profile for named control genotypes so
# the metric can be read against a known answer:
#   positive control  B3.5 -- elevated on 5 segments at once (PA/HA/NA/M/NS), the
#                             strongest multi-segment candidate in the screen
#   negative control  A6   -- 315 sequences that all carry ONE tc_geno (TC3173),
#                             so span = n_clusters = ratio = 1 in every segment;
#                             a clonal genotype must score 1, and it does
# Input is geno_tc_cluster_summary.csv (every genotype, no N threshold, so the
# clonal controls are present) rather than geno_tc_cluster_spread.csv, which
# keeps only rows with span >= 3 and would therefore drop A6 entirely.

source("_common.R")

pos.controls <- c("B3.5")
neg.controls <- c("A6")
flag.ratio   <- 2      # ratio at/above which the screen calls a segment scattered
top.n        <- 15     # how many ranked screen rows to list for context

summary.in <- newest_in(scm.dir, "_subtree_cluster_matrix_out$", "geno_tc_cluster_summary.csv")
spread.in  <- newest_in(scm.dir, "_subtree_cluster_matrix_out$", "geno_tc_cluster_spread.csv")

geno_summary <- read.csv(summary.in, check.names = FALSE)
# na.strings = "" so the spread table's `segment` value for neuraminidase stays
# the string "NA" instead of being read as missing
geno_spread  <- read.csv(spread.in, check.names = FALSE, na.strings = "")

# wide (one row per genotype, 4 columns per segment) -> long (one row per
# genotype x segment). Column names are "{SEG}_{field}" and the NA segment makes
# "NA_n_clusters" etc., so split on the LAST underscore-delimited field name
# rather than the first token.
fields <- c("n_clusters", "id_range", "span", "ratio")
to_long <- function(df) {
  df %>%
    pivot_longer(cols = matches(paste0("_(", paste(fields, collapse = "|"), ")$")),
                 names_to = "key", values_to = "value", values_transform = as.character) %>%
    mutate(field   = str_extract(key, paste0("(", paste(fields, collapse = "|"), ")$")),
           segment = str_remove(key, paste0("_", field, "$"))) %>%
    select(-key) %>%
    pivot_wider(names_from = field, values_from = value) %>%
    mutate(across(c(n_clusters, span, ratio), as.numeric),
           segment = factor(segment, levels = seg.order)) %>%
    filter(!is.na(segment)) %>%
    arrange(scheme, geno, segment)
}

long <- to_long(geno_summary)

controls <- long %>%
  filter(geno %in% c(pos.controls, neg.controls)) %>%
  mutate(control = case_when(geno %in% pos.controls ~ "positive",
                             geno %in% neg.controls ~ "negative"),
         scattered = !is.na(ratio) & ratio >= flag.ratio) %>%
  relocate(control, .after = geno)

write.csv(controls, paste0(f.out, "control_genos_spread.csv"), row.names = FALSE)

cat("\n-- control genotype per-segment profile --\n")
controls %>%
  select(control, geno, n_seqs, n_tc_geno, segment, n_clusters, id_range, span, ratio, scattered) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nsegments at/above ratio", flag.ratio, "per control:\n")
controls %>%
  group_by(control, geno, n_seqs, n_tc_geno) %>%
  summarise(n_segments_scattered = sum(scattered),
            segments = paste(segment[scattered], collapse = ";"),
            max_ratio = max(ratio, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n-- top", top.n, "screen rows overall, for context --\n")
geno_spread %>%
  arrange(desc(ratio)) %>% slice_head(n = top.n) %>%
  select(scheme, geno, n_seqs, segment, n_numeric_clusters, min_id, max_id, span, ratio) %>%
  as.data.frame() %>% print(row.names = FALSE)

# per-segment ratio for each control, dashed line at the flag threshold
p <- ggplot(controls, aes(x = segment, y = ratio, fill = segment)) +
  geom_hline(yintercept = flag.ratio, linetype = "dashed", colour = "grey40") +
  geom_col(colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = seg.colors, guide = "none") +
  facet_wrap(~paste0(geno, " (", control, " control, n=", n_seqs, ")"), ncol = 1) +
  plot.theme +
  labs(x = "segment", y = "cluster ID span / number of clusters")
ggsave(paste0(f.out, "control_genos_ratio_by_segment.tiff"),
       p, width = 6, height = 6, units = "in", dpi = 300)

cat("\nwrote", f.out, "control_genos_spread.csv + control_genos_ratio_by_segment.tiff\n")
