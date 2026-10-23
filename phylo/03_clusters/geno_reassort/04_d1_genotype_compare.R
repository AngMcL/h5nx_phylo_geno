# Which segments separate the D1 genotypes, and is the split carried by PA and NP?
#
# The action item was framed as D1.1 vs D1.8, but D1.8 is NOT in this dataset:
# GenoFlu assigns only D1.1 (3,525 sequences), D1.3 (352) and D1.2 (5), and the
# data end 2025-11-30. So the comparison runs on the D1 genotypes that exist,
# with D1.1 vs D1.3 as the main contrast, and asks the same underlying question:
# which segments carry the divergence between them, and is it PA and NP?
#
# Three views, per segment:
#   1. the GenoFlu feature itself -- which genoflu_cluster{SEG}_filled values each
#      D1 genotype carries, i.e. which segments GenoFlu uses to tell them apart
#   2. tc_cluster composition -- how many TreeCluster clusters each genotype spans
#      on that segment, their ID range, span and span/n ratio
#   3. patristic distance -- within each genotype, and BETWEEN the two genotypes,
#      on that segment's ML tree. between/within is the per-segment separation:
#      high on the segments that actually distinguish the genotypes, ~1 on the
#      segments they share.
# Distances are RAW branch-length units (subs/site), deliberately not normalised
# (see geno_reassort.Rmd for why).

source("_common.R")

focal.genos  <- c("D1.1", "D1.3", "D1.2")   # compared pairwise in this order
pair.a       <- "D1.1"                       # main contrast
pair.b       <- "D1.3"
min.n        <- 2        # genotypes with fewer sequences are reported but not measured
max.tips     <- 1500     # per-genotype cap on tips entering the distance matrix
dist.seed    <- 1

meta <- read_meta()

d1 <- meta %>% filter(genoflu_geno %in% focal.genos)

cat("\nD1 genotype sizes:\n")
d1 %>% group_by(genoflu_geno) %>%
  summarise(n = n(), first_date = min(coll_lb), last_date = max(coll_lb),
            n_tc_geno = n_distinct(tc_geno),
            countries = paste(sort(unique(country)), collapse = ";"),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# -- 1. the GenoFlu feature per segment --
gf.cols <- paste0("genoflu_cluster", seg.order, "_filled")

gf_long <- map_dfr(seg.order, function(s) {
  col <- paste0("genoflu_cluster", s, "_filled")
  d1 %>% count(geno = genoflu_geno, value = .data[[col]], name = "n") %>%
    group_by(geno) %>%
    # `values` carries the per-value counts for reading; `plain` is the bare set
    # of values, which is what the differs test has to compare -- comparing the
    # counted strings would always differ simply because the genotypes differ in
    # size (3,525 vs 352 sequences)
    summarise(values = paste0(value, "(", n, ")", collapse = ";"),
              plain  = paste(sort(unique(value)), collapse = ";"), .groups = "drop") %>%
    mutate(segment = s)
})

gf_feature <- gf_long %>%
  select(segment, geno, values) %>%
  pivot_wider(names_from = geno, values_from = values) %>%
  left_join(gf_long %>% select(segment, geno, plain) %>%
              pivot_wider(names_from = geno, values_from = plain,
                          names_prefix = "plain_"),
            by = "segment") %>%
  mutate(segment = factor(segment, levels = seg.order),
         differs = .data[[paste0("plain_", pair.a)]] !=
                   .data[[paste0("plain_", pair.b)]]) %>%
  select(-starts_with("plain_")) %>%
  arrange(segment)

write.csv(gf_feature, paste0(f.out, "d1_genoflu_feature_by_segment.csv"), row.names = FALSE)
cat("\n-- GenoFlu segment features (differs = ", pair.a, " vs ", pair.b, ") --\n", sep = "")
as.data.frame(gf_feature) %>% print(row.names = FALSE)

# -- 2. tc_cluster composition per segment --
tc_comp <- map_dfr(seg.order, function(s) {
  col <- paste0("tc_cluster", s)
  d1 %>%
    group_by(geno = genoflu_geno) %>%
    summarise(nums = list(cluster_nums(.data[[col]])),
              has_singleton = any(str_detect(.data[[col]], "-s$")),
              clusters = paste0(names(sort(table(.data[[col]]), decreasing = TRUE)),
                                "(", sort(table(.data[[col]]), decreasing = TRUE), ")",
                                collapse = ";"),
              .groups = "drop") %>%
    mutate(segment = s,
           n_clusters = lengths(nums),
           min_id = map_dbl(nums, ~ if (length(.x)) min(.x) else NA_real_),
           max_id = map_dbl(nums, ~ if (length(.x)) max(.x) else NA_real_),
           span   = max_id - min_id + 1,
           ratio  = round(span / n_clusters, 1)) %>%
    select(-nums)
}) %>%
  mutate(segment = factor(segment, levels = seg.order)) %>%
  select(geno, segment, n_clusters, min_id, max_id, span, ratio, has_singleton, clusters) %>%
  arrange(geno, segment)

write.csv(tc_comp, paste0(f.out, "d1_tc_cluster_composition.csv"), row.names = FALSE)
cat("\n-- tc_cluster composition --\n")
tc_comp %>% select(-clusters) %>% as.data.frame() %>% print(row.names = FALSE)

# -- 3. where each D1 cluster came from --
# For every cluster each D1 genotype uses, the cluster's history in the WHOLE
# dataset: how many sequences and tc_genos carry it, when and where it was first
# sampled, and -- excluding the focal genotype itself -- which other genotypes
# held it and how long before. This is what identifies a donor for a
# single-segment difference: a cluster with a long prior history in another
# genotype was acquired, not invented.
#
# BLIND SPOT: this panel is H5 only, so a segment acquired from a non-H5 IAV has
# no candidate origin here at all. Absence of a prior holder means "not in this
# panel", NOT "novel segment".

fmt_counts <- function(x, k = 4) {
  tb <- sort(table(x), decreasing = TRUE)
  tb <- tb[seq_len(min(k, length(tb)))]
  paste0(names(tb), "(", as.integer(tb), ")", collapse = ";")
}

first_by_continent <- function(df) {
  df %>%
    group_by(continent_corrected) %>%
    summarise(f = min(coll_lb), .groups = "drop") %>%
    arrange(f) %>%
    transmute(s = paste0(continent_corrected, ":", f)) %>%
    pull(s) %>% paste(collapse = ";")
}

cluster_origin <- function(geno_val) {
  focal <- d1 %>% filter(genoflu_geno == geno_val)
  focal_first <- min(focal$coll_lb)
  map_dfr(seg.order, function(s) {
    col  <- paste0("tc_cluster", s)
    used <- focal %>% count(cluster = .data[[col]], name = "n_focal") %>% arrange(desc(n_focal))
    map_dfr(seq_len(nrow(used)), function(i) {
      cl    <- used$cluster[i]
      hold  <- meta %>% filter(.data[[col]] == cl)
      other <- hold %>% filter(genoflu_geno != geno_val)
      tibble(geno = geno_val, segment = s, cluster = cl,
             n_focal = used$n_focal[i], focal_first_date = focal_first,
             n_seqs_total     = nrow(hold),
             n_tc_geno_total  = n_distinct(hold$tc_geno),
             first_date_total = min(hold$coll_lb),
             first_by_continent = first_by_continent(hold),
             n_seqs_other     = nrow(other),
             first_date_other = if (nrow(other)) min(other$coll_lb) else as.Date(NA),
             prior_days       = if (nrow(other)) as.numeric(focal_first - min(other$coll_lb))
                                else NA_real_,
             other_genos      = if (nrow(other)) fmt_counts(other$genoflu_geno) else NA_character_,
             other_continents = if (nrow(other)) fmt_counts(other$continent_corrected) else NA_character_)
    })
  })
}

origins <- map_dfr(focal.genos, cluster_origin) %>%
  mutate(segment = factor(segment, levels = seg.order)) %>%
  arrange(geno, segment, desc(n_focal))

write.csv(origins, paste0(f.out, "d1_cluster_origins.csv"), row.names = FALSE)

# the segments that actually differ between the main pair are where the origin
# question matters, so print those rows in full
diff_segs <- gf_feature %>% filter(differs) %>% pull(segment) %>% as.character()
cat("\n-- cluster origins on the segment(s) that differ:",
    paste(diff_segs, collapse = ", "), "--\n")
origins %>%
  filter(segment %in% diff_segs) %>%
  select(geno, segment, cluster, n_focal, focal_first_date, n_seqs_total,
         n_tc_geno_total, first_date_total, n_seqs_other, first_date_other,
         prior_days, other_genos) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nfirst detection by continent, differing segment(s):\n")
origins %>% filter(segment %in% diff_segs) %>%
  select(geno, segment, cluster, first_by_continent) %>%
  as.data.frame() %>% print(row.names = FALSE)

# -- 4. within- and between-genotype patristic distance per segment --
set.seed(dist.seed)

# cap a genotype's tip indices so the distance matrix stays a sane size
cap_idx <- function(idx) if (length(idx) > max.tips) sample(idx, max.tips) else idx

dist_rows <- map_dfr(seg.order, function(s) {
  sd  <- load_seg(s, meta, "genoflu_geno")
  idx <- map(set_names(focal.genos), function(g) {
    cap_idx(sd$tips$idx[which(sd$tips$genoflu_geno == g)])
  })

  # within each genotype
  within <- map_dfr(focal.genos, function(g) {
    if (length(idx[[g]]) < min.n) return(NULL)
    dv <- pw_dists(sd$tree, idx[[g]])
    tibble(comparison = "within", geno_a = g, geno_b = g, n_a = length(idx[[g]]),
           n_b = length(idx[[g]]), n_pairs = length(dv),
           d_min = min(dv), d_mean = mean(dv), d_median = median(dv), d_max = max(dv))
  })

  # between the main pair: one matrix over both tip sets, then the off-diagonal block
  between <- NULL
  if (length(idx[[pair.a]]) >= 1 && length(idx[[pair.b]]) >= 1) {
    ia <- idx[[pair.a]]; ib <- idx[[pair.b]]
    m  <- castor::get_all_pairwise_distances(sd$tree, only_clades = c(ia, ib),
                                             check_input = FALSE)
    blk <- as.vector(m[seq_along(ia), length(ia) + seq_along(ib)])
    between <- tibble(comparison = "between", geno_a = pair.a, geno_b = pair.b,
                      n_a = length(ia), n_b = length(ib), n_pairs = length(blk),
                      d_min = min(blk), d_mean = mean(blk),
                      d_median = median(blk), d_max = max(blk))
  }
  bind_rows(within, between) %>% mutate(segment = s, .before = 1)
}) %>%
  mutate(segment = factor(segment, levels = seg.order))

write.csv(dist_rows, paste0(f.out, "d1_patristic_by_segment.csv"), row.names = FALSE)

# separation = between-pair median / the larger of the two within-genotype medians
separation <- dist_rows %>%
  select(segment, comparison, geno_a, d_median) %>%
  pivot_wider(names_from = c(comparison, geno_a), values_from = d_median) %>%
  mutate(within_max = pmax(.data[[paste0("within_", pair.a)]],
                           .data[[paste0("within_", pair.b)]], na.rm = TRUE),
         between    = .data[[paste0("between_", pair.a)]],
         separation = round(between / within_max, 2)) %>%
  select(segment, !!paste0("within_", pair.a), !!paste0("within_", pair.b),
         between, separation) %>%
  arrange(desc(separation))

write.csv(separation, paste0(f.out, "d1_segment_separation.csv"), row.names = FALSE)
cat("\n-- per-segment separation,", pair.a, "vs", pair.b, "(higher = more distinct) --\n")
as.data.frame(separation) %>% print(row.names = FALSE)

# within vs between median distance per segment
p <- dist_rows %>%
  mutate(label = ifelse(comparison == "within", paste0("within ", geno_a),
                        paste0("between ", geno_a, " / ", geno_b))) %>%
  ggplot(aes(x = segment, y = d_median, fill = label)) +
  geom_col(position = position_dodge(preserve = "single"), colour = "grey30",
           linewidth = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  plot.theme +
  labs(x = "segment", y = "median pairwise patristic distance (subs/site)", fill = NULL)
ggsave(paste0(f.out, "d1_patristic_within_between_by_segment.tiff"),
       p, width = 8, height = 5, units = "in", dpi = 300)

cat("\nwrote", f.out, "d1_genoflu_feature_by_segment.csv + d1_tc_cluster_composition.csv +",
    "d1_cluster_origins.csv + d1_patristic_by_segment.csv + d1_segment_separation.csv +",
    "d1_patristic_within_between_by_segment.tiff\n")
