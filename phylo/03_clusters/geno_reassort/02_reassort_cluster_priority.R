# For the putative reassortants flagged by the span/n_clusters screen: were their
# cluster IDs found in OTHER genotypes first?
#
# A scattered genotype x segment (ratio >= flag.ratio) says the genotype's
# sequences sit in unrelated parts of that segment's tree. That is consistent
# with reassortment, but only if the segment was acquired -- i.e. if the cluster
# it sits in was already circulating in a different genotype BEFORE this one
# appeared. This script answers that per cluster, using the cluster-sharing table
# (tc_geno_segment_sharing.csv), which already records, for every tc_geno x
# segment, whether that exact cluster appears in any other tc_geno and which
# sharer is earliest (origin) vs latest (recipient).
#
# Reading the output:
#   focal_role == "recipient" / "intermediate"  the cluster existed in another
#       genotype before this one -> acquired, consistent with reassortment
#   focal_role == "origin"                      this genotype had it first
#   focal_role == "unique"                      no other genotype carries it
#   lag_days   days between the origin's first sequence and the focal's first
#   origin_diff_geno  TRUE when the origin tc_geno carries a DIFFERENT
#       genoflu/genin label than the focal genotype, i.e. the cluster was first
#       seen outside this genotype rather than in a sibling tc_geno of it
#
# Singleton "-s" buckets are never counted as shared (they pool unrelated
# singletons), so they can neither have nor be an origin.
#
# POSITIVE CONTROL -- the BC 2014 H5N2 reassortant (Pasick, Berhane et al. 2015,
# Sci Rep 5:9484, doi:10.1038/srep09484). That outbreak virus is a documented
# per-segment reassortant: PB2, PA, HA, M and NS derive from Eurasian HPAI H5N8
# and PB1, NP and NA from North American lineage waterfowl viruses. In this
# dataset it is alias TC1583 (clade 2.3.4.4c; 10 British Columbia sequences from
# Dec 2014, and TC1583 covers 285 of the 311 North American H5N2 2014-15
# sequences). Because the parentage is known per segment, it is a real test of
# the origin-tracing logic rather than of the span/ratio screen: the three North
# American segments should trace to earlier North American genotypes and the five
# Eurasian ones should not share an origin with them.

source("_common.R")

flag.ratio    <- 2         # screen threshold defining "putative reassortant"
control.alias <- "TC1583"  # BC 2014 H5N2, documented per-segment reassortant

# documented donor lineage per segment for the control (Pasick et al. 2015)
control.expected <- tibble(
  segment = seg.order,
  expected_donor = c("Eurasian H5N8", "North American waterfowl", "Eurasian H5N8",
                     "Eurasian H5N8", "North American waterfowl",
                     "North American waterfowl", "Eurasian H5N8", "Eurasian H5N8"))

spread.in  <- newest_in(scm.dir, "_subtree_cluster_matrix_out$", "geno_tc_cluster_spread.csv")
sharing.in <- newest_in(scm.dir, "_subtree_cluster_matrix_out$", "tc_geno_segment_sharing.csv")

meta    <- read_meta()
spread  <- read.csv(spread.in, na.strings = "")
sharing <- read_seg_csv(sharing.in) %>%
  mutate(earliest_date  = as.Date(earliest_date),
         origin_date    = as.Date(origin_date),
         recipient_date = as.Date(recipient_date))

geno.cols <- c(genoflu = "genoflu_geno", genin = "genin_geno")

# tc_geno alias -> the genoflu / genin labels assigned to its sequences
alias_labels <- map_dfr(names(geno.cols), function(sch) {
  col <- geno.cols[[sch]]
  meta %>%
    filter(!is.na(.data[[col]]), !.data[[col]] %in% drop.vals) %>%
    count(alias = treecluster_genotype_alias, geno = .data[[col]], name = "n") %>%
    mutate(scheme = sch)
})

# alias -> where/what that tc_geno is, used to judge an origin: for the BC 2014
# control the discriminating evidence is not WHETHER a segment had an earlier
# origin (most do) but WHOSE genotype that origin is -- the three North American
# segments should trace to North American genotypes and the five Eurasian ones
# should not.
alias_context <- meta %>%
  group_by(alias = treecluster_genotype_alias) %>%
  summarise(n_seqs    = n(),
            continent = paste0(names(sort(table(continent_corrected), decreasing = TRUE)),
                               collapse = ";"),
            subtype   = paste0(names(sort(table(subtype_filled), decreasing = TRUE)),
                               collapse = ";"),
            clade     = paste0(names(sort(table(clade), decreasing = TRUE)),
                               collapse = ";"),
            first_date = min(coll_lb),
            .groups = "drop")

# alias -> compact "Label(n);Label(n)" string, per scheme
# `labels` carries per-genotype counts for reading; `labels_plain` is the bare
# genotype set, needed wherever the values are deduplicated -- two origins of the
# same genotype with different counts ("EA-2021-AB(2)" vs "EA-2021-AB(94)") would
# otherwise list as two distinct donors.
alias_label_str <- alias_labels %>%
  group_by(scheme, alias) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  summarise(labels       = paste0(geno, "(", n, ")", collapse = ";"),
            labels_plain = paste(sort(unique(geno)), collapse = ";"),
            .groups = "drop")

# genotype -> its member tc_geno aliases, per scheme
geno_members <- alias_labels %>% select(scheme, geno, alias) %>% distinct()

flagged <- spread %>%
  filter(ratio >= flag.ratio) %>%
  select(scheme, geno, n_seqs, segment, n_numeric_clusters, span, ratio)

cat("\nflagged genotype x segment rows (ratio >=", flag.ratio, "):", nrow(flagged), "\n")

# one row per flagged genotype x segment x member cluster
priority <- flagged %>%
  left_join(geno_members, by = c("scheme", "geno"), relationship = "many-to-many") %>%
  left_join(sharing %>% select(alias, tc_geno, segment, cluster, is_singleton, shared,
                               n_other_genos, earliest_date, origin_alias, origin_date,
                               recipient_alias, recipient_date, focal_role),
            by = c("alias", "segment")) %>%
  # does the origin tc_geno carry a different label than the focal genotype?
  left_join(alias_label_str %>% rename(origin_alias = alias, origin_labels = labels,
                                       origin_labels_plain = labels_plain),
            by = c("scheme", "origin_alias")) %>%
  mutate(lag_days = as.numeric(earliest_date - origin_date),
         origin_diff_geno = case_when(
           is.na(origin_alias) ~ NA,
           is.na(origin_labels) ~ TRUE,
           TRUE ~ !str_detect(origin_labels, fixed(paste0(geno, "(")))),
         acquired = focal_role %in% c("recipient", "intermediate"),
         # an origin tc_geno that carries no label in THIS scheme (e.g. GenoFlu
         # "Not assigned") has no origin_labels, but it is still a real prior
         # holder and must stay nameable, otherwise it counts towards
         # n_origin_outside while contributing nothing to prior_genos. Fall back
         # to its alias so the two columns agree.
         origin_name = coalesce(origin_labels_plain,
                                if_else(is.na(origin_alias), NA_character_,
                                        paste0(origin_alias, "[no ", scheme, " label]")))) %>%
  select(scheme, geno, n_seqs, segment, ratio, alias, tc_geno, cluster, is_singleton,
         shared, n_other_genos, focal_role, acquired, earliest_date,
         origin_alias, origin_date, lag_days, origin_diff_geno, origin_labels,
         origin_name) %>%
  arrange(scheme, desc(ratio), geno, segment, desc(lag_days))

write.csv(priority, paste0(f.out, "reassort_cluster_priority.csv"), row.names = FALSE)

# rollup: per flagged genotype x segment, how many of its clusters were already
# circulating elsewhere, and the longest head start any origin had
rollup <- priority %>%
  group_by(scheme, geno, n_seqs, segment, ratio) %>%
  summarise(n_tc_geno        = n_distinct(alias),
            n_clusters       = n_distinct(cluster[!is_singleton]),
            n_shared         = sum(shared, na.rm = TRUE),
            n_acquired       = sum(acquired, na.rm = TRUE),
            n_origin         = sum(focal_role == "origin", na.rm = TRUE),
            n_unique         = sum(focal_role == "unique", na.rm = TRUE),
            n_origin_outside = sum(origin_diff_geno, na.rm = TRUE),
            max_lag_days     = suppressWarnings(max(lag_days, na.rm = TRUE)),
            # the distinct genotypes that held one of this genotype's clusters
            # BEFORE it did. Note this names only the single earliest holder per
            # cluster (not every prior holder), matches clusters by exact ID, and
            # orders "earliest" on dates padded by parse_date_lb -- so treat it as
            # the leading candidate donor set, not an exhaustive one.
            prior_genos      = paste(sort(unique(na.omit(
                                 origin_name[origin_diff_geno]))), collapse = " | "),
            .groups = "drop") %>%
  mutate(max_lag_days = ifelse(is.finite(max_lag_days), max_lag_days, NA_real_)) %>%
  arrange(scheme, desc(ratio))

write.csv(rollup, paste0(f.out, "reassort_cluster_priority_rollup.csv"), row.names = FALSE)

cat("\n-- rollup, top flagged rows --\n")
rollup %>% slice_head(n = 20) %>%
  select(scheme, geno, segment, ratio, n_clusters, n_shared, n_acquired,
         n_origin_outside, max_lag_days) %>%
  as.data.frame() %>% print(row.names = FALSE)

# -- positive control: the documented BC 2014 H5N2 reassortant --
# CAVEAT on lag_days here: alias TC1583 covers 285 of the 311 North American
# H5N2 2014-15 sequences, not just the 10 British Columbia ones, and its earliest
# record carries year-only precision (2014), which parse_date_lb() pads to
# 2014-01-01. So the focal date is approximate and the lags are indicative of
# direction, not exact intervals. Read origin_continent / origin_clade, not
# lag_days, as the discriminating evidence.
control <- sharing %>%
  filter(alias == control.alias) %>%
  left_join(control.expected, by = "segment") %>%
  left_join(alias_label_str %>% filter(scheme == "genoflu") %>%
              select(origin_alias = alias, origin_genoflu = labels),
            by = "origin_alias") %>%
  left_join(alias_context %>%
              select(origin_alias = alias, origin_continent = continent,
                     origin_subtype = subtype, origin_clade = clade),
            by = "origin_alias") %>%
  mutate(lag_days = as.numeric(earliest_date - origin_date),
         acquired = focal_role %in% c("recipient", "intermediate"),
         segment  = factor(segment, levels = seg.order)) %>%
  arrange(segment) %>%
  select(alias, tc_geno, segment, expected_donor, cluster, is_singleton, shared,
         n_other_genos, focal_role, acquired, earliest_date, origin_alias,
         origin_date, lag_days, origin_continent, origin_subtype, origin_clade,
         origin_genoflu, other_aliases)

write.csv(control, paste0(f.out, "control_BC2014_H5N2_cluster_priority.csv"), row.names = FALSE)

cat("\n-- positive control", control.alias, "(BC 2014 H5N2, Pasick et al. 2015) --\n")
control %>%
  select(segment, expected_donor, cluster, focal_role, origin_alias, origin_date,
         origin_continent, origin_subtype, origin_clade) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nwrote", f.out, "reassort_cluster_priority.csv + _rollup.csv +",
    "control_BC2014_H5N2_cluster_priority.csv\n")
