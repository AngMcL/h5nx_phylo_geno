# Shared setup for the small genotype-composition query scripts in this folder.
# Sourced by each query script; holds paths, theme, segment order and the
# tree/metadata loaders so the query scripts stay short.

library(tidyverse)
library(ape)
library(castor)

plot.theme <- theme(panel.grid.minor = element_blank(),
                    panel.background = element_rect("grey95"),
                    axis.line = element_line(colour = "black"),
                    text = element_text(size = 10, face = "bold"),
                    legend.text = element_text(size = 8),
                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

seg.order <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

seg.colors <- setNames(
  c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB", "#000000"),
  seg.order)

meta.in   <- "../../../data/08_clean-share/meta/metadata.csv"
tree.dir  <- "../../02_iqtree/iq_gtrfr10_rooted"
tree.tmpl <- "trim_nodup_nuc_%s_rooted_remiqtree_rooted.treefile"
resid.dir <- "../../01_fasttree/ft_rooted_outliers/2026-01-14_rem-outliers/rem_out"
resid.tmpl <- "trim_nodup_nuc_%s_rooted_resid.csv"
scm.dir   <- "../subtree_cluster_matrix"

drop.vals <- c("Not assigned", "[unassigned]", "")

today <- Sys.Date()
f.out <- paste0(today, "_geno_queries_out/")
if (!dir.exists(f.out)) dir.create(f.out)

# Newest dated output directory under `root` that actually CONTAINS `file`.
# Needed because subtree_cluster_matrix re-runs write only the chunks that were
# re-run, so the newest directory does not hold every CSV (e.g. the sharing
# tables are in 2026-09-23 while the summary/spread tables are in 2026-09-27).
# Uses sort(), NOT gtools::mixedsort(), which reverse-sorts ISO date folders
# because it reads the dashes as minus signs.
newest_in <- function(root, pattern, file) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(pattern, basename(dirs))]
  dirs <- dirs[file.exists(file.path(dirs, file))]
  if (!length(dirs)) {
    stop("no directory under ", root, " matching '", pattern, "' contains ", file)
  }
  p <- file.path(sort(dirs)[length(dirs)], file)
  cat("input:", p, "\n")
  p
}

# collection_date is mixed precision (YYYY, YYYY-MM, YYYY-MM-DD) and is never NA;
# pad to the earliest possible day so partial dates still sort for "earliest".
parse_date_lb <- function(x) {
  x <- ifelse(nchar(x) == 4, paste0(x, "-01-01"),
       ifelse(nchar(x) == 7, paste0(x, "-01"), x))
  as.Date(x)
}

read_meta <- function() {
  read.csv(meta.in) %>% mutate(coll_lb = parse_date_lb(collection_date))
}

# The `segment` column of the subtree_cluster_matrix long CSVs literally holds
# "NA" for neuraminidase, which read.csv turns into a missing value. Reading the
# whole file with na.strings = "" would keep it, but it would ALSO keep every
# genuinely-missing "NA" in the other columns as the string "NA" -- which breaks
# the date columns (as.Date("NA") errors). So read with the usual na.strings and
# rebuild `segment` from the `cluster` prefix instead: cluster IDs are
# segment-prefixed and globally unique ("NA-39", "NA-s", never a bare "NA").
read_seg_csv <- function(path) {
  df <- read.csv(path)
  if (all(c("segment", "cluster") %in% names(df))) {
    df$segment <- sub("-.*$", "", df$cluster)
  } else {
    stop("cannot rebuild the segment column of ", path, ": no `cluster` column")
  }
  df
}

# One segment's tree plus a tip -> sample -> metadata lookup, in tip order.
# `cols` are metadata columns to carry onto the tips.
load_seg <- function(s, meta, cols) {
  tr <- read.tree(file.path(tree.dir, sprintf(tree.tmpl, s)))
  rs <- read.csv(file.path(resid.dir, sprintf(resid.tmpl, s)), check.names = FALSE)
  tips <- tibble(tip = tr$tip.label, idx = seq_along(tr$tip.label)) %>%
    left_join(rs %>% select(tip = seq_new, sample), by = "tip") %>%
    left_join(meta %>% select(all_of(c("sample", cols))), by = "sample")
  cat(sprintf("%-3s: %d tips | unmatched %d\n", s, Ntip(tr), sum(is.na(tips$sample))))
  list(tree = tr, tips = tips)
}

# Lower-triangle pairwise patristic distances among a set of tip indices.
# RAW branch-length units (subs/site) -- see geno_reassort.Rmd for why these are
# deliberately not normalised.
pw_dists <- function(tree, idx) {
  d <- castor::get_all_pairwise_distances(tree, only_clades = idx, check_input = FALSE)
  d[lower.tri(d)]
}

# Numeric tc_cluster IDs for one segment (singleton "-s" buckets dropped).
# Anchored to the END of the string: an unanchored "\\d+" would grab the digit
# inside the segment name itself (PB2 -> 2, PB1 -> 1).
cluster_nums <- function(ids) {
  sort(unique(as.numeric(str_extract(unique(ids), "\\d+$"))))
}
