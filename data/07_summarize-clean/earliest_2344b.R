# Earliest detection (year / place / subtype) of H5 clade 2.3.4.4b in the final clean dataset,
# using both the phylogenetically filled-in clade (as plotted in SummarizeCleanData_v4.Rmd) and
# the pre-fill clade label from the HA metadata.
# Inputs match SummarizeCleanData_v4.Rmd / dataset_milestones.md. One HA record per sample.

library(dplyr)

f.in.meta <- "../08_all-together/nuc_seq/trim_nodup_nuc_HA_rooted_meta.csv"
f.in.fill <- "../../phylo/03_clusters/subtype_clades/2026-07-01_out/sample_filledin_N-type_HA-clade.csv"

meta <- read.csv(f.in.meta) %>%
  select(sample, isolate_id, accession, collection_year, collection_date, country, location,
         host, subtype.raw = subtype, clade.raw = clade)
stopifnot(!anyDuplicated(meta$sample))

# filled subtype = H5 + filled N type (same rule as SummarizeCleanData_v4.Rmd)
fill <- read.csv(f.in.fill) %>%
  mutate(subtype.fill = if_else(!is.na(N_type), paste0("H5", N_type), subtype)) %>%
  select(sample, clade.fill = clade, subtype.fill)

meta <- left_join(meta, fill, by = "sample")
stopifnot(!any(is.na(meta$clade.fill)))

# GISAID vs NCBI inferred from sample ID prefix
meta$source <- if_else(grepl("^EPI_ISL_", meta$sample), "GISAID", "NCBI")

# earliest 2.3.4.4b by collection_year (collection_date is NA for partial dates). Returns all
# samples from the first n.yr years with any 2.3.4.4b, so a lone early outlier (e.g. a mistyped
# date) doesn't hide the next real detection year.
earliest <- function(df, clade.col, n.yr = 2) {
  df <- df %>% filter(.data[[clade.col]] == "2.3.4.4b")
  yrs <- head(sort(unique(df$collection_year)), n.yr)
  df %>%
    filter(collection_year %in% yrs) %>%
    arrange(collection_year, collection_date) %>%
    select(sample, source, isolate_id, accession, collection_year, collection_date, country,
           location, host, subtype.raw, subtype.fill, clade.raw, clade.fill)
}

# filled-in clade
e.fill <- earliest(meta, "clade.fill")
cat("Filled-in clade: 2.3.4.4b samples in the first 2 detection years\n")
print(e.fill)

# pre-fill (raw metadata) clade
e.raw <- earliest(meta, "clade.raw")
cat("\nPre-fill clade: 2.3.4.4b samples in the first 2 detection years\n")
print(e.raw)

# per-year 2.3.4.4b counts for the first few years, filled vs pre-fill
cat("\nFirst years with 2.3.4.4b (filled vs pre-fill):\n")
print(meta %>%
        group_by(collection_year) %>%
        summarise(n.fill = sum(clade.fill == "2.3.4.4b"),
                  n.raw  = sum(clade.raw == "2.3.4.4b", na.rm = T)) %>%
        filter(n.fill > 0 | n.raw > 0) %>%
        head(8))

# FLAG: the single 2002 sample (EPI_ISL_19559114, H5N1 chicken, Snohomish County, WA,
# 2002-06-06) is almost certainly a date error -- 2.3.4.4b did not exist in 2002, and the next
# detection is 2010. Treat 2010 as the first credible 2.3.4.4b year in the dataset.
#
# NOTE: earlier 2.3.4.4b exists outside this dataset. EPI_ISL_180649 (China, 2009) is labelled
# 2.3.4.4b in GISAID but is not in our download (whole-genome isolates only), and other partial-
# genome sequences with this label likely exist. 2010 is the earliest in THIS dataset only.
#
# Citations for the earliest documented 2.3.4.4b viruses (original isolate papers and the clade
# nomenclature, incl. Lee et al. 2017 J Vet Sci 18(S1):269-280 for the "group B" division that
# became 2.3.4.4b) are collected in dataset_milestones.md, section "Earliest documented
# 2.3.4.4b: citations". Note GISAID's Publication field for these isolates lists only later
# papers that reused the sequences, not their original descriptions.
