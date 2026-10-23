library(tidyverse)

meta.in <- "../../../data/08_clean-share/meta/metadata.csv"
meta <- read.csv(meta.in)

id.cols <- c("sample_unif", "accession_id", "isolate_name", "collection_date",
             "country", "common_host", "tc_geno", "treecluster_genotype_alias",
             "genoflu_geno", "genin_geno")

today <- Sys.Date()
f.out <- paste0(today, "_missed_assignment_out/")
if (!dir.exists(f.out)) {dir.create(f.out)}

find_missed <- function(data, geno.col, unassigned.label) {
  tc_summary <- data %>%
    filter(!is.na(tc_geno), tc_geno != "") %>%
    group_by(tc_geno) %>%
    summarise(
      n_seqs = n(),
      n_unassigned = sum(.data[[geno.col]] == unassigned.label),
      genos_assigned = paste(sort(unique(.data[[geno.col]][.data[[geno.col]] != unassigned.label])), collapse = ";"),
      n_genos_assigned = n_distinct(.data[[geno.col]][.data[[geno.col]] != unassigned.label]),
      .groups = "drop"
    ) %>%
    filter(n_unassigned > 0, n_genos_assigned >= 1) %>%
    arrange(desc(n_unassigned), desc(n_seqs))

  seq_detail <- data %>%
    filter(tc_geno %in% tc_summary$tc_geno) %>%
    select(all_of(id.cols)) %>%
    left_join(tc_summary %>% select(tc_geno, n_seqs, n_unassigned, genos_assigned, n_genos_assigned),
               by = "tc_geno") %>%
    arrange(desc(n_unassigned), desc(n_seqs), tc_geno)

  list(tc_summary = tc_summary, seq_detail = seq_detail)
}

genoflu_missed <- find_missed(meta, "genoflu_geno", "Not assigned")
genin_missed   <- find_missed(meta, "genin_geno", "[unassigned]")

cat("tc_geno groups with >=1 GenoFlu-unassigned seq alongside a GenoFlu-assigned genotype:",
    nrow(genoflu_missed$tc_summary), "\n")
cat("tc_geno groups with >=1 Genin2-unassigned seq alongside a Genin2-assigned genotype:",
    nrow(genin_missed$tc_summary), "\n")

write.csv(genoflu_missed$tc_summary, paste0(f.out, "genoflu_missed_tc_summary.csv"), row.names = FALSE)
write.csv(genoflu_missed$seq_detail, paste0(f.out, "genoflu_missed_seq_detail.csv"), row.names = FALSE)
write.csv(genin_missed$tc_summary, paste0(f.out, "genin_missed_tc_summary.csv"), row.names = FALSE)
write.csv(genin_missed$seq_detail, paste0(f.out, "genin_missed_seq_detail.csv"), row.names = FALSE)

cat("\n--- GenoFlu example with most unassigned sequences ---\n")
top.tc.genoflu <- genoflu_missed$tc_summary$tc_geno[1]
print(genoflu_missed$seq_detail %>% filter(tc_geno == top.tc.genoflu))

cat("\n--- Genin2 example with most unassigned sequences ---\n")
top.tc.genin <- genin_missed$tc_summary$tc_geno[1]
print(genin_missed$seq_detail %>% filter(tc_geno == top.tc.genin))

tc_alias_lookup <- meta %>%
  distinct(tc_geno, treecluster_genotype_alias)

export_missed <- function(tc_summary, scheme) {
  tc_summary %>%
    left_join(tc_alias_lookup, by = "tc_geno") %>%
    transmute(
      scheme = scheme,
      tc_geno,
      tc_cluster = treecluster_genotype_alias,
      geno = genos_assigned,
      N_assigned = n_seqs - n_unassigned,
      N_unassigned = n_unassigned
    )
}

missed_export <- bind_rows(
  export_missed(genoflu_missed$tc_summary, "genoflu"),
  export_missed(genin_missed$tc_summary, "genin")
) %>%
  arrange(desc(N_unassigned), desc(N_assigned))

write.csv(missed_export, paste0(f.out, "missed_assignment_examples.csv"), row.names = FALSE)
cat("\nwrote", paste0(f.out, "missed_assignment_examples.csv"), "-", nrow(missed_export), "rows\n")
