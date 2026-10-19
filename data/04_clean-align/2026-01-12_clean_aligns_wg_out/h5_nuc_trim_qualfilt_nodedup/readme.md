HA nuc seqs aligned, gap trimmed, and quality filtered, but not de-duplicated 
removed seq if >10% ambiguous (low qual) or <50% coverage after trimming gaps

Main files (and original locations):
v2_h5nx_genomic_ml/data/04_clean-align/2026-01-12_clean_aligns_wg_out/cleaned_fasta/trim_filt_nuc_HA_aligned.fasta
v2_h5nx_genomic_ml/data/04_clean-align/2026-01-12_clean_aligns_wg_out/cleaned_meta/merged_cleaned_filt_meta_HA.csv

Host taxonomic labels:
v2_h5nx_genomic_ml/data/06_host/2026-01-27_out/complete_host_taxa.csv

With the 'taxa2' column that is the multi-level taxa being used in ML:
v2_h5nx_genomic_ml/data/06_host/multi-taxa/2026-03-04_out/all-taxa.csv

Note the host taxa files are connected to the de-duplicated meta, so if you are joining by sample, there may be some that are incomplete in the non-de-dup data. Should be able to fill those in based on the common host in the original metadata though