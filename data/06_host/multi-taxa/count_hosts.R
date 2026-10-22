library(tidyverse)                    
f.out<-"2026-06-19_out/"
samp.look.all<-read.csv(paste0(f.out, "sample_taxa_lookup_all.csv"))
host.counts<-samp.look.all %>% count(common, species, taxa_level, taxa1p,
                                     name="n_samples") %>% arrange(-n_samples)
write.csv(host.counts, paste0(f.out, "host_counts_ungrouped.csv"),
          row.names=FALSE)

