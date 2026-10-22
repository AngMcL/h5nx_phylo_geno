# data
Due to limitations with data sharing, we have excluded the data files from this repo. Authorized GISAID users may contact us directly for files. Accession lists and episet doi provided.

## 01a_ncbi	
* downloading and cleaning unaligned ncbi data

## 01b_gisaid	
* downloading and cleaning unaligned gisaid data

## 02_merge	
* merge all data sources together, remove redundant IDs

## 03_align	
* align segments and proteins with mafft

## 04_clean-align	
* exclude aligned seqs with >5%/10% N, <90%? length, identical seqs (within species)

## 05_annotate
* annotate the taxonomy, mutations, genetic features, ...

## 06_host
* harmonize host naming and look up full taxa set

## 07_summarize_clean
* summary plots of the final clean data analyzed
