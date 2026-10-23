# geno_reassort

Within- and between-genotype phylogenetic divergence per segment, plus small
targeted queries into the composition of individual genotypes of interest.

## Main analysis

`geno_reassort.Rmd` — within-genotype pairwise patristic distance per segment for
every `tc_geno` / `genoflu_geno` / `genin_geno`, on the rooted IQ-TREE segment
trees. Distances are RAW subs/site, deliberately not normalised; see the Rmd for
the two normalisers that were tried and rejected.

## Query scripts

Each is standalone, sources `_common.R`, and writes to `<date>_geno_queries_out/`.
Run from this folder with `Rscript <script>`. Scripts 01-03 read CSVs and
metadata only and finish in seconds; 04 reads all eight segment trees and
takes a few minutes.

- `01_control_genos_spread.R` — positive and negative controls for the cluster-ID
  span vs number-of-clusters screen; B3.5 is scattered on 5 segments (max ratio
  10.7), A6 scores 1.0 on all 8
- `02_reassort_cluster_priority.R` — for the genotype × segment rows the screen
  flags, whether each member cluster was already circulating in another genotype
  first, with origin genotype, date lag and origin continent/subtype/clade
- `03_h5n2_cohorts.R` — composition of the two H5N2 cohorts of interest, Mexican
  lineage and the North American 2014-15 epizootic
- `04_d1_genotype_compare.R` — which segments separate the D1 GenoFlu genotypes,
  by GenoFlu feature, tc_cluster composition, within/between patristic distance
  and the prior history of every cluster each genotype uses

## Findings so far

- The BC 2014 H5N2 reassortant (alias TC1583; Pasick, Berhane et al. 2015,
  doi:10.1038/srep09484) works as a documented positive control: NP and NA trace
  to North American `Am_nonGsGD` genotypes from 2009-2010 and PB2, M and NS to
  Asian H5N8/H5N2 clade 2.3.4.4b, matching the paper's split of PB1/NP/NA from
  North American waterfowl against PB2/PA/HA/M/NS from Eurasian H5N8; 7 of 8
  segments recover the published donor, PB1 being first seen in TC1583 itself
- A6 is clonal in cluster terms — 315 sequences, one `tc_geno` (TC3173) — so it
  is the natural negative control for the span/ratio screen, and any within-A6
  structure has to be measured from the trees rather than from cluster labels
- D1.1 and D1.3 differ on NA alone (`am4N1` vs `ea3`, separation 48.4×), and
  D1.2 differs from D1.1 on PA alone (`am5` vs `am4`); PA and NP are identical
  between D1.1 and D1.3, so the PA/NP hypothesis does not describe these genotypes
- the shared element between D1.1 and D1.3 is not one predominant `tc_geno` but a
  six-segment identical constellation (PB1-631, PA-156, HA-445, NP-111, M-332,
  NS-592) plus PB2 clusters from the same contiguous block, D1.3's three being a
  subset of D1.1's 21; PB2 is the only segment that varies within either
  genotype and its clusters are contiguous (ratio 1.0), so that is fine-scale
  subdivision of one PB2 lineage, not reassortment
- D1.3's NA donor is identifiable and well sampled: NA-662 appears from
  2021-10-28 in Europe, reaches North America 2022-01-01 in A3-genotype wild
  birds, and is carried by 963 of 971 A3 sequences — roughly two years before
  D1.3. D1.1's NA-334 has no prior holder at all, first appearing with D1.1
  itself, so which of the two carries the derived NA is not settled by these
  dates (both genotypes' earliest records are year-precision 2024)

## Gotchas

- The `segment` column of the subtree_cluster_matrix long CSVs literally holds
  `NA` for neuraminidase; `read_seg_csv()` rebuilds it from the `cluster` prefix
  rather than reading the whole file with `na.strings = ""`, which would leave
  genuinely-missing values as the string `"NA"` and break the date columns
- subtree_cluster_matrix re-runs write only the chunks that were re-run, so its
  newest dated directory does not hold every CSV; `newest_in()` takes the newest
  directory that actually contains the requested file
