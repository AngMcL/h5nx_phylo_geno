# 07_summarize-clean

Summarizes the final cleaned H5Nx dataset with figures for nucleotide sequence quality, host composition, geographic distribution, and subtype breakdown.

## Purpose

Reads the post-cleaned, post-outlier-removed nucleotide sequences and metadata, applies host taxonomy corrections and subtype fill-ins, then produces a suite of exploratory summary plots.

## Scripts

| Script | Status | Description |
|--------|--------|-------------|
| `SummarizeCleanData_v2.Rmd` | archived | Original version; reads directly from the phylo pipeline output (`phylo/01_fasttree/...rem_out/`). Contains bugs listed below. |
| `SummarizeCleanData_v3.Rmd` | **current** | Reads from the canonical `08_all-together/` directory; fixes host duck→mallard error, South Georgia country rename bug, and segment assignment; switches from FastTree to IQTREE trees. |

## Inputs

All inputs are from `../08_all-together/`:

| Path | Contents |
|------|----------|
| `nuc_seq/trim_nodup_nuc_*_rooted_meta.csv` | Per-segment metadata (8 files); post-temporal-outlier removal |
| `nuc_seq/trim_nodup_nuc_*_rooted_rem.fasta` | Trimmed, deduplicated nucleotide alignments (8 files) |
| `phylo_trees/trim_nodup_nuc_*_rooted_remiqtree_rooted.treefile` | IQTREE2 GTR+F+R10 rooted trees (8 files) |
| `annot/sample_taxa_lookup_all.csv` | Per-sample host taxonomy with corrected common names (duck/mallard fix) and multi-taxa groupings |
| `annot/sample_filledin_N-type_HA-clade.csv` | Phylogenetically filled-in N-type, subtype, and H5 clade assignments |

## Outputs

Date-stamped folder (`YYYY-MM-DD_out/`) containing:

| File pattern | Description |
|---|---|
| `n.clean.seq.byseg*.png` | Bar plots of sequence counts by segment, over time |
| `n.clean.samp.host.year.*.png` | Sample counts by host (thresholds: n≥10, n≥50), by year |
| `n.clean.samp.country.year.*.png` | Sample counts by country (thresholds: n≥10, n≥50), log and linear |
| `n.clean.samp.region.year.*.png` | Sample counts by global region over time |
| `n.clean.samp.FacetWrap-continent.host.year.1996.png` | Host × year faceted by continent (1996+) |
| `n.clean.samp.subtype.year.*.png` | Sample counts by H5Nx subtype |
| `prop.clean.samp.bysub.year.*.png` | Proportion of samples per subtype per year |
| `clean.coverage.bysegment.png` | Sequence coverage (% non-gap) by segment |
| `clean.width.bysegment.png` | Alignment width by segment |
| `clean.percGC.bysegment.png` | %GC content by segment |
| `clean.ambig.bysegment.png` | % ambiguous bases by segment |
| `Fig.cov.ambig.GC.bysegment.png` | Cowplot panel of the four alignment QC plots |
| `nuc.clean.align.summary.csv` | Per-sequence alignment statistics (coverage, GC, ambiguity) |
| `host_common_to_species.csv` / `.txt` | Each unique common host mapped to its species (with `has_species` flag) |
| `unique_species.csv` / `.txt` | Sorted unique species in the dataset (csv one-per-row; txt comma-separated) |
| `unique_host_common.csv` / `.txt` | Sorted unique common host names after normalization |

## Key host fix

The raw metadata `host` column incorrectly recoded several duck species as "Mallard" in `06_host/annotate_host_v2.Rmd`. This was corrected in `06_host/mallard_duck_fix.Rmd` (2026-06-18), which produced `complete_host_taxa_with_ecogroup_fixed.csv`. The fixed common names are propagated into `sample_taxa_lookup_all.csv` (from `06_host/multi-taxa/2026-06-19_out/`). In v3, the corrected `common` column from that lookup replaces the raw `host` column after joining on `sample`.

Affected host.old values that were wrongly mapped to Mallard:
- Chilean teal, Duck, Spot-billed duck, Steamer duck, Whistling duck

## Host common-name normalization (added 2026-06-29)

Two further local fixes are applied to `meta.all$host` in the `## Fix host` chunk, after the
duck/mallard correction, so all downstream host plots, tables, and reported counts use them:

1. **Latin → common** (`latin_common_fix.txt`, block format: common line / Latin line / blank
   separator). A few hosts were recorded under a scientific binomial instead of a common name;
   these are swapped to the common name (matched case-insensitively). The 4 fixed: Great knot
   (*Calidris tenuirostris*), Rufous hornero (*Furnarius rufus*), South American coati
   (*Nasua nasua*), Garganey (*Spatula querquedula*). For hosts whose common name is not keyed
   in the taxonomy lookup, the `# Host taxonomy summary` chunk backfills the species from the
   Latin name (canonical capitalized form).

2. **Common-name synonym merge** (`host_common_merge.csv`, `from` → `to`). 13 species were each
   referenced by >1 common name (synonyms / spelling / hyphen variants); each group is collapsed
   to the **most-sampled** name. 14 names merged away. Notable/lumped cases to revisit: Domestic
   ferret → European polecat; Velvet scoter → White-winged scoter; Starling → European starling
   (1-vs-1 tie). Edit `host_common_merge.csv` to change any canonical choice.

Resulting dataset counts (2026-06-29 inputs): **541 unique common hosts**, **442 unique
species** (442 hosts resolve to a species; 99 common names are absent from the taxonomy lookup).

## Bugs fixed in v3 (relative to v2)

1. **Country rename malformed**: South Georgia entry in `str_replace_all()` was an unnamed string, causing it to be erased rather than renamed.
2. **Silent name checks**: `all(names(nuc) %in% meta$seq_new)` was evaluated but not printed or asserted; replaced with `stopifnot()`.
3. **Segment by index**: `mixedsort(segments)[i]` assumed that file sort order matched segment name sort order; now extracted directly from filename with `str_extract()`.
4. **Host duck→mallard error**: Raw `host` column was never corrected in v2; v3 joins `sample_taxa_lookup_all.csv` and overwrites `host` with the fixed `common`.
5. **Input path**: v2 read from the phylo pipeline output directly, bypassing post-hoc fixes; v3 reads from `08_all-together/`.
6. **Tree type**: v2 used FastTree `.tre` files; v3 uses IQTREE2 `.treefile` (the final phylogenies).

## Known remaining issues

- `force_panelsizes(rows=c(1,1), cols=c(2.4,0.8,0.1))` in the continent facet plot is hardcoded for a 2×3 layout — will need adjustment if the number of continents in the data changes.
- **TODO (future 06_host build):** the Latin→common and synonym-merge fixes (`latin_common_fix.txt`, `host_common_merge.csv`) are applied locally here only. They should be adopted upstream in the 06_host taxonomy lookup (and thus propagate to multi-taxa, 08_all-together, 00_partition, and data-prep) so every stage shares one canonical common-name set. Deferred because these are all low-frequency hosts that get grouped together in the multi-taxa step anyway.
