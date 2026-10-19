# clean_ncbi (v2)

* Clean ncbi nucleotide sequences and metadata:
	- Inspect and clean metadata
	- Format genbank title and locations to exclude special characters
	- Parse information from GenBankTitle: strain/host/location/year(subtype) segment
	- Fill in missing collection dates from parsed IDs (n=6300)
	- Harmonize segment naming fill in missing from parsed IDs
	- Identify duplicate segments for each organism name, biosample, or assembly ID. 
		- Exclude most complete, shortest, then oldest duplicate, or at random
	- Exclude samples with incomplete segments
	- Exclude vaccine strains, passaged, clones, escape mutants
	- Exclude if no collection year after parsing
	- Summarize n, geo, host, subtype, date
	- Filter sequences to exclude meta removals
	- Split nucleotide, amino acid sequences into segments; export unaligned
	
* Clean ncbi protein sequences and metadata:
	- Parsed titles: accession, strain, protein, assembly, isolate, biosample
	- Initially 235 unique proteins, but spelling mistakes/abbr/…
	- Harmonized protein naming
	- Excluded proteins with no matching assembly, isolate, sample, or strain in nucleotide data
	- Removed duplicate proteins within a sample 
	- Renamed sequences and exported unaligned proteins