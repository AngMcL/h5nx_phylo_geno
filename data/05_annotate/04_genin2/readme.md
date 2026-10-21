## Genin2

genin2 version 2.1.6, by Alessandro Sartori (asartori@izsvenezie.it)
https://github.com/izsvenezie-virology/genin2

Requires: sequence name starts with the > character and ends with an undersore (_) followed by the name of the segment, e.g.:
>any_text|any_string/seq_name_PB1
                             ^^^^
cite the GitHub repository:
https://github.com/izsvenezie-virology/genin2

- pip install genin2
- Input data were formatted in 05_annotate/01_flumut, which used as input the cleaned alignments (filtered of low qual, aligned, and trimmed)
- Data was formatted with <sample>_<segment> and output in batches of 20,000 sequences
- Copied the entire output folder from flumut/01_format into the genin2 folder and then renamed in input_genin2
- Renamed the files to replace flumut with genin in the fasta and meta using a shell script
- then ran the iterate_genin.sh script across the batches of 20,000 seqs



