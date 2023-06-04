# simulate sequencing reads from the Ig heavy chain locus for evaluating sterile/productive transcript distinction

Steps:

Before running alignments, use `hisat2_index.bash` and `star_index.bash` to index the chr14 genomic sequence for preparation of the alignment step. We downloaded the GRCh38 reference chr14 sequence and GTF files from [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index) (release 105) but you may choose to use the entire genome - just that it would be slower.

1. Use the `simulate_IGHC_read.Rmd` code to apply the polyester R package to sample reads.
2. Run `hisat2_alignment.bash` and `star_alignment.bash` to use HISAT2/STAR to align the sequences.
3. Run `add_tag.bash` to add origin cell/UMI barcodes to the BAM.
4. `get_headers.bash` to generate a text file which record the ground-truth location of sampling the read (this is stored in the FASTQ headers - this script takes them out and write out a separate text file). 
