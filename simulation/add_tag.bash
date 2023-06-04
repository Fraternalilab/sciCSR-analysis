#!/bin/bash

# add tag of sample ID to BAM entries

add_tag () {
        samtools view -H $1 > `echo ${1%.bam}`_tagged.sam
        samtools view $1 | awk 'BEGIN{FS="\t"}; { split($1, a, "/"); print $0"\tCB:Z:"a[2]"\tUB:Z:"a[1] }' >> `echo ${1%.bam}`_tagged.sam
        samtools view -bT /ref_genome/Homo_sapiens.GRCh38.dna.chromosome.14.fa `echo ${1%.bam}`_tagged.sam > `echo ${1%.bam}`_tagged.bam
        samtools index `echo ${1%.bam}`_tagged.bam
}

for folder in $(ls -d case*/); do
	name=${folder%/}
	cd $folder
	echo $name
	add_tag STARaln_fiveprime_`echo $name`_sorted.bam
	add_tag STARaln_threeprime_`echo $name`_sorted.bam
        add_tag HISAT2_fiveprime_`echo $name`_sorted.bam
        add_tag HISAT2_threeprime_`echo $name`_sorted.bam
	rm *.sam
	cd ../
done
