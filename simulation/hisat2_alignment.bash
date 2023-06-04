#!/bin/bash

# perform alignment using HISAT2

for folder in $(ls -d case*/); do
	name=${folder%/}
	echo $name
	echo "5' library ..."
	/docker_main/hisat2-2.2.1/hisat2 -x /ref_genome/hisat2_indices/chr14 \
		-U `echo $name`/fiveprime_`echo $name`.fastq.gz \
		-S `echo $name`/HISAT2_fiveprime_`echo $name`.sam
	samtools view -bT /ref_genome/Homo_sapiens.GRCh38.dna.chromosome.14.fa \
		`echo $name`/HISAT2_fiveprime_`echo $name`.sam > `echo $name`/HISAT2_fiveprime_`echo $name`.bam
	samtools sort `echo $name`/HISAT2_fiveprime_`echo $name`.bam `echo $name`/HISAT2_fiveprime_`echo $name`_sorted
	samtools index `echo $name`/HISAT2_fiveprime_`echo $name`_sorted.bam
	rm `echo $name`/HISAT2_fiveprime_`echo $name`.[bs]am
        echo "3' library ..."
        /docker_main/hisat2-2.2.1/hisat2 -x /ref_genome/hisat2_indices/chr14 \
                -U `echo $name`/threeprime_`echo $name`.fastq.gz \
                -S `echo $name`/HISAT2_threeprime_`echo $name`.sam
        samtools view -bT /ref_genome/Homo_sapiens.GRCh38.dna.chromosome.14.fa \
                `echo $name`/HISAT2_threeprime_`echo $name`.sam > `echo $name`/HISAT2_threeprime_`echo $name`.bam	
        samtools sort `echo $name`/HISAT2_threeprime_`echo $name`.bam `echo $name`/HISAT2_threeprime_`echo $name`_sorted
        samtools index `echo $name`/HISAT2_threeprime_`echo $name`_sorted.bam
        rm `echo $name`/HISAT2_threeprime_`echo $name`.[bs]am
done
