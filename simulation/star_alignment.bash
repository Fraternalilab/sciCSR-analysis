#!/bin/bash

# perform alignment using STAR

for folder in $(ls -d case*); do
	cd $folder
	name=${folder%/}
	echo $name
	echo "5' library ..."
	STAR --genomeDir /ref_genome/star_indices --outSAMmultNmax -1 --runThreadN 4 \
		--readNameSeparator space --outSAMunmapped Within KeepPairs \
		--outSAMtype BAM SortedByCoordinate --readFilesIn fiveprime_`echo $name`.fastq.gz \
		 --outFileNamePrefix STARaln_fiveprime_`echo $name`_ --readFilesCommand zcat
	mv STARaln_fiveprime_`echo $name`_Aligned.sortedByCoord.out.bam STARaln_fiveprime_`echo $name`_sorted.bam
	samtools index STARaln_fiveprime_`echo $name`_sorted.bam
        echo "3' library ..."
        STAR --genomeDir /ref_genome/star_indices --outSAMmultNmax -1 --runThreadN 4 \
                --readNameSeparator space --outSAMunmapped Within KeepPairs \
                --outSAMtype BAM SortedByCoordinate --readFilesIn threeprime_`echo $name`.fastq.gz \
                 --outFileNamePrefix STARaln_threeprime_`echo $name`_ --readFilesCommand zcat
        mv STARaln_threeprime_`echo $name`_Aligned.sortedByCoord.out.bam STARaln_threeprime_`echo $name`_sorted.bam
        samtools index STARaln_threeprime_`echo $name`_sorted.bam
	cd ../
done
