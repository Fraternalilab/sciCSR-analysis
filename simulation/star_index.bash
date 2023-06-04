#!/bin/bash

# generate index of chr14 for the STAR aligner

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles Homo_sapiens.GRCh38.dna.chromosome.14.fa --sjdbGTFfile Homo_sapiens.GRCh38.105.chr14.gtf --sjdbOverhang 99 
