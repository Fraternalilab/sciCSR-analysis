#!/bin/bash

# get sequence headers from FASTQ files and output separately

for folder in $(ls -d case*/); do 
	name=${folder%/} 
	zcat $name/fiveprime_`echo $name`.fastq.gz | grep "@" > $name/fiveprime_`echo $name`.fastq_header
	zcat $name/threeprime_`echo $name`.fastq.gz | grep "@" > $name/threeprime_`echo $name`.fastq_header
done
