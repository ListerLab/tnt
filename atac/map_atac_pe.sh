#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

R1="$1"
R2="$2"
prefix=$(echo "$R1" | sed 's/.fastq.gz//g')
cores=22
index="/home/sbuckberry/working_data_01/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"

# Trim the nextera adapters
sh /home/sbuckberry/working_data_01/bin/bbmap/bbduk.sh \
in="$R1" \
in2="$R2" \
out="$prefix".tmp_R1.fq \
out2="$prefix".tmp_R2.fq \
literal=GCGATCGAGGACGGCAGATGTGTATAAGAGACAG,CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG \
ktrim=r \
mink=3 \
threads="$cores" \
overwrite=true

## bowtie2 alignment
(bowtie2 -q --threads "$cores" -X2000 \
-x "$index" \
-1 "$prefix".tmp_R1.fq -2 "$prefix".tmp_R2.fq -S "$prefix".sam) 2> "$prefix".log

## Convert to BAM file and indexing 
samtools view -bS "$prefix".sam > "$prefix".bam
samtools sort "$prefix".bam "$prefix"_sorted.bam 

