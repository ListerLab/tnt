#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

R1="$1"
R2="$2"

prefix=$(basename "$R1" .fastq.gz) && 
prefix2=$(basename "$R2" .fastq.gz) &&

cores=48
index="/home/sbuckberry/working_data_01/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"

# Adapter and quality trimming
fastp -i "$R1" -I "$R2" -o "$prefix"_trim.fastq.gz -O "$prefix2"_trim.fastq.gz

# Map using bowtie
(bowtie2 -q --threads "$cores" \
-x "$index" -X 2000 \
-1 "$prefix"_trim.fastq.gz -2 "$prefix2"_trim.fastq.gz | \
samtools view -bSu - | samtools sort -T "$prefix" - > "$prefix".bam) 2> "$prefix".log

# Create md5 sum for output bam file
md5sum "$prefix".bam > "$prefix".md5

# make the index
sambamba index "$prefix".bam

# Remove temp files
rm "$prefix"_trim.fastq.gz




