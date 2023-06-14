#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

# Fastq input
readR1="$1"
readR2="$2"
cores="$3"

# Set base name for output
baseName=$(basename "$readR1" .fastq.gz)

# Argument 5: Path to the index folder
indexPath="/home/sbuckberry/working_data_01/genomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2_ERCC_Index/hisat2_hg19_ercc_index"

# Trim the reads
fastp -i "$readR1" -I "$readR2" -o "$baseName"_trimmed_R1.fastq.gz -O "$baseName"_trimmed_R2.fastq.gz
# Align the reads using hisat2
hisat2 --time \
--rna-strandness RF \
--threads "$cores" \
-k 5 \
-x "$indexPath" \
-1 "$baseName"_trimmed_R1.fastq.gz -2 "$baseName"_trimmed_R2.fastq.gz \
-S "$baseName".sam

# convert to bam, sort and index alignment
samtools view -bSu "$baseName".sam | samtools sort -T "$baseName"_sorted - > "$baseName"_5MM.bam &&

samtools index "$baseName"_5MM.bam &
md5sum "$baseName"_5MM.bam > "$baseName".bam.md5

# Cleanup files
rm "$baseName"_trimmed_R1.fastq.gz "$baseName"_trimmed_R2.fastq.gz "$baseName".sam
