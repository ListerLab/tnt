#!/bin/bash

bam=$1
bed="/home/sbuckberry/working_data_02/polo_project/human_ips/resources/pastor_icr_regions.bed"

out_base=$(basename "$bam" .bam)

# Outputs chr, left-most mapping base, cigar string and mC call string

intersectBed -abam "$bam" -b "$bed" | samtools view | cut -f3,4,6,15 | sed 's/XM:Z://g' | awk -F'\t' '$3!=""' > "$out_base"_per_read_ICR_mC_calls.txt



