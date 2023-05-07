#!/bin/bash

bam=$1
bed=$2
intersectBed -abam "$bam" -b "$bed" | samtools view chrL | cut -f15 | sed 's/XM:Z://g' | sed 's/-//g' | sed '/^[[:space:]]*$/d' | pigz > "$bam"_naive_partial_mC_bins_5kb_read_calls.txt.gz

