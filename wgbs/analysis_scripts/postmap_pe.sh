#!/bin/bash

# REquired tools to be installed in PATH
# paleomix (https://paleomix.readthedocs.io/en/latest/installation.html)
# sambamba
# samtools
# CGmap tools
# 
# Note that all contigs require the chr prefix
# 

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

# Argument 1: Path to a sorted BAM file
file="$1"

prefix=$(basename "$file" .bam)

# Argument 2: allocated cores
cores="$2"

# Path to genome fasta file (must be uncompressed)
genome="/home/sbuckberry/working_data_01/genomes/hg19_bsseeker_bt2_index/hg19_L_PhiX.fa"

# Get the chromosome/contig id's for splitting. grep -V avoids bug with samtools output
samtools idxstats "$file" | cut -f 1 | grep -v \* | grep chr > "$prefix".chromNames

# Split the bamfile into chromosome specific files
cat "$prefix".chromNames | parallel -j$cores samtools view -b -o "$prefix".{1}.temp.bam "$file" {1}

# Split the BAM file into single-end and paired end reads for removing duplicates
for i in "$prefix".chr*.temp.bam;
	do sambamba view -t "$cores" -f bam -F "paired" -o "$i"_paired.bam "$i";
done

for i in "$prefix".chr*.temp.bam;
	do sambamba view -t "$cores" -f bam -F "not paired" -o "$i"_singletons.bam "$i";
done

## Deduplicate the single-end merged reads with 5' and 3' matching ends
parallel -j"$cores" paleomix rmdup_collapsed --remove-duplicates {} ">" {}_dedup.bam ::: "$prefix".chr*.temp.bam_singletons.bam

# Remove the PCR duplicates for the paired reads
for i in "$prefix".chr*.temp.bam_paired.bam;
	do sambamba markdup -r -t "$cores" -p --tmpdir=/scratchfs/sbuckberry/tmp "$i" "$i"_dedup.bam;
done

# Merge the deduplicated data for each chromosome
cat "$prefix".chromNames | parallel -j1 sambamba merge -t "$cores" "$prefix"_{1}_all_dedup_temp.bam "$prefix".{1}.temp.bam_paired.bam_dedup.bam "$prefix".{1}.temp.bam_singletons.bam_dedup.bam

# Make a final BAM file of deduplicated files
sambamba merge -t "$cores" "$prefix"_dedup.bam "$prefix"_chr*_all_dedup_temp.bam

# Call methylation for each chromosome (parallel)
parallel -j "$cores" /home/sbuckberry/working_data_01/bin/cgmaptools/cgmaptools convert bam2cgmap \
--bam {} --genome "$genome" -o {.} ::: "$prefix"_chr*_all_dedup_temp.bam

# Recombine the output files and clean up
ls "$prefix"_chr*_all_dedup_temp.ATCGmap.gz | xargs zcat > "$prefix".ATCGmap && pigz -f -p $cores "$prefix".ATCGmap
ls "$prefix"_chr*_all_dedup_temp.CGmap.gz | xargs zcat > "$prefix".CGmap && pigz -f -p $cores "$prefix".CGmap

# clean up postmap temp files
rm "$prefix"*temp* \
"$prefix".chromNames

# Calculate stats
Rscript /home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/analysis_scripts/calc_wgbs_stats.R \
"$file" "$prefix"_dedup.bam "$prefix".CGmap.gz &

Rscript /home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/analysis_scripts/calc_mC_stats.R \
"$prefix".CGmap.gz

# Make bigwig files
Rscript /home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/analysis_scripts/make_browser_tracks.R \
"$prefix".CGmap.gz "$prefix"
