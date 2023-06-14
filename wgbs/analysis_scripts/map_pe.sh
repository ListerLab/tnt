#!/bin/bash

### Script requires the following software in your PATH
# samtools
# sambamba
# fastp
# bbmap
# bs_seeker2-align.py
# bamUtils

### All of the file splitting in this script was optimied for our servers.
### The number of files and cores allocated appeared to be the fastest for the disk I/O and cores available on this system 

# If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

# Get the start time
start=`date +%s`

# Fastq input
readR1="$1"
readR2="$2"

# Number of cores. Must be a multiple of 6 as it is currently setup. 
cores="$3"

# Set number of files to be split into. In testing, 6 files ran the fastest, more slowed likely due to I/O limitations
files=6

constant=2

bt_cores=$(echo "($cores/$files)/$constant" | bc) &&

# Set base name for output
baseR1=$(basename "$readR1" .fastq.gz) && 
baseR2=$(basename "$readR2" .fastq.gz) &&

### Catch stdout
exec > >(tee -i "$baseR1"_logfile.txt)

#Name of the genome file you want to use
genome="hg19_L_PhiX.fa" # For human hg19 and chrL

#Path to the directory containing the BSSeeker index directory
indexBS="/home/sbuckberry/working_data_01/genomes/hg19_bsseeker_bt2_index/" # For human hg19

# Pre-trim reads report
fastp -Q -A -G -L --in1="$readR1" --in2="$readR2" --thread=16 \
--html="$baseR1"_fastp_report.html --json="$baseR1"_fastp_report.json &&

known_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

# Remove adapters using bbduk script in the bbmap toolkit
bbduk.sh \
in="$readR1" in2="$readR2" \
out="$baseR1"_trimmed.fastq.gz out2="$baseR2"_trimmed.fastq.gz \
literal="$known_adapter" threads="$cores" overwrite=t -Xmx20g \
overwrite=true &&

# Merge overlapping read pairs using the bbmerge script in the bbmap toolkit
bbmerge.sh in1="$baseR1"_trimmed.fastq.gz in2="$baseR2"_trimmed.fastq.gz qtrim=r \
out="$baseR1"_merged.fastq.gz outu1="$baseR1"_unmerged.fastq.gz outu2="$baseR2"_unmerged.fastq.gz \
-Xmx20g &&


##### Here is where we split the files for run time optimisation.

###### File split pairs and generate reports
fastp -A --split="$files" --split_prefix_digits=2 --thread="$files" \
--in1="$baseR1"_unmerged.fastq.gz --in2="$baseR2"_unmerged.fastq.gz \
--out1="$baseR1"_unmerged_split.fastq.gz --out2="$baseR2"_unmerged_split.fastq.gz \
--json "$baseR1"_unmerged_report.json --html "$baseR1"_unmerged_report.html &&

fastp -A --split="$files" --split_prefix_digits=2 --thread="$files" \
--in1="$baseR1"_merged.fastq.gz \
--out1="$baseR1"_merged_split.fastq.gz \
--json "$baseR1"_merged_report.json --html "$baseR1"_merged_report.html &&

############# Map the unmerged paired-end reads

# Setup the trimmed file manifest for mapping
ls 0{1..6}."$baseR1"_unmerged_split.fastq.gz > "$baseR1"_files &&
ls 0{1..6}."$baseR2"_unmerged_split.fastq.gz > "$baseR2"_files &&
sed 's/.fastq.gz/.bam/g' "$baseR1"_files > "$baseR1"_out &&
paste "$baseR1"_files "$baseR2"_files "$baseR1"_out > "$baseR1"_map_manifest &&

###### Paired-end alignment. The --split_line argument has not undergone optimisation testing. Larger number may be more efficient!
cat "$baseR1"_map_manifest | parallel --colsep="\t" -j"$files" bs_seeker2-align.py \
--aligner=bowtie2 --bt2--end-to-end --bt2-p "$bt_cores" -e 300 -X 2000 \
--temp_dir=/scratchfs/sbuckberry/tmp \
--split_line=4000000 \
-1 {1} -2 {2} -o {3} \
-d "$indexBS" \
-g "$genome" &&

###### Sort the output bam files
while read i; do sambamba sort -t "$cores" $i; done < "$baseR1"_out &&

### Clip the overlaps in the sorted BAM files using bamUtil
parallel -j"$files" bam clipOverlap --stats --in {} --out {.}_clipped.bam ::: 0{1..6}."$baseR1"_unmerged_split.sorted.bam &&

###### Merge the bam files
sambamba merge -t "$cores" "$baseR1"_pairs.bam 0{1..6}."$baseR1"_unmerged_split.sorted_clipped.bam &&

###### Create a merged log file 
cat 0{1..6}."$baseR1"_unmerged_split.bam.bs_seeker2_log > "$baseR1"_pairs_bs_seeker2.log &&

############# Map the merged pairs
ls 0{1..6}."$baseR1"_merged_split.fastq.gz > "$baseR1"_merged_fq && 
ls 0{1..6}."$baseR1"_merged_split.fastq.gz | sed 's/.fastq.gz/.bam/g' > "$baseR1"_merged_read_bams &&

paste "$baseR1"_merged_fq "$baseR1"_merged_read_bams | parallel -j"$files" --colsep="\t" bs_seeker2-align.py \
--aligner=bowtie2 --bt2--end-to-end --bt2-p "$bt_cores" -e 400 \
--temp_dir=/scratchfs/sbuckberry/tmp \
-i {1} -o {2} \
-d "$indexBS" \
-g "$genome" &&

###### Create a merged log file
cat 0{1..6}."$baseR1"_merged_split.bam.bs_seeker2_log > "$baseR1"_merged_pairs_bs_seeker2.log &&

###### Sort the output bam files
while read i; do sambamba sort -t "$cores" $i; done < "$baseR1"_merged_read_bams &&

###### Merge the single end map bam files
sambamba merge -t "$cores" "$baseR1"_merged_reads.bam 0{1..6}."$baseR1"_merged_split.sorted.bam &&

#### Merge all bam files
sambamba merge -t "$cores" "$baseR1".bam "$baseR1"_merged_reads.bam "$baseR1"_pairs.bam &&

# Create md5 sum for output bam file
md5sum "$baseR1".bam > "$baseR1".bam.md5 &&

# Clean up the temp files
#rm 0{1..6}."$baseR1"* \
#0{1..6}."$baseR2"* \
#"$baseR1"_trimmed.fastq.gz \
#"$baseR2"_trimmed.fastq.gz \
#"$baseR1"_merged.fastq.gz \
#"$baseR1"_unmerged.fastq.gz \
#"$baseR2"_unmerged.fastq.gz \
#"$baseR1"_files \
#"$baseR2"_files \
#"$baseR1"_out \
#"$baseR1"_map_manifest \
#"$baseR1"_merged_read_bams \
#"$baseR1"_merged_fq

end=`date +%s`

runtime=$((end-start))

echo $runtime > "$baseR1"_runtime_seconds.txt
