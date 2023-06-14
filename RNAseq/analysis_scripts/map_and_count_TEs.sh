parallel -j10 TEcount --sortByPos --format BAM --mode uniq --stranded reverse \
--project {.} \
--GTF /home/sbuckberry/working_data_02/polo_project/human_ips/resources/genes_and_ERCC92.gtf \
--TE /home/sbuckberry/working_data_02/polo_project/human_ips/resources/hg19_rmsk_TE.gtf \
-b {} ::: /home/sbuckberry/working_data_02/polo_project/human_ips/RNAseq/polyA/*_merged_lanes.sorted.bam