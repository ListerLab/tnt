computeMatrix scale-regions -p 32 -S atac/atac_bigwigs/RL2436_coverage.bigwig \
atac/atac_bigwigs/RL2442_coverage.bigwig \
atac/atac_bigwigs/RL2448_coverage.bigwig \
atac/atac_bigwigs/RL2438_coverage.bigwig \
atac/atac_bigwigs/RL2454_coverage.bigwig \
atac/atac_bigwigs/RL2458_coverage.bigwig \
atac/atac_bigwigs/RL2437_coverage.bigwig \
atac/atac_bigwigs/RL2449_coverage.bigwig \
atac/atac_bigwigs/RL2457_coverage.bigwig \
atac/atac_bigwigs/RL2439_coverage.bigwig \
atac/atac_bigwigs/RL2443_coverage.bigwig \
atac/atac_bigwigs/RL2450_coverage.bigwig \
-R RNAseq/processed_data/MEL1_primed_down_gene_body.bed \
RNAseq/processed_data/MEL1_primed_up_gene_body.bed \
--beforeRegionStartLength 5000 \
--regionBodyLength 5000 \
--afterRegionStartLength 5000 \
--skipZeros \
-o atac/de_gene_atac_plots_matrix.mat.gz

plotHeatmap -m atac/de_gene_atac_plots_matrix.mat.gz \
-out atac/plots/de_gene_atac_heatmaps.pdf

git add atac/plots/de_gene_atac_heatmaps.pdf
git commit -m 'plot updates'
git pull && git push


computeMatrix scale-regions -p 32 -S atac/atac_bigwigs/RL2436_coverage.bigwig \
atac/atac_bigwigs/RL2442_coverage.bigwig \
atac/atac_bigwigs/RL2448_coverage.bigwig \
atac/atac_bigwigs/RL2438_coverage.bigwig \
atac/atac_bigwigs/RL2454_coverage.bigwig \
atac/atac_bigwigs/RL2458_coverage.bigwig \
atac/atac_bigwigs/RL2437_coverage.bigwig \
atac/atac_bigwigs/RL2449_coverage.bigwig \
atac/atac_bigwigs/RL2457_coverage.bigwig \
atac/atac_bigwigs/RL2439_coverage.bigwig \
atac/atac_bigwigs/RL2443_coverage.bigwig \
atac/atac_bigwigs/RL2450_coverage.bigwig \
-R atac/all_differential_peaks_primed_up.bed \
atac/all_differential_peaks_primed_down.bed \
--beforeRegionStartLength 3000 \
--regionBodyLength 1500 \
--afterRegionStartLength 3000 \
--skipZeros \
-o atac/de_atac_peaks_matrix.mat.gz

plotHeatmap -m atac/de_atac_peaks_matrix.mat.gz \
--heatmapWidth 2 \
--heatmapHeight 14 \
--colorMap 'Reds' \
-o atac/plots/de_atac_peaks_heatmaps.pdf

git add atac/plots/de_atac_peaks_heatmaps.pdf
git commit -m 'plot updates'
git pull && git push
