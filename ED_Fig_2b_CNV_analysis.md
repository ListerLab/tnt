Extended Data Figure 2b CNV analysis
================
Sam Buckberry
2023-06-23

``` r
source("R/project_functions.R")
```

## Experimental aim

The aim of this experiment is to determine if naive treatment of cells
during reprogramming causes any genomic instability, compared with
conventional primed reprogramming.

We performed Primed, TNT and NtP reprogramming for two independent
fibroblast lines (32F and 38F). To test for genomic deletions or
insertions, we extracted DNA and performed genome-wide genotyping using
the Illumina Global Screening array v2 with 663,320 SNPs surveyed. Data
were processed using GSGT Version 1.9.4.

    ##                 ID                                            File
    ## 1:      TNT_32F_P3      19C106621_KY1_MA1764-68_FinalReport.txt.gz
    ## 2:      TNT_38F_P3      19C106622_KY2_MA1764-72_FinalReport.txt.gz
    ## 3: NtP_38F_P12_P16      19C106623_KY3_MA1764-78_FinalReport.txt.gz
    ## 4:  Primed_32F_P33   19C108849_32FP33_MA1796-78_FinalReport.txt.gz
    ## 5:  Primed_38F_P22   19C108850_38FP22_MA1796-79_FinalReport.txt.gz
    ## 6: NtP_32F_P17_P23 19C108851_32FP1723_MA1796-85_FinalReport.txt.gz

## Data analysis

List the data files and check ids.

``` r
snp_fls <- list.files(path = "snp_array/",
                      pattern = "FinalReport.txt.gz",
                      full.names = TRUE)


if (all(basename(snp_fls) == ids)){
    names(snp_fls) <- names(ids)
} else {
    stop()
}
```

Function to read snp chip data and format for analysis

``` r
read_snp <- function(path, chroms=c(1:22, "X")){
    
    df <- fread(input = path, skip = 10, header = TRUE)
    
    colnames(df) <- str_replace_all(colnames(df),
                                    pattern = " ",
                                    replacement = "_")
    df <- df[df$Chr %in% chroms, ]
    
    ind <- match(basename(path), ids)
    
    df$Sample_ID <- names(ids)[ind]
    
    return(df)
}
```

Function to plot allele frequencies for each chromosome

``` r
plot_BAF <- function(snp_fls, chrom_list=chroms){
    
    df <- lapply(snp_fls, read_snp, chroms=chrom_list) %>% do.call(rbind, .)
    df$Chr <- factor(df$Chr, levels=c(1:22, "X", "Y"))
    
    gg_freq <- ggplot(data = df, aes(y = B_Allele_Freq, x = Position)) +
        geom_point(size=0.1, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        xlab("Position") + ylab("BAF") +
        facet_grid(Sample_ID~Chr, scales = "free", space = "free") +
        sams_pub_theme() +
        theme(panel.spacing.x = unit(0, "lines")) +
        theme(panel.border = element_rect(fill = NA))
    
    return(gg_freq)
}
```

Function to plot Log Ratio of alleles

``` r
plot_LR <- function(snp_fls, chrom_list=chroms){
    
    df <- lapply(snp_fls, read_snp, chroms=chrom_list) %>% do.call(rbind, .)
    df$Chr <- factor(df$Chr, levels=c(1:22, "X", "Y"))
    
    gg_log_ratio <- ggplot(df, aes(x = Position, y = Log_R_Ratio)) +
        geom_point(size=0.1, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0), limits = c(-3,3)) +
        #geom_smooth(size=0.5) +
        xlab("Position") + ylab("Log R ratio") +
        facet_grid(Sample_ID~Chr, scales = "free", space = "free") +
        sams_pub_theme() +
        theme(panel.spacing.x = unit(0, "lines")) +
        theme(panel.border = element_rect(fill = NA))
        
    return(gg_log_ratio)
}
```

Plot BAF

``` r
gg_baf <- plot_BAF(snp_fls = snp_fls, chrom_list = 1:22) +
    theme(axis.text.x = element_blank(), strip.text.y = element_text(angle = 0))
gg_baf
```

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Plot Log Ratios

``` r
gg_lr <- plot_LR(snp_fls = snp_fls, chrom_list = 1:22) +
    theme(axis.text.x = element_blank(), strip.text.y = element_text(angle = 0))
gg_lr
```

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Write the plots to PNG files to make publication figures. There are too
many data points for PDF files.

``` r
png("snp_array/all_snp_array_baf_plots.png", width = 12, height = 8,
    units = "in", res = 300)
gg_baf
```

    ## Warning: Removed 5875 rows containing missing values (`geom_point()`).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("snp_array/all_snp_array_log_rato_plots.png", width = 12, height = 8,
    units = "in", res = 300)
gg_lr
```

    ## Warning: Removed 5956 rows containing missing values (`geom_point()`).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

For each sample, plot the chromosomes separately

``` r
plot_LR_single_png <- function(x, chrom_list=chroms){
    
    snp_fl <- snp_fls[x]
    
    df <- read_snp(path = snp_fl, chroms = chrom_list)
    df$Chr <- factor(df$Chr, levels=c(1:22, "X", "Y"))

    gg_log_ratio <- ggplot(df, aes(x = Position, y = Log_R_Ratio)) +
        geom_point(size=0.1, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0), limits = c(-2,2)) +
        #geom_smooth(size=0.5) +
        xlab("Position") + ylab("Log R ratio") +
        facet_grid(Chr~., scales = "free", space = "free") +
        sams_pub_theme() +
        theme(panel.spacing.x = unit(0, "lines")) +
        theme(panel.border = element_rect(fill = NA))
    
    out <- str_c("snp_array/", names(snp_fls)[x], "_LR_contig_plots.png")
    
    ggsave(plot = gg_log_ratio, filename = out, device = "png",
           width = 12, height = 8, units = "in", dpi = 300, bg = "white")

}

lapply(X = 1:length(snp_fls), FUN = plot_LR_single_png)
```

    ## Warning: Removed 570 rows containing missing values (`geom_point()`).

    ## Warning: Removed 627 rows containing missing values (`geom_point()`).

    ## Warning: Removed 574 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1620 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1334 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1547 rows containing missing values (`geom_point()`).

    ## [[1]]
    ## [1] "snp_array/TNT_32F_P3_LR_contig_plots.png"
    ## 
    ## [[2]]
    ## [1] "snp_array/TNT_38F_P3_LR_contig_plots.png"
    ## 
    ## [[3]]
    ## [1] "snp_array/NtP_38F_P12_P16_LR_contig_plots.png"
    ## 
    ## [[4]]
    ## [1] "snp_array/Primed_32F_P33_LR_contig_plots.png"
    ## 
    ## [[5]]
    ## [1] "snp_array/Primed_38F_P22_LR_contig_plots.png"
    ## 
    ## [[6]]
    ## [1] "snp_array/NtP_32F_P17_P23_LR_contig_plots.png"

``` r
plot_BAFsingle_png <- function(x, chrom_list=chroms){
    
    snp_fl <- snp_fls[x]
    
    df <- lapply(snp_fl, read_snp, chroms=chrom_list) %>% do.call(rbind, .)
    df$Chr <- factor(df$Chr, levels=c(1:22, "X", "Y"))
    
    gg_freq <- ggplot(data = df, aes(y = B_Allele_Freq, x = Position)) +
        geom_point(size=0.1, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        xlab("Position") + ylab("BAF") +
        facet_grid(Chr~., scales = "free", space = "free") +
        sams_pub_theme() +
        theme(panel.spacing.x = unit(0, "lines")) +
        theme(panel.border = element_rect(fill = NA))
    
    out <- str_c("snp_array/", names(snp_fls)[x], "_BAF_contig_plots.png")
    
    ggsave(plot = gg_freq, filename = out, device = "png",
           width = 12, height = 8, units = "in", dpi = 300, bg = "white")
}

lapply(X = 1:length(snp_fls), FUN = plot_BAFsingle_png)
```

    ## Warning: Removed 549 rows containing missing values (`geom_point()`).

    ## Warning: Removed 601 rows containing missing values (`geom_point()`).

    ## Warning: Removed 561 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1596 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1311 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1522 rows containing missing values (`geom_point()`).

    ## [[1]]
    ## [1] "snp_array/TNT_32F_P3_BAF_contig_plots.png"
    ## 
    ## [[2]]
    ## [1] "snp_array/TNT_38F_P3_BAF_contig_plots.png"
    ## 
    ## [[3]]
    ## [1] "snp_array/NtP_38F_P12_P16_BAF_contig_plots.png"
    ## 
    ## [[4]]
    ## [1] "snp_array/Primed_32F_P33_BAF_contig_plots.png"
    ## 
    ## [[5]]
    ## [1] "snp_array/Primed_38F_P22_BAF_contig_plots.png"
    ## 
    ## [[6]]
    ## [1] "snp_array/NtP_32F_P17_P23_BAF_contig_plots.png"

Plot regions of interest

``` r
plot_snp_region <- function(snp_fl, chrom, range){
    
    df <- lapply(snp_fl, read_snp, chroms=chrom) %>% do.call(rbind, .)
    df <- df[df$Position >= as.numeric(range[1]) & df$Position <= as.numeric(range[2]), ]
    
    gg_freq <- ggplot(data = df, aes(y = B_Allele_Freq, x = Position)) +
        geom_point(size=0.2, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        xlab("Position") + ylab("BAF") +
        sams_pub_theme()
    gg_freq <- gg_freq + ggtitle(str_c(basename(snp_fl), "   chr", chrom, ":",
                                       as.integer(range[1]), "-",
                                       as.integer(range[2])))
    
    gg_log_ratio <- ggplot(df, aes(x = Position, y = Log_R_Ratio)) +
        geom_point(size=0.2, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        geom_smooth(size=0.5) +
        xlab("Position") + ylab("Log R ratio") +
        sams_pub_theme()
    
    gg <- cowplot::plot_grid(plotlist = list(gg_freq,gg_log_ratio), nrow = 2)
    return(gg)
}

pdf("snp_array/region_plots_chr6.pdf")
lapply(snp_fls, plot_snp_region, chrom=6, range=c(1e8, 1.25e8))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

    ## Warning: Removed 2 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 2 rows containing non-finite values (`stat_smooth()`).
    ## Removed 2 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 1 rows containing non-finite values (`stat_smooth()`).
    ## Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 3 rows containing non-finite values (`stat_smooth()`).
    ## Removed 3 rows containing missing values (`geom_point()`).

    ## Warning: Removed 6 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 6 rows containing non-finite values (`stat_smooth()`).
    ## Removed 6 rows containing missing values (`geom_point()`).

    ## Warning: Removed 12 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 12 rows containing non-finite values (`stat_smooth()`).
    ## Removed 12 rows containing missing values (`geom_point()`).

    ## Warning: Removed 7 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 7 rows containing non-finite values (`stat_smooth()`).
    ## Removed 7 rows containing missing values (`geom_point()`).

    ## $TNT_32F_P3

    ## 
    ## $TNT_38F_P3

    ## 
    ## $NtP_38F_P12_P16

    ## 
    ## $Primed_32F_P33

    ## 
    ## $Primed_38F_P22

    ## 
    ## $NtP_32F_P17_P23

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pdf("snp_array/region_plots_chr4.pdf")
lapply(snp_fls, plot_snp_region, chrom=4, range=c(0, 5e7))
```

    ## Warning: Removed 5 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 5 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 5 rows containing missing values (`geom_point()`).

    ## Warning: Removed 4 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 4 rows containing non-finite values (`stat_smooth()`).
    ## Removed 4 rows containing missing values (`geom_point()`).

    ## Warning: Removed 6 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 6 rows containing non-finite values (`stat_smooth()`).
    ## Removed 6 rows containing missing values (`geom_point()`).

    ## Warning: Removed 29 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 29 rows containing non-finite values (`stat_smooth()`).
    ## Removed 29 rows containing missing values (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 15 rows containing non-finite values (`stat_smooth()`).
    ## Removed 15 rows containing missing values (`geom_point()`).

    ## Warning: Removed 20 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'

    ## Warning: Removed 20 rows containing non-finite values (`stat_smooth()`).
    ## Removed 20 rows containing missing values (`geom_point()`).

    ## $TNT_32F_P3

    ## 
    ## $TNT_38F_P3

    ## 
    ## $NtP_38F_P12_P16

    ## 
    ## $Primed_32F_P33

    ## 
    ## $Primed_38F_P22

    ## 
    ## $NtP_32F_P17_P23

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Perform CNV detection

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("copynumber")
```

``` r
library(copynumber)
```

Make matricies of BAF and LRR

``` r
dat <- read_snp(path = snp_fls[1])

# measure must be one of "Log_R_Ratio" or "B_Allele_Freq"

make_mat <- function(file_list, measure="Log_R_Ratio"){
  
  dat_list <- lapply(file_list, read_snp)
  
  get_measure <- function(x){
    dat <- dat_list[[x]]
    val <- dat[ ,..measure]
    return(val)
    }
  
  dat_measure <- lapply(1:length(dat_list), get_measure)
  dat_measure <- do.call(cbind, dat_measure)
  colnames(dat_measure) <- names(file_list)
  
  loci <- data.frame(chr=dat_list[[1]]$Chr, pos=dat_list[[1]]$Position)
  
  df <- cbind(loci, dat_measure)
  
  df_sorted <- df[order(df$chr, df$pos, decreasing = FALSE), ]
  
  return(df_sorted)
}

mat_lrr <- make_mat(file_list = snp_fls, measure = "Log_R_Ratio")
mat_baf <- make_mat(file_list = snp_fls, measure = "B_Allele_Freq")
```

Filter rows with missing data

``` r
keep <- complete.cases(mat_lrr)
mat_lrr <- mat_lrr[keep, ]
mat_baf <- mat_baf[keep, ]
```

Normalisation using winsorize

``` r
lrr_win <- winsorize(mat_lrr, assembly = "hg19")
```

    ## winsorize finished for chromosome arm 1p 
    ## winsorize finished for chromosome arm 1q 
    ## winsorize finished for chromosome arm 10p 
    ## winsorize finished for chromosome arm 10q 
    ## winsorize finished for chromosome arm 11p 
    ## winsorize finished for chromosome arm 11q 
    ## winsorize finished for chromosome arm 12p 
    ## winsorize finished for chromosome arm 12q 
    ## winsorize finished for chromosome arm 13q 
    ## winsorize finished for chromosome arm 14q 
    ## winsorize finished for chromosome arm 15q 
    ## winsorize finished for chromosome arm 16p 
    ## winsorize finished for chromosome arm 16q 
    ## winsorize finished for chromosome arm 17p 
    ## winsorize finished for chromosome arm 17q 
    ## winsorize finished for chromosome arm 18p 
    ## winsorize finished for chromosome arm 18q 
    ## winsorize finished for chromosome arm 19p 
    ## winsorize finished for chromosome arm 19q 
    ## winsorize finished for chromosome arm 2p 
    ## winsorize finished for chromosome arm 2q 
    ## winsorize finished for chromosome arm 20p 
    ## winsorize finished for chromosome arm 20q 
    ## winsorize finished for chromosome arm 21p 
    ## winsorize finished for chromosome arm 21q 
    ## winsorize finished for chromosome arm 22q 
    ## winsorize finished for chromosome arm 3p 
    ## winsorize finished for chromosome arm 3q 
    ## winsorize finished for chromosome arm 4p 
    ## winsorize finished for chromosome arm 4q 
    ## winsorize finished for chromosome arm 5p 
    ## winsorize finished for chromosome arm 5q 
    ## winsorize finished for chromosome arm 6p 
    ## winsorize finished for chromosome arm 6q 
    ## winsorize finished for chromosome arm 7p 
    ## winsorize finished for chromosome arm 7q 
    ## winsorize finished for chromosome arm 8p 
    ## winsorize finished for chromosome arm 8q 
    ## winsorize finished for chromosome arm 9p 
    ## winsorize finished for chromosome arm 9q 
    ## winsorize finished for chromosome arm Xp 
    ## winsorize finished for chromosome arm Xq

``` r
#baf_win <- winsorize(mat_baf, assembly = "hg19")
```

``` r
plotGamma(mat_lrr, sample = 3, chrom=4, cex=3, cv = TRUE)
```

    ## cv progress:  20 % 
    ## cv progress:  40 % 
    ## cv progress:  60 % 
    ## cv progress:  80 % 
    ## cv progress:  100 %

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

    ## $gamma
    ##  [1]  10  20  30  40  50  60  70  80  90 100
    ## 
    ## $pred.error
    ##  [1] 80.83934 80.38810 80.63362 81.08527 81.86999 82.62941 83.09135 83.43638
    ##  [9] 83.79579 84.59579
    ## 
    ## $optGamma
    ## [1] 20

Call CNV segments

``` r
allele_seg <- aspcf(logR = lrr_win, BAF = mat_baf, kmin = 10, assembly = "hg19")
```

    ## aspcf finished for chromosome arm 1p 
    ## aspcf finished for chromosome arm 1q 
    ## aspcf finished for chromosome arm 10p 
    ## aspcf finished for chromosome arm 10q 
    ## aspcf finished for chromosome arm 11p 
    ## aspcf finished for chromosome arm 11q 
    ## aspcf finished for chromosome arm 12p 
    ## aspcf finished for chromosome arm 12q 
    ## aspcf finished for chromosome arm 13q 
    ## aspcf finished for chromosome arm 14q 
    ## aspcf finished for chromosome arm 15q 
    ## aspcf finished for chromosome arm 16p 
    ## aspcf finished for chromosome arm 16q 
    ## aspcf finished for chromosome arm 17p 
    ## aspcf finished for chromosome arm 17q 
    ## aspcf finished for chromosome arm 18p 
    ## aspcf finished for chromosome arm 18q 
    ## aspcf finished for chromosome arm 19p 
    ## aspcf finished for chromosome arm 19q 
    ## aspcf finished for chromosome arm 2p 
    ## aspcf finished for chromosome arm 2q 
    ## aspcf finished for chromosome arm 20p 
    ## aspcf finished for chromosome arm 20q

    ## aspcf is not run for TNT_32F_P3 in chromosome arm 21p because all of the BAF-values are outside the threshold values. Mean is returned for logR.

    ## aspcf is not run for Primed_32F_P33 in chromosome arm 21p because all of the BAF-values are outside the threshold values. Mean is returned for logR.

    ## aspcf is not run for Primed_38F_P22 in chromosome arm 21p because all of the BAF-values are outside the threshold values. Mean is returned for logR.

    ## aspcf is not run for NtP_32F_P17_P23 in chromosome arm 21p because all of the BAF-values are outside the threshold values. Mean is returned for logR.

    ## aspcf finished for chromosome arm 21p 
    ## aspcf finished for chromosome arm 21q 
    ## aspcf finished for chromosome arm 22q 
    ## aspcf finished for chromosome arm 3p 
    ## aspcf finished for chromosome arm 3q 
    ## aspcf finished for chromosome arm 4p 
    ## aspcf finished for chromosome arm 4q 
    ## aspcf finished for chromosome arm 5p 
    ## aspcf finished for chromosome arm 5q 
    ## aspcf finished for chromosome arm 6p 
    ## aspcf finished for chromosome arm 6q 
    ## aspcf finished for chromosome arm 7p 
    ## aspcf finished for chromosome arm 7q 
    ## aspcf finished for chromosome arm 8p 
    ## aspcf finished for chromosome arm 8q 
    ## aspcf finished for chromosome arm 9p 
    ## aspcf finished for chromosome arm 9q 
    ## aspcf finished for chromosome arm Xp 
    ## aspcf finished for chromosome arm Xq

Plot segmentation

``` r
plotAllele(logR = lrr_win, BAF = mat_baf, segments = allele_seg, sample=c(1:6),chrom=c(20), ylim=c(-1,1))
```

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in max(seg.lim2[2, ], na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-20-5.png)<!-- -->![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-20-6.png)<!-- -->

``` r
plotHeatmap(segments = allele_seg, upper.lim = 0.25)
```

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
pdf("snp_array/test.pdf", height = 16)
plotAberration(segments=allele_seg, thres.gain = 0.2,
               thres.loss = -0.25, layout = c(11,2), chrom = 1:22)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
allele_seg$width <- allele_seg$end.pos - allele_seg$start.pos
allele_seg[allele_seg$logR.mean < -0.25 & allele_seg$n.probes > 10, ]
```

    ##             sampleID chrom arm start.pos   end.pos n.probes logR.mean BAF.mean
    ## 2294  Primed_38F_P22    20   p  14761710  15368224      124   -0.3234   0.7028
    ## 2857 NtP_38F_P12_P16     4   p   5641753   7715181      967   -0.2863   0.8426
    ## 3758 NtP_32F_P17_P23     6   q 105609468 108281713      637   -0.2788   0.7763
    ##        width
    ## 2294  606514
    ## 2857 2073428
    ## 3758 2672245

Plot segment statistics

``` r
gg_seg_size <- ggplot(data = allele_seg, aes(x = log(width))) +
  geom_density() +
  facet_grid(sampleID~.)
gg_seg_size
```

    ## Warning: Removed 2 rows containing non-finite values (`stat_density()`).

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
gg_seg_lrr <- ggplot(data = allele_seg, aes(x = logR.mean)) +
  geom_density() +
  facet_grid(sampleID~.)
gg_seg_lrr
```

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
gg_seg_stat <- ggplot(data = allele_seg, aes(x = logR.mean, y = BAF.mean)) +
  geom_point() +
  facet_wrap(sampleID~.) +
  geom_vline(xintercept = c(-.2, .2)) +
  geom_hline(yintercept = 0.7) +
  xlim(c(-0.5, 0.5)) +
  ylim(c(0.5, 1))
gg_seg_stat
```

    ## Warning: Removed 4 rows containing missing values (`geom_point()`).

![](ED_Fig_2b_CNV_analysis_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

Plot InDel hits as matrix

``` r
# Set heatmap thresholds
lrr_threshold <- 0.25
baf_threshold <- 0.7
#snp_count <- 10

allele_seg_filt <- allele_seg[(abs(allele_seg$logR.mean) > lrr_threshold) &
                                (allele_seg$BAF.mean > baf_threshold), ]

hits_mat <- matrix(data = 0, nrow = length(unique(allele_seg$sampleID)),
                   ncol = length(unique(allele_seg$chrom)))

rownames(hits_mat) <- c("Primed_32F_P33", "Primed_38F_P22",
                        "TNT_32F_P3", "TNT_38F_P3",
                        "NtP_32F_P17_P23", "NtP_38F_P12_P16")

colnames(hits_mat) <- c(1:22, "X")

hit_samples <- allele_seg_filt$sampleID %>% unique()


for (x in 1:length(hit_samples)) {
  
  hits_mat[hit_samples[x],
           allele_seg_filt$chrom[allele_seg_filt$sampleID == hit_samples[x]]] <- 1
}

pdf("snp_array/indel_hit_heatmap.pdf", width = 6, height = 1.5)
pheatmap(hits_mat, color = c("white", "red"),
         cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE, angle_col = 0)
dev.off()
```

    ## pdf 
    ##   3

``` r
pheatmap(hits_mat, color = c("white", "red"),
         cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE, angle_col = 0)
```

Plot indel segments

``` r
#logR = lrr_win, BAF = mat_baf, segments = allele_seg

plot_snp_region <- function(segment=allele_seg_filt[1, ],
                            lrr=lrr_win, baf=mat_baf, flank_multiple=2){
    
    plot_chrom <- segment$chrom
    segment_width <- segment$end.pos - segment$start.pos
    plot_range <- c((segment$start.pos - (segment_width*flank_multiple)),
                      segment$end.pos + (segment_width*flank_multiple))
      
    lrr_sub <- lrr[(lrr$chrom == plot_chrom) &
                     (lrr$pos > min(plot_range)) &
                     (lrr$pos < max(plot_range)), ]
    
    lrr_sub <- reshape2::melt(lrr_sub, id.vars=c("chrom", "pos"))
    colnames(lrr_sub)[1] <- "chr"
    lrr_sub$measure <- "LRR"
    
    baf_sub <- baf[(baf$chr == plot_chrom) &
                     (baf$pos > min(plot_range)) &
                     (baf$pos < max(plot_range)), ]
      
    baf_sub <- reshape2::melt(baf_sub, id.vars=c("chr", "pos"))
    baf_sub$measure <- "BAF"
    
    df <- rbind(lrr_sub, baf_sub)
    df$variable <- factor(df$variable,
                          levels=c("Primed_32F_P33", "Primed_38F_P22",
                        "TNT_32F_P3", "TNT_38F_P3",
                        "NtP_32F_P17_P23", "NtP_38F_P12_P16"))
    
    gg <- ggplot(data = df, aes(y = value, x = pos)) +
        geom_point(size=0.1, alpha=0.5) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        facet_grid(measure~variable, scales = "free_y", space = "free") +
        xlab("Position") + 
        sams_pub_theme()
    
    return(gg)
}

gg1 <- plot_snp_region(segment = allele_seg_filt[1, ])
gg2 <- plot_snp_region(segment = allele_seg_filt[2, ])
gg3 <- plot_snp_region(segment = allele_seg_filt[3, ])

pdf("snp_array/indel_hit_plots.pdf", width = 8, height = 6)
cowplot::plot_grid(
  gg1, gg2, gg3, nrow = 3, align = "v")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Export source data for manuscript

``` r
all_dat <- rbind(gg1$data, gg2$data, gg3$data)

wb_ed_fig3b <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_ed_fig3b, sheetName = "ED_Fig_3b")
openxlsx::writeData(wb = wb_ed_fig3b, sheet = "ED_Fig_3b",
                    x = all_dat)
openxlsx::saveWorkbook(wb = wb_ed_fig3b,
                       file = "ED_Figure_3b_source_data.xlsx", overwrite = TRUE)
```
