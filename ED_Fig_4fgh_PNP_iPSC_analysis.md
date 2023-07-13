Extended Data Figures 4f,g,h: PNP-hiPSC WGBS analysis
================
Sam Buckberry
2022-07-13

``` r
source("R/project_functions.R")
```

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: ggplot2

    ## Loading required package: lattice

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     subtract

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     shift

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     shift

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     shift

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, second

    ## Loading required package: BSgenome

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: rtracklayer

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'ggthemes'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     theme_map

    ## Loading required package: Rsamtools

    ## 
    ## Attaching package: 'VariantAnnotation'

    ## The following object is masked from 'package:stringr':
    ## 
    ##     fixed

    ## The following object is masked from 'package:base':
    ## 
    ##     tabulate

    ## 
    ## Attaching package: 'ChIPpeakAnno'

    ## The following object is masked from 'package:VariantAnnotation':
    ## 
    ##     info

    ## 
    ## Attaching package: 'gtools'

    ## The following object is masked from 'package:e1071':
    ## 
    ##     permutations

    ## 
    ## Attaching package: 'UpSetR'

    ## The following object is masked from 'package:lattice':
    ## 
    ##     histogram

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    ## Loading required package: grid

    ## 
    ## Attaching package: 'grid'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     pattern

Load the metadata

``` r
mdat <- read.csv("wgbs/metadata/wgbs_metadata_local.csv")
```

List relevant library ids

``` r
## Primed lines
pr <- c("RL417", "RL418")

## Primed-Naive lines
pn <- c("RL699", "RL700")

## Primed-Naive-Primed lines
pnp <- c("RL3073", "RL3074", "RL3075")

## ESC lines
esc <- c("RL2351_merge", "RL2352_merge",
         "SRR1561745_merge", "SRS004213", "SRS114877",
         "SRS606777", "SRS606778")

libs <- c(pr, pn, pnp, esc)
```

## Plot CG DMRs

Load DMRs

``` r
ips_esc_dmrs <- readRDS("wgbs/processed_data/classified_dmrs_granges.Rds")
```

Calculate mCG/CG for relevant samples

``` r
dmr_mCG <- make_mC_matrix(obj_fls = mdat$BSseq_CG[mdat$Library_id %in% libs],
                          gr = ips_esc_dmrs, cores = 3)
```

    ## Making matrix of mC levels for regions...

``` r
colnames(dmr_mCG) <- mdat$Library_id[mdat$Library_id %in% libs]
ind <- match(colnames(dmr_mCG), mdat$Library_id)

coldat <- data.frame(row.names = colnames(dmr_mCG), Group=mdat$Group[ind])

pdf("wgbs/plots/pnp-cg-dmr-heatmap.pdf",
    width = 4, height = 4)

pheatmap(dmr_mCG[complete.cases(dmr_mCG), ], show_rownames = FALSE,
         clustering_distance_cols = "correlation",
         annotation_col = coldat, fontsize = 6,
         labels_col = mdat$Manuscript.Name[ind])

dev.off()
```

    ## pdf 
    ##   3

``` r
pheatmap(dmr_mCG[complete.cases(dmr_mCG), ], show_rownames = FALSE,
         clustering_distance_cols = "correlation",
         annotation_col = coldat, fontsize = 6,
         labels_col = mdat$Manuscript.Name[ind])
```

![](ED_Fig_4fgh_PNP_iPSC_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Plot ICRs

``` r
icr_gr <- readRDS("resources/imprint_control_regions_granges_hg19.Rds")

icr_mCG <- make_mC_matrix(obj_fls = mdat$BSseq_CG[mdat$Library_id %in% c(pr, pn, pnp)],
                          gr = icr_gr, cores = 3)
```

    ## Making matrix of mC levels for regions...

``` r
colnames(icr_mCG) <- mdat$Library_id[mdat$Library_id %in% c(pr, pn, pnp)]


icr_mCG <- reshape2::melt(icr_mCG)

ind2 <- match(icr_mCG$Var2, mdat$Library_id)

icr_mCG$group <- mdat$Group[ind2]
icr_mCG$id <- mdat$Manuscript.Name[ind2]
icr_mCG$group <- factor(icr_mCG$group, levels = c("Primed-hiPSC",
                                                  "PtN-hiPSC",
                                                  "PNP-hiPSC"))

gg_icr <- ggplot(icr_mCG, aes(x = id, y = value, fill=group)) +
    geom_boxplot(lwd=line_mm) + ylab("ICR mCG/CG") + xlab("") +
    facet_grid(.~group, space = "free", scales = "free", drop = TRUE) + 
    sams_pub_theme()
```

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.

``` r
pdf("wgbs/plots/pnp-icr-boxplots.pdf", width = 2, height = 2)
gg_icr
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
gg_icr
```

![](ED_Fig_4fgh_PNP_iPSC_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
mdat <- read.csv("wgbs/metadata/wgbs_metadata_local.csv")
dfl <- readRDS("wgbs/processed_data/mCH_DMR_window_dat_all.Rds")

id_ind <- match(dfl$id, basename(mdat$BSseq_CA))

dfl$id <- mdat$Library_id[id_ind]

CH_dat <- readxl::read_excel(path = "resources/nature13551-s3.xlsx")
CH_dmr <- GRanges(seqnames = CH_dat$chr,
                  ranges = IRanges(start = as.numeric(CH_dat$start),
                                   end = as.numeric(CH_dat$end)))
CH_dmr$presence <- CH_dat$presence
CH_dmr$loci <- gr_to_loci(CH_dmr)

dfl$loci <- rownames(dfl)

dfl$loci <- rep(gr_to_loci(CH_dmr), times=length(unique(dfl$id)))

keep <- dfl$loci %in% CH_dmr$loci[CH_dmr$presence == "all iPSCs"]

table(keep)
```

    ## keep
    ## FALSE  TRUE 
    ##  7811  3139

``` r
dfl <- dfl[keep, ]

### Normalise to max (flank normalisation)
max_norm <- function(x){
    vec <- dfl[x, 1:30] / max(dfl[x, 1:30], na.rm = TRUE)
    return(vec)
}

#flank_max_norm <- function(x){
#    vec <- dfl[x, 1:30] / max(dfl[x, c(1:10, 21:30)], na.rm = TRUE)
#    return(vec)
#}

dfl_norm <- lapply(X = 1:nrow(dfl), max_norm) %>%
    do.call(rbind, .)

dfl_norm$id <- dfl$id
dfl_norm$loci <- dfl$loci

ind <- match(dfl_norm$id, mdat$Library_id)

dfl_norm$group <- mdat$Group[ind]
dfl_norm$progenitor <- mdat$Progenitor[ind]
dfl_norm$batch <- mdat$Batch[ind]
dfl_norm$lab <- mdat$Lab[ind]
dfl_norm$background <- mdat$Background[ind]


norm_dat <- reshape2::melt(dfl_norm)
```

    ## Using id, loci, group, progenitor, batch, lab, background as id variables

``` r
## Primed lines
pr <- c("RL417", "RL418")

## Primed-Naive lines
pn <- c("RL699", "RL700")

## Primed-Naive-Primed lines
pnp <- c("RL3073", "RL3074", "RL3075")

    df <- norm_dat[norm_dat$id %in% c(pr, pn, pnp), ]
    df <- rbind(df, norm_dat[norm_dat$group == "hESC", ])
    
    colnames(df)[colnames(df) == "variable"] <- "bin"
    
    plot_dat <- reshape2::melt(df)
```

    ## Using id, loci, group, progenitor, batch, lab, background, bin as id variables

``` r
    hm_dat <- plot_dat[plot_dat$bin %in% 11:20, ]

    hm_dat <- hm_dat %>% dplyr::group_by(loci, id) %>%
    dplyr::summarise(mean=mean(value, na.rm=TRUE))
```

    ## `summarise()` has grouped output by 'loci'. You can override using the
    ## `.groups` argument.

``` r
    hm_dat <- hm_dat %>% tidyr::spread(id, mean) %>% data.frame()
    rownames(hm_dat) <- hm_dat$loci
    hm_dat$loci <- NULL

    indx <- match(colnames(hm_dat), mdat$Library_id)
    coldat <- mdat[indx, c("Library_id", "Group")]
    rownames(coldat) <- coldat$Library_id
    coldat$Library_id <- NULL
    
    # hm <- pheatmap(hm_dat[complete.cases(hm_dat), ],
    #          annotation_col = coldat, main = title,
    #          border_color = NA, show_rownames = FALSE)

    ## Calculate stats on flank norm mean for CH-DMRs for each group

    stats_groups <- factor(coldat$Group)

    design <- model.matrix(~ 0+stats_groups)
    colnames(design) <- str_remove(string = colnames(design),
                                      pattern = "stats_groups") %>% make.names()

    cont <- makeContrasts(hESC - Primed.hiPSC,
                            levels = design)

    ch_fit <- lmFit(hm_dat, design)
    ch_fit2 <- contrasts.fit(ch_fit, cont)

    ch_fit2 <- eBayes(ch_fit2)

            tt <- topTable(ch_fit2, adjust="BH",
                      number = nrow(hm_dat), sort="none")
            tt$loci <- rownames(tt)
            tt$contrast <- dimnames(cont)$Contrasts
            tt$progenitor <- "Fibroblast"
            tt$batch <- "PNP"
            tt$lab <- "Lister"
            tt$background <- "32F"

    #tt <- lapply(1:ncol(cont), get_tt) %>% do.call(rbind, .)

    tt$diff <- "NS"
    tt$diff[(tt$logFC > 0) & (tt$adj.P.Val < 0.05)] <- "hiPSC-hypo"
    tt$diff[(tt$logFC < 0) & (tt$adj.P.Val < 0.05)] <- "hiPSC-hyper"
    
    ## Get the loci that are differential for up and down
    ips_hypo <- tt$loci[tt$diff == "hiPSC-hypo" & tt$contrast == "hESC - Primed.hiPSC"]
    ips_hyper <- tt$loci[tt$diff == "hiPSC-hyper" & tt$contrast == "hESC - Primed.hiPSC"]
    
    ## Flag loci direction for plotting
    df$direction <- "NS"
    df$direction[df$loci %in% ips_hyper] <- "hiPSC-hyper"
    df$direction[df$loci %in% ips_hypo] <- "hiPSC-hypo"

    ## Aggregate the data for plotting
    df$bin <- as.numeric(df$bin)
    
    df_grp <- df %>% dplyr::group_by(id, group, bin, direction) %>%
        dplyr::summarise(mean=mean(value, na.rm = TRUE))
```

    ## `summarise()` has grouped output by 'id', 'group', 'bin'. You can override
    ## using the `.groups` argument.

``` r
    scale_max <- function(x){ x / max(x, na.rm = TRUE)}
    
    scale_group_max <- function(id){

        df0 <- df_grp[df_grp$id == id, ]
        df0$norm <- NA
        df0$norm[df0$direction == "hiPSC-hyper"] <- 
            scale_max(df0$mean[df0$direction == "hiPSC-hyper"])
        df0$norm[df0$direction == "hiPSC-hypo"] <- 
            scale_max(df0$mean[df0$direction == "hiPSC-hypo"])
        df0$norm[df0$direction == "NS"] <- 
            scale_max(df0$mean[df0$direction == "NS"])

        return(df0)
    }              
    
    df_norm <- lapply(unique(df_grp$id), scale_group_max) %>%
        do.call(rbind, .)
```

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

    ## Warning in max(x, na.rm = TRUE): no non-missing arguments to max; returning -Inf

``` r
    # Set line width for plots
    line_mm <- 0.25
    
    pp01_keep <- ggplot(df_norm[df_norm$direction == "hiPSC-hypo", ],
                        aes(x=bin, y = norm, group=id,
                               fill=group, colour=group)) +
        geom_line(size=0.5, alpha = 0.7) +
        ylab("Normalised mCA/CA") +
    xlab("") +
    geom_vline(xintercept = c(11,20), linetype="dashed", size=line_mm) +
    theme_bw() + 
    theme(plot.background = element_blank(),
          panel.grid.minor.y = element_line(),
          panel.grid.major.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border = element_blank(),
          strip.text.y = element_text(size = 6),
          text = element_text(size=6),
          strip.background = element_blank(),
          legend.position = "right",
          axis.line.x = element_line(color = 'black', size = line_mm),
          axis.text.y = element_text(color = 'black'),
          axis.line.y = element_line(color = 'black', size = line_mm),
          axis.ticks.y = element_line(color = 'black', size = line_mm),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())    
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

``` r
pdf("wgbs/plots/pnp-ch-dmr-profile-plots.pdf", width = 3, height = 3)
pp01_keep
```

    ## Warning: Removed 180 rows containing missing values (`geom_line()`).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pp01_keep
```

    ## Warning: Removed 180 rows containing missing values (`geom_line()`).

![](ED_Fig_4fgh_PNP_iPSC_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
wb_ed_fig4fgh <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb_ed_fig4fgh, sheetName = "ED_Fig_4f")
openxlsx::writeData(wb = wb_ed_fig4fgh, sheet = "ED_Fig_4f",
                    x = dmr_mCG[complete.cases(dmr_mCG), ])

openxlsx::addWorksheet(wb_ed_fig4fgh, sheetName = "ED_Fig_4g")
openxlsx::writeData(wb = wb_ed_fig4fgh, sheet = "ED_Fig_4g",
                    x = pp01_keep$data)

openxlsx::addWorksheet(wb_ed_fig4fgh, sheetName = "ED_Fig_4h")
openxlsx::writeData(wb = wb_ed_fig4fgh, sheet = "ED_Fig_4h",
                    x = gg_icr$data)

openxlsx::saveWorkbook(wb = wb_ed_fig4fgh,
                       file = "ED_Figure_4fgh_source_data.xlsx", overwrite = TRUE)
```

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] RColorBrewer_1.1-3                     
    ##  [2] XML_3.99-0.12                          
    ##  [3] ggExtra_0.10.0                         
    ##  [4] gprofiler2_0.2.1                       
    ##  [5] gt_0.8.0                               
    ##  [6] Gviz_1.40.1                            
    ##  [7] edgeR_3.38.4                           
    ##  [8] limma_3.52.4                           
    ##  [9] UpSetR_1.4.0                           
    ## [10] gtools_3.9.4                           
    ## [11] ggdendro_0.1.23                        
    ## [12] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
    ## [13] ChIPpeakAnno_3.30.1                    
    ## [14] ggridges_0.5.4                         
    ## [15] ggalluvial_0.12.3                      
    ## [16] alluvial_0.1-2                         
    ## [17] VariantAnnotation_1.42.1               
    ## [18] Rsamtools_2.12.0                       
    ## [19] ggthemes_4.2.4                         
    ## [20] cowplot_1.1.1                          
    ## [21] ggrepel_0.9.2                          
    ## [22] ggfortify_0.4.15                       
    ## [23] pheatmap_1.0.12                        
    ## [24] GenomicFeatures_1.48.4                 
    ## [25] AnnotationDbi_1.58.0                   
    ## [26] BSgenome.Hsapiens.UCSC.hg19_1.4.3      
    ## [27] BSgenome_1.64.0                        
    ## [28] rtracklayer_1.56.1                     
    ## [29] Biostrings_2.64.1                      
    ## [30] XVector_0.36.0                         
    ## [31] data.table_1.14.6                      
    ## [32] readxl_1.4.1                           
    ## [33] openxlsx_4.2.5.1                       
    ## [34] stringr_1.5.0                          
    ## [35] magrittr_2.0.3                         
    ## [36] bsseq_1.32.0                           
    ## [37] SummarizedExperiment_1.26.1            
    ## [38] MatrixGenerics_1.8.1                   
    ## [39] matrixStats_0.63.0                     
    ## [40] GenomicRanges_1.48.0                   
    ## [41] GenomeInfoDb_1.32.4                    
    ## [42] IRanges_2.30.1                         
    ## [43] S4Vectors_0.34.0                       
    ## [44] e1071_1.7-12                           
    ## [45] caret_6.0-93                           
    ## [46] lattice_0.20-45                        
    ## [47] ggplot2_3.4.1                          
    ## [48] Biobase_2.56.0                         
    ## [49] BiocGenerics_0.42.0                    
    ## [50] preprocessCore_1.58.0                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rappdirs_0.3.3            ModelMetrics_1.2.2.2     
    ##   [3] R.methodsS3_1.8.2         tidyr_1.2.1              
    ##   [5] bit64_4.0.5               knitr_1.41               
    ##   [7] DelayedArray_0.22.0       R.utils_2.12.2           
    ##   [9] rpart_4.1.19              KEGGREST_1.36.3          
    ##  [11] hardhat_1.2.0             RCurl_1.98-1.9           
    ##  [13] AnnotationFilter_1.20.0   generics_0.1.3           
    ##  [15] lambda.r_1.2.4            RSQLite_2.2.19           
    ##  [17] proxy_0.4-27              future_1.29.0            
    ##  [19] bit_4.0.5                 xml2_1.3.3               
    ##  [21] lubridate_1.9.0           httpuv_1.6.6             
    ##  [23] assertthat_0.2.1          gower_1.0.0              
    ##  [25] xfun_0.35                 hms_1.1.2                
    ##  [27] evaluate_0.18             promises_1.2.0.1         
    ##  [29] fansi_1.0.4               restfulr_0.0.15          
    ##  [31] progress_1.2.2            dbplyr_2.2.1             
    ##  [33] DBI_1.1.3                 htmlwidgets_1.5.4        
    ##  [35] futile.logger_1.4.3       purrr_0.3.5              
    ##  [37] ellipsis_0.3.2            dplyr_1.0.10             
    ##  [39] backports_1.4.1           permute_0.9-7            
    ##  [41] biomaRt_2.52.0            deldir_1.0-6             
    ##  [43] sparseMatrixStats_1.8.0   vctrs_0.5.2              
    ##  [45] ensembldb_2.20.2          cachem_1.0.6             
    ##  [47] withr_2.5.0               checkmate_2.1.0          
    ##  [49] GenomicAlignments_1.32.1  prettyunits_1.1.1        
    ##  [51] cluster_2.1.4             lazyeval_0.2.2           
    ##  [53] crayon_1.5.2              labeling_0.4.2           
    ##  [55] recipes_1.0.3             pkgconfig_2.0.3          
    ##  [57] nlme_3.1-160              ProtGenerics_1.28.0      
    ##  [59] nnet_7.3-18               rlang_1.0.6              
    ##  [61] globals_0.16.2            lifecycle_1.0.3          
    ##  [63] miniUI_0.1.1.1            filelock_1.0.2           
    ##  [65] BiocFileCache_2.4.0       dichromat_2.0-0.1        
    ##  [67] VennDiagram_1.7.3         cellranger_1.1.0         
    ##  [69] graph_1.74.0              Matrix_1.5-3             
    ##  [71] Rhdf5lib_1.18.2           base64enc_0.1-3          
    ##  [73] png_0.1-8                 viridisLite_0.4.1        
    ##  [75] rjson_0.2.21              bitops_1.0-7             
    ##  [77] R.oo_1.25.0               rhdf5filters_1.8.0       
    ##  [79] pROC_1.18.0               blob_1.2.3               
    ##  [81] DelayedMatrixStats_1.18.2 regioneR_1.28.0          
    ##  [83] parallelly_1.32.1         jpeg_0.1-10              
    ##  [85] scales_1.2.1              memoise_2.0.1            
    ##  [87] plyr_1.8.8                zlibbioc_1.42.0          
    ##  [89] compiler_4.2.1            BiocIO_1.6.0             
    ##  [91] cli_3.6.0                 listenv_0.8.0            
    ##  [93] htmlTable_2.4.1           formatR_1.12             
    ##  [95] Formula_1.2-4             MASS_7.3-58.1            
    ##  [97] tidyselect_1.2.0          stringi_1.7.12           
    ##  [99] highr_0.9                 yaml_2.3.6               
    ## [101] locfit_1.5-9.6            latticeExtra_0.6-30      
    ## [103] tools_4.2.1               timechange_0.1.1         
    ## [105] future.apply_1.10.0       rstudioapi_0.14          
    ## [107] foreach_1.5.2             foreign_0.8-83           
    ## [109] gridExtra_2.3             prodlim_2019.11.13       
    ## [111] farver_2.1.1              digest_0.6.30            
    ## [113] shiny_1.7.3               lava_1.7.0               
    ## [115] Rcpp_1.0.9                later_1.3.0              
    ## [117] httr_1.4.4                biovizBase_1.44.0        
    ## [119] colorspace_2.1-0          splines_4.2.1            
    ## [121] RBGL_1.72.0               multtest_2.52.0          
    ## [123] plotly_4.10.1             xtable_1.8-4             
    ## [125] jsonlite_1.8.3            futile.options_1.0.1     
    ## [127] timeDate_4021.106         ipred_0.9-13             
    ## [129] R6_2.5.1                  Hmisc_4.7-2              
    ## [131] pillar_1.8.1              htmltools_0.5.3          
    ## [133] mime_0.12                 glue_1.6.2               
    ## [135] fastmap_1.1.0             BiocParallel_1.30.4      
    ## [137] class_7.3-20              codetools_0.2-18         
    ## [139] utf8_1.2.3                tibble_3.1.8             
    ## [141] curl_4.3.3                zip_2.2.2                
    ## [143] interp_1.1-3              survival_3.4-0           
    ## [145] rmarkdown_2.18            InteractionSet_1.24.0    
    ## [147] munsell_0.5.0             rhdf5_2.40.0             
    ## [149] GenomeInfoDbData_1.2.8    iterators_1.0.14         
    ## [151] HDF5Array_1.24.2          reshape2_1.4.4           
    ## [153] gtable_0.3.1
