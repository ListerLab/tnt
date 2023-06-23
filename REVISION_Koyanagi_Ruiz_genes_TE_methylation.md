Evaluate published criteria of hiPSC quality
================
Sam Buckberry
2023-06-23

Analysis of data in context of findings in

Koyanagi-Aoi, M. et al. Differentiation-defective phenotypes revealed by
large-scale analyses of human pluripotent stem cells. Proc. Natl. Acad.
Sci. U. S. A. 110, 20569–20574 (2013).

Ruiz, S. et al. Identification of a specific reprogramming-associated
epigenetic signature in human induced pluripotent stem cells. Proc.
Natl. Acad. Sci. U. S. A. 109, 16196–16201 (2012).

as per reviewer request

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

``` r
gr_te <- readRDS("resources/hg19_rmsk_TE_granges.Rds")

gtfPath <- "resources/genes.gtf.gz"
txdb <- makeTxDbFromGFF(file = gtfPath,
                        format = "gtf")
```

    ## Import genomic features from the file as a GRanges object ...

    ## OK

    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ... OK

``` r
gene_gr <- genes(txdb)
```

    ##   643 genes were dropped because they have exons located on both strands
    ##   of the same reference sequence or on more than one reference sequence,
    ##   so cannot be represented by a single genomic range.
    ##   Use 'single.strand.genes.only=FALSE' to get all the genes in a
    ##   GRangesList object, or use suppressMessages() to suppress this message.

``` r
gene_sub <- gene_gr[gene_gr$gene_id %in% c("HHLA1", "ABHD12B", "C4orf51")]

te_sub <- gr_te[overlapsAny(gr_te, gene_sub)]

te_sub <- te_sub[te_sub$gene == "LTR7"]

mdat <- read.csv("wgbs/metadata/wgbs_metadata_local.csv")
mdat <- mdat[mdat$Group %in% c("Primed-hiPSC", "TNT-hiPSC", "NtP-hiPSC", "hESC"), ]
mdat <- mdat[mdat$Media != "t2iLGoY", ]
mdat <- mdat[mdat$Progenitor %in% c("ESC", "Fibroblast"), ]
mdat <- mdat[mdat$Lab == "Lister", ]
mdat <- mdat[mdat$Batch != "C", ]
mdat <- mdat[!grepl("merge", mdat$Library_id), ]

file.exists(mdat$BSseq_CG)
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE

``` r
te_mCG <- make_mC_matrix(obj_fls = mdat$BSseq_CG, gr = te_sub, cores = 3)
```

    ## Making matrix of mC levels for regions...

``` r
colnames(te_mCG) <- mdat$Library_id
te_mCG <- data.frame(te_mCG)

te_mCG$gene <- c("ABHD12B", "ABHD12B", "C4orf51", "C4orf51", "HHLA1", "HHLA1")

te_mCG_melt <- reshape2::melt(te_mCG)
```

    ## Using gene as id variables

``` r
ind <- match(te_mCG_melt$variable, mdat$Library_id)

te_mCG_melt$group <- mdat$Group[ind]
te_mCG_melt$lab <- mdat$Lab[ind]
te_mCG_melt$background <- mdat$Background[ind]

te_mCG_melt <- te_mCG_melt[te_mCG_melt$lab == "Lister" &
                               (te_mCG_melt$group %in% c("Primed-hiPSC", "TNT-hiPSC", "NtP-hiPSC", "hESC")), ]

te_mCG_melt$group <- factor(te_mCG_melt$group,
                            levels=c("Primed-hiPSC", "TNT-hiPSC", "NtP-hiPSC", "hESC"))

te_mCG_melt$gene <- factor(te_mCG_melt$gene, levels=c("HHLA1", "ABHD12B", "C4orf51"))

pdf("wgbs/plots/Koyanagi_genes_te_mCG_boxplots.pdf", width = 3, height = 2)
ggplot(te_mCG_melt, aes(x = group, y = value, group=group)) +
    geom_boxplot() + ylab("LTR7 mCG/CG") +
    facet_grid(.~gene, scales = "free", space = "free", drop = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6,
                                     colour = 'black'),
          axis.text.y = element_text(size=6, colour='black', angle = 0),
          strip.text.y = element_text(size = 6),
          text = element_text(size=6),
          strip.background = element_blank(),
          axis.line.x = element_line(color = 'black', size = line_mm),
          axis.line.y = element_line(color = 'black', size = line_mm),
          axis.ticks = element_line(color = 'black', size = line_mm))
```

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
wb_ed_fig9c <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_ed_fig9c, sheetName = "ED_Fig_9c")
openxlsx::writeData(wb = wb_ed_fig9c, sheet = "ED_Fig_9c",
                    x = te_mCG_melt)
openxlsx::saveWorkbook(wb = wb_ed_fig9c,
                       file = "ED_Figure_9c_source_data.xlsx", overwrite = TRUE)


ruiz_genes <- c("PTPRT", "TMEM132C", "TMEM132D", "TCERG1L",
                "DPP6", "FAM19A5", "RBFOX1", "CSMD1", "C22orf34")


ruiz_gene_gr <- gene_gr[gene_gr$gene_id %in% ruiz_genes]

## Authors add 10kb to gene start and end
ruiz_gene_gr_expand <- ruiz_gene_gr + 10000

## Get mCG
ruiz_mCG <- make_mC_matrix(obj_fls = mdat$BSseq_CG,
                           gr = ruiz_gene_gr_expand, cores = 4)
```

    ## Making matrix of mC levels for regions...

``` r
colnames(ruiz_mCG) <- mdat$Library_id
ruiz_mCG <- data.frame(ruiz_mCG)

ruiz_mCG$gene <- as.character(ruiz_gene_gr$gene_id)

ruiz_mCG_melt <- reshape2::melt(ruiz_mCG)
```

    ## Using gene as id variables

``` r
ind2 <- match(ruiz_mCG_melt$variable, mdat$Library_id)

ruiz_mCG_melt$group <- mdat$Group[ind2]
ruiz_mCG_melt$lab <- mdat$Lab[ind2]
ruiz_mCG_melt$background <- mdat$Background[ind2]

ruiz_mCG_melt <- ruiz_mCG_melt[ruiz_mCG_melt$lab == "Lister" &
                               (ruiz_mCG_melt$group %in% c("Primed-hiPSC",
                                                           "TNT-hiPSC", "NtP-hiPSC", "hESC")), ]

ruiz_mCG_melt$group <- factor(ruiz_mCG_melt$group,
                            levels=c("Primed-hiPSC", "TNT-hiPSC",
                                     "NtP-hiPSC", "hESC"))

pdf("wgbs/plots/ruiz_gene_methylation_boxplots.pdf", width = 7, height = 2.5)
ggplot(ruiz_mCG_melt, aes(x = group, y = value, group=group, )) +
    geom_boxplot() + ylab("LTR7 mCG/CG") +
    facet_grid(.~gene, scales = "free", space = "free", drop = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6,
                                     colour = 'black'),
          axis.text.y = element_text(size=6, colour='black', angle = 0),
          strip.text.y = element_text(size = 6),
          text = element_text(size=6),
          strip.background = element_blank(),
          axis.line.x = element_line(color = 'black', size = line_mm),
          axis.line.y = element_line(color = 'black', size = line_mm),
          axis.ticks = element_line(color = 'black', size = line_mm))
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
wb_ed_fig9e <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_ed_fig9e, sheetName = "ED_Fig_9e")
openxlsx::writeData(wb = wb_ed_fig9e, sheet = "ED_Fig_9e",
                    x = ruiz_mCG_melt)
openxlsx::saveWorkbook(wb = wb_ed_fig9e,
                       file = "ED_Figure_9e_source_data.xlsx", overwrite = TRUE)
```
