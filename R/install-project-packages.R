if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

bioc_packages <-  c("Biobase", "preprocessCore", "bsseq",
                    "BSgenome.Hsapiens.UCSC.hg19",
                   "GenomicRanges", "GenomicFeatures", "Biostrings",
                   "rtracklayer", "VariantAnnotation", "ChIPpeakAnno",
                   "TxDb.Hsapiens.UCSC.hg19.knownGene", "edgeR", "limma",
                   "Gviz", "gt", "Glimma")

cran_packages <- c("statmod", "caret", "e1071", "magrittr", "stringr",
                   "data.table", "gprofiler2",
                   "pheatmap", "ggplot2", "ggfortify", "ggrepel",
                   "parallel", "cowplot", "ggthemes", "alluvial",
                   "ggalluvial", "readxl", "ggExtra",
                   "ggridges", "ggdendro", "gtools", "UpSetR", "gt", "xlsx")

BiocManager::install(pkgs = cran_packages, site_repository = "CRAN")
BiocManager::install(pkgs = bioc_packages)

