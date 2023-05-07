source("R/project_functions.R")

tpm <- read.table("RNAseq/processed_data/timcourse_RNAseq_hg19_UCSC_tpm.txt")
mdat <- read.csv("RNAseq/rna-timecourse-sample-sheet.csv", header = TRUE)

tpm <- tpm[ ,mdat$id]

stopifnot(all(colnames(tpm)==mdat$id))

# Quantile normalise expression
#tpm <- normalizeBetweenArrays(log2(tpm+1), method = "quantile")

plot_gene <- function(geneID, y_min=NA, y_max=NA){

    #dat <- tpm[rownames(tpm) == geneID, ] %>% as.numeric()
    dat <- log2(tpm[rownames(tpm) == geneID, ]+1) %>% as.numeric()
    dat <- data.frame(timepoints=mdat$timepoint, tpm=dat, media=mdat$media)

    # Add in naive and primed d7 pseudo-timepoints for plotting
    dat2 <- data.frame(timepoints=c("d7", "d7"),
                      tpm=dat$tpm[dat$timepoints == "d7"],
                      media=c("Primed", "Naive"))
    dat <- rbind(dat2, dat)

    dat$timepoints <- str_replace(string = dat$timepoints,
                                  pattern = "^d", replacement = "day_")

    dat$timepoints <- factor(dat$timepoints, levels = c("day_0", "day_3", "day_7",
                                                        "day_13", "day_21", "P3",
                                                        "P10"))

    dat$media <- as.character(dat$media)
    dat$media[dat$media == "Fibroblast"] <- "Early"
    dat$media <- factor(dat$media, levels=c("Early", "Naive", "Primed"))

    gg <- ggplot2::ggplot(dat, aes(x=timepoints, group=media,
                                   col=media, fill=media, y=tpm)) +
        scale_x_discrete(limits=levels(dat$timepoints)) +
        scale_y_continuous(limits = c(y_min, y_max)) +
        geom_point(size=1) +
        stat_summary(mapping = aes(group=media, colour=media), geom="line",
                     fun = mean, size=0.5) +
        scale_color_manual(values = reprog_pal[c(3, 2, 1)]) +
        scale_fill_manual(values = reprog_pal[c(3, 2, 1)]) +
        ggtitle(geneID) +
        xlab("") + ylab("") +
        sams_pub_theme(legend_pos = "NA")
    gg
}

### Cluster genes

# Naive_2
p1 <- plot_gene(geneID = "LMNA") + ylab("Cluster 1")
p2 <- plot_gene(geneID = "NFIX")

# Primed_4
p3 <- plot_gene(geneID = "LGALS3") + ylab("Cluster 2")
p4 <- plot_gene(geneID = "COL6A3")

# Primed_1
p5 <- plot_gene(geneID = "STARD13") + ylab("Cluster 3")
p6 <- plot_gene(geneID = "KLF6")

# Primed_2 cluster
p7 <- plot_gene(geneID = "SNHG1") + ylab("Cluster 4")
p8 <- plot_gene(geneID = "TUBB2A")

# Primed_3 transient cluster
p9 <- plot_gene(geneID = "RNASET2") + ylab("Cluster 5")
p10 <- plot_gene(geneID = "ZMAT3")

# Plot all together
pdf("RNAseq/plots/Fig1_cluster_example_genes.pdf", height = 6, width = 2.6)
plot_grid(plotlist = list(p1, p2, p3, p4, p5,
                          p6, p7, p8, p9, p10),
          align = "hv", nrow = 5, ncol = 2)
dev.off()


### Pluripotency genes
pl_gene_plots <- list(plot_gene(geneID = "ANPEP", y_min = 0),
                    plot_gene(geneID = "SNAI2", y_min = 0),
                    plot_gene(geneID = "SERPINE1", y_min = 0),
                    plot_gene(geneID = "MMP1", y_min = 0),
                    plot_gene(geneID = "ESM1", y_min = 0),
                    plot_gene(geneID = "ZIC2", y_min = 0),
                    plot_gene(geneID = "SOX11", y_min = 0),
                    plot_gene(geneID = "ZIC3", y_min = 0),
                    plot_gene(geneID = "SFRP2", y_min = 0),
                    plot_gene(geneID = "SALL2", y_min = 0),
                    plot_gene(geneID = "KLF4", y_min = 0),
                    plot_gene(geneID = "KLF5", y_min = 0),
                    plot_gene(geneID = "KLF17", y_min = 0),
                    plot_gene(geneID = "TFCP2L1", y_min = 0),
                    plot_gene(geneID = "DPPA5", y_min = 0))

pdf("RNAseq/plots/Fig1_pluripotency_genes.pdf", height = 6, width = 3.9)
plot_grid(plotlist = pl_gene_plots, align = "hv", nrow = 5, ncol = 3,
          byrow = FALSE)
dev.off()

### DNA methylation genes
mC_gene_plots <- list(plot_gene(geneID = "DNMT1", y_min = 0),
                      plot_gene(geneID = "DNMT3A", y_min = 0),
                      plot_gene(geneID = "DNMT3B", y_min = 0),
                      plot_gene(geneID = "DNMT3L", y_min = 0),
                      plot_gene(geneID = "TET1", y_min = 0),
                      plot_gene(geneID = "TET2", y_min = 0),
                      plot_gene(geneID = "TET3", y_min = 0),
                      plot_gene(geneID = "DPPA3", y_min = 0))

pdf("RNAseq/plots/Fig1_methylation_genes.pdf", height = 5, width = 2.6)
plot_grid(plotlist = mC_gene_plots, align = "hv", nrow = 4, ncol = 2,
          byrow = FALSE)
dev.off()



