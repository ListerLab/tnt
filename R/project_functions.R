
library(preprocessCore)
library(Biobase)
library(caret)
library(e1071)
library(bsseq)
library(magrittr)
library(stringr)
library(openxlsx)
library(readxl)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(rtracklayer)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(stringr)
library(parallel)
library(pheatmap)
library(cowplot)
library(ggthemes)
#library(googlesheets4)
library(VariantAnnotation)
library(alluvial)
library(ggalluvial)
library(ggridges)
#library(AllelicImbalance)
#library(vcfR)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(ggbio)
#library(DSS)
library(ggdendro)
library(gtools)
library(UpSetR)
library(edgeR)
library(Gviz)
library(gt)
library(gprofiler2)
library(ggExtra)
library(XML)
library(RColorBrewer)


line_mm <- 0.5 / 2.835

vitC <- c("#004358", "#1695A3", "#BEDB39", "#FFE11A", "#1F8A70", "#FD7400")
#barplot(rep(1, times=length(vitC)), col=vitC)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#000000",
               "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

reprog_pal <- c(Primed="#009E73", Naive="#0072B2", HDF="#CC79A7",
                TNT="#EEBC4C", NtP="#55B3EA", "#5B7876", "#AFB4B7")

#barplot(rep(1, times=length(cbPalette)), col=cbPalette)

reprog_pal2 <- c(TNT="#eebc4cff", Primed="#009e73ff",
                 NtP="#55b3eaff", ESC="#a3a3a3ff",
                 Naive="#0072b2ff", HDF="#cc79a7ff")


#barplot(rep(1, times=length(reprog_pal2)), col=reprog_pal2)

calc_cv <- function(x) {100 * sd(x) / mean(x)}


standardise_mat <- function(mat){


        set_z <- function(x) {
                 z <- (mat[x, ] - mean(mat[x, ], na.rm=TRUE)) / sd(mat[x, ], na.rm=TRUE)
                 return(z)
                 }

        z_mat <- lapply(1:nrow(mat), set_z) %>% do.call(rbind, .)

        return(z_mat)
}


### ggplot theme for publication figures
sams_pub_theme <- function(legend_pos = "NA", x.text.angle = 45, hjust=1,
                           y.text.angle=0, line_point=0.5){
        line_mm <- line_point / 2.835
        custom_theme <- theme_bw() +
                theme(plot.background = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_line(size = line_mm),
                      panel.border = element_blank(),
                      axis.text.x = element_text(angle = x.text.angle, hjust = hjust, size = 6,
                                                 colour = 'black'),
                      axis.text.y = element_text(size=6, colour='black', angle = y.text.angle),
                      strip.text.y = element_text(size = 6),
                      text = element_text(size=6),
                      strip.background = element_blank(),
                      legend.position = legend_pos,
                      axis.line.x = element_line(color = 'black', size = line_mm),
                      axis.line.y = element_line(color = 'black', size = line_mm),
                      axis.ticks = element_line(color = 'black', size = line_mm))
        return(custom_theme)
}



plot_dendrogram <- function(df, title="", method="spearman",
                            hclust_method="average"){

        dissimilarity <- 1 - cor(df, method = method)
        distance <- as.dist(dissimilarity)
        hc <- hclust(distance, method = hclust_method)
        ggd <- ggdendrogram(hc, rotate = TRUE)
        ggd <- ggd + ggtitle(label = title)
        return(ggd)
}


# Correlation plot function
plot_correlation <- function(x, y, xlab="x", ylab="y", pch=19, lineCol="#004358",
                    legendPos="topleft"){
        fit <- lm(y~x)
        pValue <- signif(summary(fit)$coefficients["x","Pr(>|t|)"], digits=3)
        rSquared <- signif(summary(fit)$r.squared, digits=2)

        plot(x, y, pch=pch, xlab=xlab, ylab=ylab)
        abline(lm(y~x), col=lineCol, lwd=2)

        # Set up the legend
        rp <- vector('expression',2)
        rp[1] <- substitute(expression(italic(R)^2 == r2), list(r2 = format(rSquared, dig=3)))[2]
        rp[2] <- substitute(expression(italic(P) == pV), list(pV = format(pValue, digits = 2)))[2]
        legend(legendPos, rp, bty="n")
}

plot_pca <- function(mat, dim1=1, dim2=2, scale=TRUE){

        # Remove incomplete cases
        mat <- mat[complete.cases(mat), ]

        # Transpose
        mat <- t(mat)

        # Remove low variance features
        lowVar <- nearZeroVar(mat, saveMetrics = TRUE)
        message(str_c("Low variance features removed = ", sum(lowVar$nzv)))
        mat <- mat[ ,!lowVar$nzv]

        # Calculate PC's
        pr <- prcomp(x = mat, scale.=scale)
        pc1 <- (summary(pr)$importance[2, dim1] * 100) %>% round(digits = 1)
        pc2 <- (summary(pr)$importance[2, dim2] * 100) %>% round(digits = 1)

        pc1_dat <- pr$x[ ,dim1]
        pc2_dat <- pr$x[ ,dim2]
        samples <- rownames(pr$x)

        pca_df <- data.frame(Sample=samples, PC1=pc1_dat, PC2=pc2_dat)

        gg_pca <-  ggplot(data = pca_df,
                          mapping = aes(x = PC2, y = PC1, label=Sample)) +
                geom_point(alpha=0.8, size=4) +
                theme_linedraw() +
                theme(panel.grid = element_line(colour = 'grey')) +
                scale_color_manual(values = reprog_pal) +
                geom_text_repel(data = subset(pca_df, samples %in% samples),
                                point.padding = unit(1, "lines"), size=3) +
                xlab(str_c("PC", dim2, " (", pc2, "%)")) +
                ylab(str_c("PC", dim1, " (", pc1, "%)"))
        gg_pca
}



gr_to_loci <- function(gr){
        str_c(str_c(as.character(seqnames(gr)), start(gr), sep = ":") %>%
                      str_c(end(gr), sep = "-"))
}

# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path){
        stopifnot(file.exists(rds_path))
        message(str_c("Reading ", rds_path))
        message(Sys.time())
        bs_obj <- readRDS(file = rds_path)

        return(bs_obj)
}

# Calculate C coverage and methylation for ranges. Returns GRanges object
calc_mC_window <- function(rds_path, gr){

        message(str_c("Reading ", rds_path, " ", Sys.time()))
        bs_obj <- read_bs_obj(rds_path)

        message("Calculating coverage...")
        gr$Cov <- getCoverage(BSseq = bs_obj, regions = gr, type = "Cov",
                              what = "perRegionTotal")

        message("Calculating M...")
        gr$M <- getCoverage(BSseq = bs_obj, regions = gr, type = "M",
                            what = "perRegionTotal")

        message("Calculating methylation percentage...")
        gr$pc <- gr$M / gr$Cov

        return(gr)
}

# Calculate mC levels (%) for ranges from a list of bs_seq objects. Returns matrix.
make_mC_matrix <- function(obj_fls, gr, cores=1){

        grl <- mclapply(X = obj_fls, FUN = calc_mC_window,
                        gr = gr, mc.cores = cores) %>%
            GRangesList()

        loci <- str_c(as.character(seqnames(grl[[1]])),
                      start(grl[[1]]), sep = ":") %>%
                str_c(end(grl[[1]]), sep = "-")

        get_mC <- function(x){
                mC <- grl[[x]]$pc %>% c()
                return(mC)
        }

        message("Making matrix of mC levels for regions...")
        dat <- lapply(X = 1:length(grl), FUN = get_mC)
        dat <- do.call(cbind, dat)
        colnames(dat) <- basename(obj_fls)

        rownames(dat) <- loci

        return(dat)
}


#### Mfuzz functions====================================
acore <-  function(eset,cl,min.acore=0.5){

        atmp <- list(NULL)

        for (i in 1:dim(cl[[1]])[[1]]){
                index <- (cl[[3]]==i & (cl[[4]][,i] > min.acore)) # selection of genes
                atmp[[i]] <- data.frame(NAME=dimnames(exprs(eset))[[1]][index],
                                        MEM.SHIP=cl[[4]][index,i])
        }

        atmp
}

mfuzz <- function(eset,centers,m,...){

        cl<-cmeans(exprs(eset),centers=centers,method="cmeans",m=m,...)


}





standardise <- function(eset){
        data <- exprs(eset)
        for (i in 1:dim(data)[[1]]){
                data[i,] <- (data[i,] - mean(data[i,],na.rm=TRUE))/sd(data[i,],na.rm=TRUE)
        }
        exprs(eset) <- data
        eset
}

mestimate<- function(eset){
        N <-  dim(exprs(eset))[[1]]
        D <- dim(exprs(eset))[[2]]
        m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
        return(m.sj)
}




##### A function that imports a bed file to a GRanges object
# Only takes the first 3 columns from a bed file
bed_to_gr <- function(bed_path){
        dat <- fread(bed_path, sep = "\t", header = FALSE)
        gr <- GRanges(seqnames = dat$V1,
                      ranges = IRanges(start = dat$V2,
                                       end = dat$V3))
        return(gr)

}


# Read DMR's from DSS and HOME format to GRanges objects
read_dmr <- function(dmr_file){
        dmr <- fread(dmr_file)
        gr_dmr <- GRanges(seqnames = dmr$chr,
                          ranges = IRanges(start=dmr$start, end = dmr$end))
        return(gr_dmr)
}


# Add loci UCSC format to GRanges object
add_loci <- function(gr){
        gr$loci <- str_c(seqnames(gr), start(gr), sep = ":") %>%
                str_c(end(gr), sep = "-") %>% as.character()
        return(gr)
}

## Scale 0-1 function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

##### A function that writes a bed file from a genomic ranges object
# Only takes the first 3 columns from a bed file
gr_to_bed <- function(gr, out_path){
        dat <- as.data.frame(gr)[ ,1:3]
        dat$name <- "."
        dat$score <- "."
        dat$strand <- strand(gr)
        write.table(x = dat, file = out_path, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = FALSE)
}
