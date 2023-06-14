library(data.table)
library(R.utils)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(magrittr)
library(stringr)
library(GenomicRanges)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

CGmap <- args[1]

# Calc non-conversion stats
dat <- data.table::fread(input = paste('gzip -dc ', CGmap), select = c(1,2,3,5,7,8),
                         col.names = c("chr", "base", "position", "context", 
                                       "C_reads", "CT_reads"))

chrL <- dat[dat$chr == "chrL", ,drop=TRUE]

non_conv_rate <- (sum(chrL$C_reads)/sum(chrL$CT_reads))*100

nc_mCG <- sum(as.numeric(chrL[chrL$context == "CG", ]$C_reads)) / sum(as.numeric(chrL[chrL$context == "CG"]$CT_reads))*100
nc_mCA <- sum(as.numeric(chrL[chrL$context == "CA", ]$C_reads)) / sum(as.numeric(chrL[chrL$context == "CA"]$CT_reads))*100
nc_mCC <- sum(as.numeric(chrL[chrL$context == "CC", ]$C_reads)) / sum(as.numeric(chrL[chrL$context == "CC"]$CT_reads))*100
nc_mCT <- sum(as.numeric(chrL[chrL$context == "CT", ]$C_reads)) / sum(as.numeric(chrL[chrL$context == "CT"]$CT_reads))*100


chrL_stat <- round(data.frame(nc_mC = non_conv_rate,
                              nc_mCG = nc_mCG, nc_mCA = nc_mCA, nc_mCC = nc_mCC, nc_mCT = nc_mCT,
                              positions_covered = nrow(chrL),
                              median_x_coverage = median(chrL$CT_reads),
                              mC_min_pc = min((chrL$C_reads / chrL$CT_reads)*100),
                              mC_max_pc = max((chrL$C_reads / chrL$CT_reads)*100),
                              mean_mC_level_pc = mean((chrL$C_reads / chrL$CT_reads)*100), 
                              mC_std_dev = sd((chrL$C_reads / chrL$CT_reads)*100)), digits = 2)


regions <- GRanges(seqnames = seqnames(BSgenome.Hsapiens.UCSC.hg19),
                   ranges = IRanges(start = 1, end = seqlengths(BSgenome.Hsapiens.UCSC.hg19)))

chroms <- paste("chr", c(1:22, "X"), sep = "")
regions <- regions[seqnames(regions) %in% chroms]

nuc_counts <- getSeq(BSgenome.Hsapiens.UCSC.hg19, regions) %>%
        alphabetFrequency() %>% data.frame()

c_count <- sum(nuc_counts$C) + sum(nuc_counts$G)

dat <- dat[dat$chr %in% chroms, ,drop=TRUE]

covered_c <- nrow(dat) / c_count

c_summary <- summary(dat$CT_reads)
names(c_summary) <- str_c("Nuclear_genome_", names(c_summary)) %>%
        str_replace(pattern = " ", "_")

c_var <- var(dat$CT_reads)


# Caclulate per chromosome stats
calc_mC <- function(chr){
        
        keep <- dat$chr == chr
        chr_dat <- dat[keep, ]
        
        mCG <- sum(as.numeric(chr_dat[chr_dat$context == "CG", ]$C_reads)) / sum(as.numeric(chr_dat[chr_dat$context == "CG"]$CT_reads))*100
        mCA <- sum(as.numeric(chr_dat[chr_dat$context == "CA", ]$C_reads)) / sum(as.numeric(chr_dat[chr_dat$context == "CA"]$CT_reads))*100
        mCC <- sum(as.numeric(chr_dat[chr_dat$context == "CC", ]$C_reads)) / sum(as.numeric(chr_dat[chr_dat$context == "CC"]$CT_reads))*100
        mCT <- sum(as.numeric(chr_dat[chr_dat$context == "CT", ]$C_reads)) / sum(as.numeric(chr_dat[chr_dat$context == "CT"]$CT_reads))*100
        
        df <- data.frame(Chr=chr, mCG=mCG, mCA=mCA, mCC=mCC, mCT=mCT)
        return(df)
}

mC_chrom_dat <- lapply(X = chroms, FUN = calc_mC)
mC_chrom_dat <- do.call(rbind, mC_chrom_dat)

# Global statistics for autosomes
dat <- dat[dat$chr != "chrX", ]

# Calc mC percentages genome wide
mCG <- sum(as.numeric(dat[dat$context == "CG", ]$C_reads)) / sum(as.numeric(dat[dat$context == "CG"]$CT_reads))*100
mCA <- sum(as.numeric(dat[dat$context == "CA", ]$C_reads)) / sum(as.numeric(dat[dat$context == "CA"]$CT_reads))*100
mCC <- sum(as.numeric(dat[dat$context == "CC", ]$C_reads)) / sum(as.numeric(dat[dat$context == "CC"]$CT_reads))*100
mCT <- sum(as.numeric(dat[dat$context == "CT", ]$C_reads)) / sum(as.numeric(dat[dat$context == "CT"]$CT_reads))*100


df <- data.frame(chrL_non_conv_pc = non_conv_rate,
                 nc_mCG = nc_mCG, nc_mCA = nc_mCA, nc_mCC = nc_mCC, nc_mCT = nc_mCT,
                 chrL_positions_covered = nrow(chrL),
                 chrL_median_coverage = median(chrL$CT_reads),
                 mCG = mCG, mCA = mCA, mCC = mCC, mCT = mCT,
                 Nuclear_genome_covered_c_pc = covered_c * 100,
                 Nuclear_genome_min = as.numeric(c_summary[1]),
                 Nuclear_genome_1stQu = as.numeric(c_summary[2]),
                 Nuclear_genome_median = as.numeric(c_summary[3]),
                 Nuclear_genome_mean = as.numeric(c_summary[4]),
                 Nuclear_genome_3rdQu = as.numeric(c_summary[5]),
                 Nuclear_genome_max = as.numeric(c_summary[6]),
                 Nuclear_genome_cov_variance = c_var) %>% round(digits = 2)

write.table(x = df, file = str_c(basename(CGmap), "_postmap_stats.txt"), quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

write.table(x = mC_chrom_dat, file = str_c(basename(CGmap), "_mC_contig_stats.txt"), quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

