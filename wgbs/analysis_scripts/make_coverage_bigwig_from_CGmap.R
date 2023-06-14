#!/usr/bin/env Rscript

library(bsseq)
library(magrittr)
library(stringr)
library(data.table)
library(R.utils)
library(rtracklayer)
library(GenomicRanges)

# Load library for hg19. This will need to change for other genome assemblies
library(BSgenome.Hsapiens.UCSC.hg19)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

# Path to a CGmap file from BSseeker2
path <- args[1]

# Outfile path
#out_path <- args[2]

# Out path
file_ext <- str_c("_c_cov.bigwig")
out_path <- str_replace(string = path, pattern = ".CGmap.gz$",
                        replacement = file_ext)


#---- Check inputs
stopifnot(file.exists(path))

# Function to create bigwig file
make_coverage_bigiwg <- function(path, out_path){
    
    message(str_c("Reading ", path, "..."))
    
    # Read the data
    dat <- data.table::fread(path, header = FALSE,
                             select = c(1,3,8), sep = "\t",
                             col.names = c("chr", "position", "CT_reads"))
    gr <- GenomicRanges::GRanges(seqnames = dat$chr, 
                                 ranges = IRanges(dat$position, dat$position))
    
    gr$score <- dat$CT_reads

    # Clean up to free up memory
    rm(dat)
    
    #Add seqinfo
    
    hg19_info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
    gr <- gr[seqnames(gr) %in% hg19_info@seqnames]
    gr <- keepSeqlevels(gr, value = hg19_info@seqnames[hg19_info@seqnames %in% seqnames(gr)])
    info <- intersect(seqinfo(gr), hg19_info)    
    seqinfo(gr) <- info
    
    # Write the bigwig file
    message(str_c("Saving ", out_path))
    rtracklayer::export.bw(object = gr, con = out_path)
    
    message("Done!")
}

# Execute the function with specified arguments
make_coverage_bigiwg(path = path, out_path = out_path)
