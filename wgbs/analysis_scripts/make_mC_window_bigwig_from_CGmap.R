#!/usr/bin/env Rscript

library(bsseq)
library(magrittr)
library(stringr)
library(data.table)
library(R.utils)
library(rtracklayer)
library(GenomicRanges)
library(tools)

# Load library for hg19. This will need to change for other genome assemblies
library(BSgenome.Hsapiens.UCSC.hg19)

# Do not output scientific notation
options(scipen=999)

# Set command line arguments
args <- commandArgs(TRUE)

# Path to a CGmap file from BSseeker2
cgmap_path <- args[1]

# Function to create bigwig file
CGmap_to_window_bigwig <- function(path, out_folder=".", context="CA",
                                   win_size=5000, win_step=1000,
                                   nc_correct=TRUE, nc_contig="chrL",
                                   contigs=str_c("chr", c(1:22, "X")),
                                   bs_genome="BSgenome.Hsapiens.UCSC.hg19"){
    
    #---- Check inputs
    stopifnot(file.exists(path))
    
    # Ensure context is upper case
    context <- casefold(context, upper = TRUE)
    stopifnot(context %in% c("CG", "CA", "CC", "CT", "CH"))
    
    #---- Setup out path
    file_ext <- str_c("_", context, "_", "Window",
                      win_size, "_Step", win_step, ".bigwig")
    
    out_path <- str_c(file_path_sans_ext(path), file_ext)
    
    out_path <- str_c(out_folder, "/", basename(out_path))
    
    #out_path <- str_replace(string = path, pattern = ".CGmap.gz$",
    #                        replacement = file_ext)
    
    #---- Read CGmap file to bsseq object and subset for context
    
    # Read the data
    message(str_c("Decompressing and reading ", path))
    dat <- data.table::fread(path, header = FALSE,
                             select = c(1,2,3,5,7,8), sep = "\t",
                             col.names = c("chr", "base", "position",
                                           "diContext", "C_reads", "CT_reads"))
    
    # subset data by context
    message(str_c("Subsetting data for dinucleotide context ", context, "..."))
    
    keep <- NA
    
    if (context == "CH"){
        keep <- dat$diContext != "CG"
    } else {
        keep <- dat$diContext == context
    }
    
    dat <- subset(dat, keep)
    
    # Add strand info
    dat$base <- ifelse(test = dat$base == "C", yes = "+", no = "-")
    
    #Load data into BSseq object
    bs_obj <- bsseq::BSseq(chr = dat$chr, pos=dat$position,
                           M = matrix(dat$C_reads),
                           Cov = matrix(dat$CT_reads))
    
    # Add strand to bsseq object
    bsseq::strand(bs_obj) <- dat$base
    
    rm(dat); gc()
    
    # Collapse strand if context is CG.
    # This helps reduce file size through aggregating symmetrical CGs
    if (context=="CG"){
        message("Collapsing strands for CG context...")
        bs_obj <- bsseq::strandCollapse(bs_obj)
    }
    
    #---- Make windows
    message("Making windows for genome...")
    gr_genome <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg19))
    gr_genome <- gr_genome[seqnames(gr_genome) %in% contigs]
    
    gr_windows <- slidingWindows(x = gr_genome, width = win_size, step = win_step)
    gr_windows <- unlist(gr_windows)

    # Only keep windows that have coverage
    gr_windows <- gr_windows[overlapsAny(gr_windows, bs_obj@rowRanges)]
    
    #---- Calculate methylation for windows per chromosome
    # Without looping through chromosomes, this would eat a lot of memory
    
    message("Calculating methylation levels for windows...")
    calc_mC_window_contig <- function(contig){
        
        message(contig)
        
        bs_sub <- chrSelectBSseq(BSseq = bs_obj, seqnames = contig)
        
        gr_sub <- gr_windows[seqnames(gr_windows) == contig]
        gr_sub <- dropSeqlevels(gr_sub, (!seqlevels(gr_sub) %in% contig))
        
        Cov <- getCoverage(BSseq = bs_sub, regions = gr_sub, type = "Cov",
                              what = "perRegionTotal")
        
        M <- getCoverage(BSseq = bs_sub, regions = gr_sub, type = "M",
                            what = "perRegionTotal")
        
        gr_sub$pc <- M / Cov
        
        gr_sub <- gr_sub[!is.na(as.vector(Cov))]

        return(gr_sub)
        
    }
    
    gr_mc <- lapply(contigs, calc_mC_window_contig) %>% GRangesList() %>% unlist()

    #---- Calculate non-conversion rate and subtract from windows
    # Note that this calculates the non-conversion for the context only
    
    correct_nc <- function(gr){
        
        nc_obj <- chrSelectBSseq(BSseq = bs_obj, seqnames = nc_contig)
        nc_cov <- getCoverage(BSseq = nc_obj, type = "Cov",
                              what = "perBase")
        nc_m <- getCoverage(BSseq = nc_obj, type = "M",
                            what = "perBase")
        
        nc_rate <- sum(nc_m) / sum(nc_cov)
        
        # Subtract non-conversion rate with a floor of zero
        gr$pc <- gr$pc - nc_rate
        gr$pc[gr$pc < 0] <- 0
        
        return(gr)
    }
    
    if (nc_correct == TRUE){
        message("Correcting for non-conversion...")
        gr_mc <- correct_nc(gr_mc)
    }
    
    #---- Extract last range to resize seperately or else it will overlap
    message("Preparing data for bigwig output...")
    
    ind <- seqnames(gr_mc)
    ends <- ind@lengths %>% cumsum()
    gr_ends <- gr_mc[ends] 
    gr_mc <- gr_mc[-ends]
    gr_mc <- GenomicRanges::resize(gr_mc, width = win_step, fix="center")
    
    ind2 <- seqnames(gr_mc)
    ends2 <- ind2@lengths %>% cumsum()
    max_val <- end(gr_mc[ends2])
    start(gr_ends) <- max_val + 1
    gr_mc <- c(gr_mc, gr_ends) %>% sort()
    
    gr_mc <- keepSeqlevels(x = gr_mc, value = contigs)
    
    all_seqlengths <- seqlengths(x = BSgenome.Hsapiens.UCSC.hg19)
    all_seqlengths <- all_seqlengths[names(all_seqlengths) %in% seqnames(gr_mc)]
    seqlengths(gr_mc) <- all_seqlengths
    
    gr_mc <- GenomicRanges::trim(gr_mc)
    strand(gr_mc) <- "+"

    gr_mc$score <- as.numeric(gr_mc$pc)
    gr_mc$pc <- NULL
    
    # Write the output file
    message("Writing output file...")
    export.bw(object = gr_mc, con = out_path)
    
    message("Done!")
}

CGmap_to_window_bigwig(path = cgmap_path)

fls <- list.files(path = ".", pattern = ".CGmap.gz$", full.names = TRUE)
fls <- fls[!grepl(x=fls, pattern="ATCGmap.gz")]

lapply(X = fls, FUN = CGmap_to_window_bigwig)

