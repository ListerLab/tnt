
installed <- installed.packages()

# List CRAN packages that must be installed and then check if already installed
cran_packages <- c("stringr", "magrittr", "data.table")

cran_packages <- cran_packages[!cran_packages %in% installed[ ,1]]

# Install the required packages
if (length(cran_packages) >= 1) {
  lapply(cran_packages, install.packages, repos="https://cloud.r-project.org")
}

# List bioc packages that must be installed and then check if already installed
bioc_packages <- c("bsseq")

bioc_packages <- bioc_packages[!bioc_packages %in% installed[ ,1]]

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (length(bioc_packages) >=1)
    BiocManager::install(bioc_packages)

library(bsseq)
library(magrittr)
library(stringr)
library(data.table)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)


# Path to a CGmap file from BSseeker2
path <- args[1]

make_bsseq_obj <- function(CGmap_path, context="CG"){
        
        gc()
        
        id <- basename(CGmap_path) %>% stringr::str_replace(pattern = ".CGmap.gz",
                                                            replacement = "")
            
        message(str_c("Loading ", CGmap_path, "..."))
        
        dat <- data.table::fread(CGmap_path, header = FALSE,
                                 select = c(1,2,3,5,7,8), sep = "\t",
                                 col.names = c("chr", "base", "position",
                                               "diContext", "C_reads", "CT_reads"))
        
        # subset context
        dat <- dat[dat$diContext == context, ]
        dim(dat)
        # Add strand
        dat$base <- ifelse(test = dat$base == "C", yes = "+", no = "-")
        
        #Load data into BSseq object
        message(str_c("Making BSseq object for ", id, "..."))
        
        bs_obj <- bsseq::BSseq(chr = dat$chr, pos=dat$position,
                        sampleNames = id,
                        M = matrix(dat$C_reads),
                        Cov = matrix(dat$CT_reads))
        
        # Add strand
        bsseq::strand(bs_obj) <- dat$base
        
        # Collapse strand for CG
        if (context=="CG"){
            bs_obj <- bsseq::strandCollapse(bs_obj)
        }

        message("Saving BSseq object...")
        saveRDS(object = bs_obj, file = str_c(id, "_", context, "_bsseq_obj.Rds"))
        message("Done!")
}

make_bsseq_obj(CGmap_path=path, context="CG")

