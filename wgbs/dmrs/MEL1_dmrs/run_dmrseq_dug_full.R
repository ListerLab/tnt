.libPaths('/p9/mcc_hg19ips/sw/R/4.0.2-gcc-blas/library')

message("==== Loading libraries ====")
library(bsseq)
library(dmrseq)
library(magrittr)
library(stringr)
library(BiocParallel)
message("==== Finished loading libraries ====")

args = commandArgs(TRUE)

# Path to DMR manifest file. This is a csv file with the column names id, group, rds_path 
manifest <- args[1] # must be absolute path
out_path <- args[2] # must be folder path ending in /
prefix <- args[3] # Out file prefix

#-------- Read the DMR manifest file

message("==== Check that the manifest file exists ====")
# Check that the manifest file exists
stopifnot(file.exists(manifest) == TRUE)

message("==== Read the manifest file ====")
# Read the manifest file
man_dat <- read.table(manifest, header = TRUE, sep = ",",
                      colClasses = "character")

message("==== Check that Rds files exist ====")
# Check that Rds files exist
rds_check <- lapply(man_dat$rds_path, file.exists) %>% unlist()
stopifnot(all(rds_check) == TRUE)

message("==== Get groups and test that only two groups exist ====")
# Get groups and test that only two groups exist
uniq_groups <- man_dat$group %>% unique()
stopifnot(length(uniq_groups) == 2)

message("==== Check that library ids are unique ====")
# Check that library ids are unique
stopifnot(length(unique(man_dat$id)) == nrow(man_dat))

# Chromosomes to test for DMRs
chrom_list <- str_c("chr", c(1:22, "X", "Y"))

message("==== Read a Bs_seq boject from .Rds file ====")
# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path, chroms){
    
   bs_obj <- readRDS(file = rds_path)
   bs_obj <- chrSelectBSseq(BSseq = bs_obj,
                            seqnames = chroms,
                            order = TRUE)
   bs_obj <- strandCollapse(bs_obj)
   return(bs_obj)

}

message("==== Load the data, and sub-select targeted chromosomes ====")
# Load the data, and sub-select targeted chromosomes
obj_list <- lapply(X = man_dat$rds_path, read_bs_obj,
                   chroms = chrom_list)

message("==== Combine all of the bsseq objects into one ====")
# Combine all of the bsseq objects into one
obj_list <- bsseq::combineList(x = obj_list)

message("==== Remove CpG with no coverage ====")
# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

message("==== Setup the groups ====")
# Setup the groups
pData(obj_list)$CellType <- factor(man_dat$group)
colnames(obj_list) <- man_dat$id

message("==== setting up BiocParallel ====")

# The line below will need to be removed to test the above 
BiocParallel::register(BPPARAM = MulticoreParam(workers = 30))

registered('MulticoreParam')

chunks <- 1

message(str_c("Chunks = ",chunks))

message("==== Running dmrseq ====")
dmrs <- dmrseq(obj_list, testCovariate = "CellType",
               maxPerms = 10,
               bpSpan = 500,
               maxGap = 500,
               chrsPerChunk = chunks)

message("==== Output results ====")
# Output results
out_list <- list(manifest_data = man_dat, dmr_granges = dmrs)

out_file <- str_c(out_path, prefix, "_", uniq_groups[1],
                 "_vs_", uniq_groups[2],
                 "_dmrseq_dmrs.Rds")

saveRDS(object = out_list, file = out_file)

message(str_c("DMR file saved to ", out_file))

