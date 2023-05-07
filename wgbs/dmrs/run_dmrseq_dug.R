library(bsseq)
library(dmrseq)
library(magrittr)
library(stringr)
library(BiocParallel)
library(batchtools)

args = commandArgs(TRUE)

# Path to DMR manifest file. This is a csv file with the column names id, group, rds_path
manifest <- args[1] # must be absolute path
out_path <- args[2] # must be folder path ending in /
prefix <- args[3] # Out file prefix

#-------- Read the DMR manifest file

# Check that the manifest file exists
message("=== Checking manifest exists ===")
stopifnot(file.exists(manifest) == TRUE)

# Read the manifest file
man_dat <- read.table(manifest, header = TRUE, sep = ",",
                      colClasses = "character")

# Check that Rds files exist
message("=== Checking files in manifest exist ===")
rds_check <- lapply(man_dat$rds_path, file.exists) %>% unlist()
stopifnot(all(rds_check) == TRUE)

# Get groups and test that only two groups exist
message("=== Checking there are only two groups in manifest ===")
uniq_groups <- man_dat$group %>% unique()
stopifnot(length(uniq_groups) == 2)

# Check that library ids are unique
message("=== Checking ids are unique ===")
stopifnot(length(unique(man_dat$id)) == nrow(man_dat))

# Chromosomes to test for DMRs
message("=== Setting chromosomes to test for DMRs ===")
chrom_list <- str_c("chr", 1:22)
message(chrom_list)

message("=== Reading BSseq objects ===")
# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path, chroms){

    bs_obj <- readRDS(file = rds_path)
    bs_obj <- chrSelectBSseq(BSseq = bs_obj,
                             seqnames = chroms,
                             order = TRUE)
    bs_obj <- strandCollapse(bs_obj)
    return(bs_obj)

}

# Load the data, and sub-select targeted chromosomes
obj_list <- lapply(X = man_dat$rds_path, read_bs_obj,
                   chroms = chrom_list)

# Combine all of the bsseq objects into one
message("=== Combining BSseq objects ===")
obj_list <- bsseq::combineList(x = obj_list)

# Remove CpG with no coverage
message("=== Removing CpG with no coverage ===")
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

# Setup the groups
pData(obj_list)$CellType <- factor(man_dat$group)
colnames(obj_list) <- man_dat$id

#-------- Call DMRs
BiocParallel::register(BPPARAM = MulticoreParam(workers = 14))

message("=== Calling DMRs ===")
dmrs <- dmrseq(obj_list, testCovariate = "CellType",
               bpSpan = 500,
               maxGap = 500,
               maxPerms = 20)

# Output results
message("=== Formatting results ===")
out_list <- list(manifest_data = man_dat, dmr_granges = dmrs)

out_file <- str_c(out_path, prefix, "_", uniq_groups[1],
                  "_vs_", uniq_groups[2],
                  "_dmrseq_dmrs.Rds")

message("=== Saving output file ===")
saveRDS(object = out_list, file = out_file)

message(str_c("DMR file saved to ", out_file))
