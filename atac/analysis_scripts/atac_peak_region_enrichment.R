source("R/server_libraries_and_functions.R")

all_elements_gr <- readRDS("resources/all_genomic_elements_granges_final.Rds")

#### Peak enrichment
test_atac_enrichment <- function(target_element="Fibroblast_H3K9me3", test_bed,
                                 permutations=200, cores=4, alt="auto"){

    stopifnot(file.exists(test_bed))
    gc()
    message(Sys.time())+
        message(str_c("Running ", target_element, "..."))
    target_gr <- all_elements_gr[all_elements_gr$class == target_element]
    message(str_c("n = "), length(target_gr))

    test_gr <- bed_to_gr(test_bed)

    test_gr <- resize(test_gr, width = 100, fix = "center")

    set.seed(123)
    perm_test_results <- regioneR::permTest(A=test_gr, B=target_gr,
                                            alternative = "auto",
                                            #alternative = alt,
                                            randomize.function=regioneR::randomizeRegions,
                                            ntimes = permutations,
                                            evaluate.function=regioneR::numOverlaps,
                                            count.once = TRUE,
                                            genome="hg19",
                                            #mask = unmappable_gr,
                                            allow.overlaps=FALSE,
                                            per.chromosome=TRUE,
                                            mc.cores=cores,
                                            mc.set.seed=FALSE, force.parallel=TRUE)

    lz <- regioneR::localZScore(pt=perm_test_results, A=test_gr, B=target_gr,
                      step = mean(width(test_gr)/2),
                      window = 100000)

    results <- list(pt=perm_test_results, lz=lz, element=target_element,
                    test_bed=test_bed)

    message("Done!")
    return(results)
}

primed_down_bed <- "atac/all_differential_peaks_primed_down.bed"

test_elements <- c("Fibroblast_LAD", "ESC_LAD", "Constitutive_LAD", "Fibroblast_H3K9me3",
  "ESC_H3K9me3", "Constitutive_H3K9me3",
  "Fibroblast_PMD")

#down_fib_h3k9 <- test_atac_enrichment(target_element = "Fibroblast_H3K9me3",
#                                      test_bed = primed_down_bed)

down_primed_enrich <- lapply(test_elements, test_atac_enrichment,
                             test_bed = primed_down_bed, cores=10)
saveRDS(down_primed_enrich, "atac/primed_down_peak_region_enrichments.Rds")

