#### plot mofif enrichments
library(XML)

source("R/project_functions.R")

#### Plot the heatmap for enrichment of known factors
homer_known <- list.files("atac/",
                          pattern = "knownResults.html", recursive = TRUE,
                          full.names = TRUE)

primed_down <- "atac/all_differential_peaks_primed_downhomer/knownResults.html"
primed_up <- "atac/all_differential_peaks_primed_uphomer/knownResults.html"
file.exists(primed_up)

#html <- homer_known[1]

read_homer_known <- function(html, p_threshold=0.01, fdr_threshold=0.01,
                             fg_bg_cut=1.2){

    dat <- readHTMLTable(html, header = TRUE, as.data.frame = TRUE, stringsAsFactors = FALSE)[[1]]

    ## Clean up the table
    dat <- dat[ ,c(3:10)]

    colnames(dat) <- c("Name", "P", "logP", "q", "nTarget", "pcTarget", "nBackground", "pcBackground")
    tf<- str_split(dat$Name, pattern = "/", simplify = TRUE)[ ,1]

    # Get the stats
    pc_target <- str_replace(string = dat$pcTarget,
                             pattern = "%",
                             replacement = "") %>% as.numeric()

    pc_background <- str_replace(string = dat$pcBackground,
                                 pattern = "%",
                                 replacement = "") %>% as.numeric()

    pval <- as.character(dat$logP) %>% as.numeric()
    pval <- 10^pval

    fdr <- as.numeric(dat$q)

    # Get the cluster ID
    cluster <- dirname(html) %>% str_split(pattern = "/", simplify = TRUE)
    cluster <- cluster[ ,ncol(cluster)]

    clean_dat <- data.frame(TF=tf, pc_target=pc_target,
                            pc_background=pc_background,
                            pval=pval, fdr=fdr,
                            cluster=cluster)

    ## Remove duplicate TFs
    clean_dat <- clean_dat[!duplicated(clean_dat$TF), ]

    ## Keep only those at present in more than 10% of sequences
    clean_dat <- clean_dat[clean_dat$pc_target > 10, ]

    ## Filter on FG:BG
    clean_dat$fg_bg <- clean_dat$pc_target / clean_dat$pc_background
    clean_dat <- clean_dat[clean_dat$fg_bg >= fg_bg_cut, ]

    # Filter on P-value
    clean_dat <- clean_dat[clean_dat$fdr <= fdr_threshold, ]
    clean_dat <- clean_dat[clean_dat$pval <= p_threshold, ]

    return(clean_dat)
}

results <- lapply(c(primed_down, primed_up), read_homer_known)

results <- do.call(rbind, results)

# Change name for reformating
res_df <- results

# Convert the pvalues
res_df$pval <- -log10(res_df$pval)
res_df$pval[is.infinite(res_df$pval)] <- max(res_df$pval[!is.infinite(res_df$pval)]) + 1


# Add TF id
tf_ids <- as.character(res_df$TF) %>%
    str_split(pattern = stringr::fixed("("), simplify = TRUE) %>% data.frame()
res_df$tf_id <- tf_ids$X1
res_df$tf_class <- str_replace(string = tf_ids$X2, pattern = stringr::fixed(")"), replacement = "")

res_df$cluster %<>% str_remove_all(pattern = "all_differential_peaks_") %>%
    str_remove_all("homer")

gg_motif <- ggplot(res_df, aes(y = factor(tf_id),
                               x = factor(cluster), group=tf_class)) +        ## global aes
    #geom_tile(aes(fill = pval)) +         ## to get the rect filled
    geom_point(aes(fill = pval, size = pc_target),
               colour="black", pch=21, stroke=0.176389)  +    ## geom_point for circle illusion
    scale_fill_gradient(low = "grey", high = "red") +       ## color of the corresponding aes
    scale_size(range = c(1, 3))+             ## to tune the size of circles
    xlab("") +
    ylab("") +
    facet_grid(tf_class~., scales = "free", space = "free") +
    #scale_y_discrete(label="", limits=levels(res_df$tf_id)) +
    sams_pub_theme(legend_pos = "right")
gg_motif


pdf("atac/plots/atac_primed_differential_motif_plots.pdf", width = 3, height = 5)
gg_motif
dev.off()


wb_ed_fig8c <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb_ed_fig8c, sheetName = "ED_Fig_8c")
openxlsx::writeData(wb = wb_ed_fig8c, sheet = "ED_Fig_8c",
                    x = gg_motif$data)
openxlsx::saveWorkbook(wb = wb_ed_fig8c,
                       file = "ED_Figure_8c_source_data.xlsx", overwrite = TRUE)

