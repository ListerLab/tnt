# Epigenetic and functional correction of human iPS cells by transient naive reprogramming

Analysis code for Buckberry, S., Liu, X., Poppe, D. et al. Transient naive reprogramming corrects hiPS cells functionally and epigenetically. Nature 620, 863–872 (2023). https://doi.org/10.1038/s41586-023-06424-7

Raw and processed data for this study are accessible under [GEO SuperSeries GSE159297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159297).

### Markdown generated HTML files with analysis code for figures
- [Figure 1 and associated Extended Data Figures](Figure_1.md)
- [Figure 2 and associated Extended Data Figures](Fig_2.md)
- [Figure 3 and associated Extended Data Figures](Fig_3.md)
- [Extended Data Figure 2b: CNV analysis](ED_Fig_2b_CNV_analysis.md)
- [Extended Data Figure 2g,h,i: Ma et al. re-analysis](ED_Fig_3ghi_Ma_et_al_analysis.md)
- [Extended Data Figure 3f,g,h: Primed-Naive-Primed analysis](ED_Fig_4fgh_PNP_iPSC_analysis.md)
- [Extended Data Figure 3i: Lentiviral insertion analysis](ED_Fig_4i_lenti_insertion_analysis.md)
- [Figure 4 and associated Extended Data Figures](Fig_4.md)
- [Extended Data Figure 8k: Choi et al. reanalysis](Choi_ESC_iPSC_differential_expression.md)
- [Extended Data Figure 8l: Ma et al. reanalysis](SCNT_differential_expression.md)
- [Extended Data Figure 9c,d,e: Evaluation of published criteria](REVISION_Koyanagi_Ruiz_genes_TE_methylation.md)
- [Figure 5 and associated Extended Data Figures](REVISION_differentiation_quantifications.md)  
- [CG-DMR analyses](CG_DMR_analysis.md)
- [CH-DMR analyses](CH_DMR_analysis.md)
- [Imprinting analyses](REVISION_imprinting_analyses.md)

### Installing project packages and loading custom functions
Run the following R scripts to install and load R packages and functions required to run the analyses detailed in this repo. 
`R/install-project-packages.R`  
`R/project_functions.R`

### How to make the processed data files for analyses
Details of analyses and figures from this study are in the markdown files listed above. Additional information on data processing can be found in data type sub-directories `/wgbs/analysis_scripts`, `RNAseq/analysis_scripts`, `/ChIPseq/`, and `/ATACseq/analysis_scrpits`.   
A vast majority of analyses and figures presented in this study were produced with processed WGBS data in CGmap file format that can be downloaded from GEO. These CGmap files can be used to create [BSseq objects](https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html#3_Using_objects_of_class_BSseq) for data analysis in R.   

For each sample, we generated two BSseq objects: One for CpG methylation, and another for CpA methylation using the scripts `wgbs/analysis_scripts/make_BSseq_CG_context_obj.R` `wgbs/analysis_scripts/make_BSseq_CA_context_obj.R`.  

Downloading the CGmap files for WGBS data and processing with these scripts should generate `.Rds` files containing BSseq objects with the same naming convention in the sample metadata file `wgbs/metadata/wgbs_metadata_local.csv`.  

This metadata file has the information required to reproduce the analyses detailed in the R markdown files in this repository.  

All additional files for running the analyses in these Rmarkdown documents will be availble if you clone this repository. 
