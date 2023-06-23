# Epigenetic and functional correction of human iPS cells by transient naive reprogramming

Analysis code for Buckberry, Liu, Poppe, Tan et al. "Epigenetic and functional correction of human iPS cells by transient naive reprogramming"

Sequence data for this study are accessible under GEO SuperSeries GSE159297.

### Markdown generated HTML files with analysis code for figures
- [Figure 1 and associated Extended Data Figures](Figure_1.html)
- [Figure 2 and associated Extended Data Figures](Fig_2.html)
- [Figure 3 and associated Extended Data Figures](Fig_3.html)
    - [Extended Data Figure 2b: CNV analysis](ED_Fig_2b.html)
    - [Extended Data Figure 2g,h,i: Ma et al. re-analysis]()
    - [Extended Data Figure 3f,g,h: Primed-Naive-Primed analysis]()
    - [Extended Data Figure 3i: Lentiviral insertion analysis]()
    - [Extended Data Figure 4: DMR reproducability analysis]()
- [Figure 4 and associated Extended Data Figures]()
    - [Extended Data Figure 8k: Choi et al. reanalysis](Choi_ESC_iPSC_differential_expression.html)
    - [Extended Data Figure 8l: Ma et al. reanalysis](SCNT_differential_expression.html)
- [Figure 5 and associated Extended Data Figures]()
- [CG-DMR analyses](CG_DMR_analysis.html)
- [CH-DMR analyses](CH_DMR_analysis.nb.html)
- [Imprinting analyses](REVISION_imprinting_analyses.html)

### How to make the processed data files for analyses
Details of analyses and figures from this study are in the markdown HTML files listed above. Additional information on data processing can be found in data type sub-directories `/wgbs/analysis_scripts`, `RNAseq/analysis_scripts`, `/ChIPseq/`, and `/ATACseq/analysis_scrpits`.   

A vast majority of analyses and figures presented in this study were produced with processed WGBS data in CGmap files that can be downloaded from GEO. These CGmap files can be useed to create [BSseq objects](https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html#3_Using_objects_of_class_BSseq) for data analysis in R.   

For each sample, we generated two BSseq objects: One for CpG methylation, and another for CpA methylation using the scripts `wgbs/analysis_scripts/make_BSseq_CG_context_obj.R` `wgbs/analysis_scripts/make_BSseq_CA_context_obj.R`.  

Downloading the CGmap files for WGBS data and processing with these scripts should generate `.Rds` files containing BSseq objects with the same naming convention in the sample metadata file `wgbs/metadata/wgbs_metadata_local.csv`.  

This metadata file has the information required to reproduce the analyses detailed in the R markdown files in this repository.  











