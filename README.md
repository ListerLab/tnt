# Epigenetic and functional correction of human iPS cells by transient naive reprogramming

Analysis code for Buckberry, Liu, Poppe, Tan et al. "Epigenetic and functional correction of human iPS cells by transient naive reprogramming"

Sequence data for this study are accessible under GEO SuperSeries GSE159297.

Details of analyses and figures from this study are in the markdown HTML files in the top-level directory of this repository. Additional information on data processing can be found in data type sub-directories `/wgbs/analysis_scripts`, `RNAseq/analysis_scripts`, `/ChIPseq/`, and `/ATACseq/analysis_scrpits`.   

A vast majority of analyses and figures presented in this study were produced with processed WGBS data in CGmap files that can be downloaded from GEO. We used these CGmap files to create [BSseq objects](https://www.bioconductor.org/packages/devel/bioc/vignettes/bsseq/inst/doc/bsseq.html#3_Using_objects_of_class_BSseq) for data analysis in R.  
For each sample, we generated two BSseq objects: One for CpG methylation, and another for CpA methylation using the scripts `wgbs/analysis_scripts/make_BSseq_CG_context_obj.R` `wgbs/analysis_scripts/make_BSseq_CA_context_obj.R`.    
Downloading the CGmap files for WGBS data and processing with these scripts should generate `.Rds` files containing BSseq objects with the same naming convention in the sample metadata file `wgbs/metadata/wgbs_metadata_local.csv`.  
This metadata file has the information required to reproduce the analyses detailed in the R markdown files in this repository.  










