#!/bin/bash
#rj name=dmrseq queue=hg19ips features=knl,centos7,fastio logdir=logs schema=input.schema runtime=24
set -euo pipefail

# Add the modules needed for the R script
module add hpc
module use /p9/mcc_hg19ips/sw/modules
module add bashutils
module add R

logvars manifest outpath prefix

which Rscript

Rscript ./run_dmrseq_dug_full.R "$manifest" "$outpath" "$prefix"
