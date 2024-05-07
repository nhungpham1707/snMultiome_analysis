# install if needed 
if(!require(dotenv, quietly = TRUE)){
  install.packages('dotenv')
}

if(!require(devtools)){
  install.packages("devtools")
}
if(!require('clustermq', quietly = TRUE)){
  install.packages("clustermq")
}
if(!require(SCutils, quietly = TRUE)){
  devtools::install_bitbucket("princessmaximacenter/scutils")
}

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# test why it still install scater
# if (!require("scater", quietly = TRUE)){
# BiocManager::install("scater")
# }

if (!require("scROSHI", quietly = TRUE)){
  devtools::install_github("ETH-NEXUS/scROSHI")
} # may not work, if fail try 
# devtools::install_github("ETH-NEXUS/scROSHI")
if (!require("clustermq", quietly = TRUE)){
  install.packages("clustermq")
}

suppressPackageStartupMessages({
  # packages for drake pipeline management
  library(conflicted)
  library(dotenv)
  library(drake)
  # library(dflow)
  
  # packages for the analysis
  library(Signac)  
  library(Seurat)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(tidyverse)
  library(celldex)
  library(SingleR)
  library(regioneR)
  library(SCutils) # from maxima bitbucket for genelist
  library(scater) 
  library(reshape2)
  library(ggplot2)
  library(SingleR)
  library(celldex)
  library(scROSHI) # dev github version  
  library(scales) # for default colors 
  library(tools)
  

  library(SingleCellExperiment)
  library(cowplot)
  library(limma)
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(CellMixS)
  library(SeuratDisk)
  library(SCpubr)
  conflict_prefer('select', 'clusterProfiler')
  conflict_prefer('filter', 'clusterProfiler')
  
  conflict_prefer('intersect', 'base')
  conflict_prefer('setdiff', 'base')
  conflict_prefer('union', 'base')

})
  conflicts_prefer(base::unname)

suppressPackageStartupMessages({
  library(harmony)
  library(SingleCellExperiment)
  library(cowplot)
  library(limma)
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(scater)
  library(CellMixS)
  library(bluster) # for clustering behavior 
  library(cowplot)
  library(lisi)
  library(EnsDb.Hsapiens.v86)
})