#!/usr/bin/env Rscript
setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code')
## Load your packages, e.g. library(drake).
library(infercnv)
library(drake)
library(dplyr)
library(purrr)
## Load your R files ----
source('./R/global_variables_infercnv.R')
source('./R/infercnv_function.R')
source('./R/utils.R')

## read metadata ----
filename <- '30_11_2023_all_multiome_libraries.csv' 
metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
metadata <- metadata[c(1:3,5:6),] # to test code 
lbLst <- unique(metadata$name)
idLst <- lbLst %>% map(splitName)
normal_cells <- c('B_cell', 'T_cells', 'Macrophage', 'Monocyte', 'NK_cell')
## define plan ----
inLink <- MergAtacInferInputDir
outLink <- cellAtacInferDir
geneOderLink <- paste0(base_data_dir, '/gencode_v19_gene_pos.txt')
infercnv_plan <- drake_plan(
    AtacIfcnvOb= target(make_infercnvObj(lib, normal_cells, geneOderLink, inLink),
                transform = map(lib = !!lbLst,
                            id.vars = !!idLst,
                            .id = id.vars)),
    AtacIfcnvRes = target(run_infercnv(AtacIfcnvOb, outLink),
                transform = map(AtacIfcnvOb,
                            id.vars = !!idLst,
                            .id = id.vars))
)

vis_drake_graph(infercnv_plan, targets_only = TRUE, lock_cache = FALSE, file = 'infercnv.png', font_size = 20 )