#!/usr/bin/env Rscript
setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
library(infercnv)
library(drake)
library(dplyr)
library(purrr)
## Load your R files ----
source('./R/global_variables_infercnv.R')
source('./R/infercnv_function.R')
source('./R/utils.R')

# define normal cells 
normal_cells <- c('B_cell', 'T_cells', 'Macrophage', 'Monocyte', 'NK_cell')
## define plan ----
inLink <- AtacnoRInferInputDir
outLink <- cellAtacnoRInferDir
filename <- '30012024_remove_relapse.csv'
metadata <- getData(filename, delim = ',')

lib <- c(unique(metadata$name), 'merge')
geneOderLink <- paste0(base_data_dir, '/gencode_v19_gene_pos.txt')
mrginfercnv_plan <- drake_plan(
    AtacIfcnvOb= make_infercnvObj(lib, normal_cells, geneOderLink, inLink),
    AtacIfcnvRes = run_infercnv(AtacIfcnvOb, outLink)
)

vis_drake_graph(mrginfercnv_plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_mrginfercnv.png', font_size = 20 )