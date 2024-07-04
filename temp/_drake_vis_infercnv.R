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

## read metadata ----
# filename <- '30_11_2023_all_multiome_libraries.csv' 
# filename <- '25012024_all_multiome_lib.csv'
filename <- '30012024_remove_relapse.csv'

metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
specialLib <- c("LX189_LX190_an_325") # fail sample
specialLibInd <- grep(specialLib, metadata$name)
nospecialMet <- metadata[-specialLibInd,]
lbLst <- unique(nospecialMet$name) 
idLst <- lbLst %>% map(splitName)
lbLst <- c(lbLst, 'merge')
idLst <- c(idLst, 'merge')
# define normal cells 
normal_cells <- c('B_cell', 'T_cells', 'Macrophage', 'Monocyte', 'NK_cell')
## define plan ----
inLink <- AtacnoRInferInputDir
outLink <- cellAtacnoRInferDir 
geneOderLink <- paste0(base_data_dir, '/gencode_v19_gene_pos.txt')
infercnv_plan <- drake_plan(
    AtacnoRIfcnvOb= target(make_infercnvObj(lib, normal_cells, geneOderLink, inLink),
                transform = map(lib = !!lbLst,
                            id.vars = !!idLst,
                            .id = id.vars)),
    AtacnoRIfcnvRes = target(run_infercnv(AtacnoRIfcnvOb, outLink),
                transform = map(AtacnoRIfcnvOb,
                            id.vars = !!idLst,
                            .id = id.vars))
)

vis_drake_graph(infercnv_plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_infercnv.png', font_size = 20 )