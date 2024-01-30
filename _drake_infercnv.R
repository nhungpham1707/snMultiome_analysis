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
filename <- '25012024_all_multiome_lib.csv'
metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
specialLib <- c("LX093_LX094_an_163")
specialLibInd <- grep(specialLib, metadata$name)
nospecialMet <- metadata[-specialLibInd,]
# get libraries that only demultiplex 
# with souporcell
soclId <- specialLib %>% map(splitName)
# get multiplex libraries list 
mulLib <- unique(nospecialMet$name[nchar(nospecialMet$souporcell_link) > 0])
# mulLib <- mulLib[1]
mulId <- mulLib %>% map(splitName)
# get single libraries list 
sngLib <- unique(nospecialMet$name[nchar(nospecialMet$souporcell_link) == 0])
# sngLib <- sngLib[1]
snglId <- sngLib %>% map(splitName)
# get all samples to make combine peaks

lbLst <- unique(c(specialLib, mulLib, sngLib)) 
idLst <- lbLst %>% map(splitName)
# add merge sample
lbLst <- c(lbLst, 'merge')
idLst <- c(idLst, 'merge')
# define normal cells 
normal_cells <- c('B_cell', 'T_cells', 'Macrophage', 'Monocyte', 'NK_cell')
## define plan ----
inLink <- AtacInferInputDir
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
make(infercnv_plan,lock_envir = TRUE, lock_cache = FALSE, verbose = 0)

# options(clustermq.scheduler = "multicore") # nolint
# make(infercnv_plan, parallelism = "clustermq", jobs = 2, lock_cache = FALSE)

vis_drake_graph(infercnv_plan, targets_only = TRUE, lock_cache = FALSE, file = 'infercnv.png', font_size = 20 )