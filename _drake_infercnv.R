#!/usr/bin/env Rscript
setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
library(infercnv)
library(drake)
library(dplyr)
library(purrr)
library(ggplot2)
library(futile.logger)
## Load your R files ----
source('./global_variables_infercnv.R')
source('./R/infercnv_function.R')
source('./R/utils.R')

## read metadata ----
# filename <- '30_11_2023_all_multiome_libraries.csv' 
# filename <- '25012024_all_multiome_lib.csv'
filename <- '30012024_remove_relapse.csv'

metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
specialLib <- c("LX189_LX190_an_325", "LX099_LX100_an_166") # fail samples

# specialLib <- c("LX189_LX190_an_325", "LX099_LX100_an_166", 'LX051_LX052_an_128') # fail samples w new normal cells
nospecialMet <- metadata[!metadata$name %in% specialLib,]
lbLst <- unique(nospecialMet$name) 
idLst <- lbLst %>% map(splitName)
# lbLst <- c(lbLst, 'merge')
# idLst <- c(idLst, 'merge')

# define normal cells 
normal_cells <- c('B_cell', 'T_cells', 'Macrophage', 'Monocyte', 'NK_cell')

# normal_cells_HM <- c('B_cell', 'Macrophage', 'Monocyte', 'NK_cell') # remove T cell as there is only 1 cell in LX049, cannot run HMM 
# normal_cells <- c('B_cell', 'BM & Prog.', 'DC', 'CMP', 'Monocyte', 'Neutrophils', 'NK_cell', 'Pro-B_cell_CD34+', 'T_cells', 'Endothelial_cells', 'Erythroblast', 'Gametocytes', 'GMP', 'HSC_-G-CSF', 'HSC_CD34+', 'Macrophage', 'MEP', 'MSC', 'Myelocyte', 'Osteoblasts', 'Platelets', 'Pre-B_cell_CD34−', 'Pro-Myelocyte')
## define plan ----
inLink <- AtacInferInputDir
outLink <- cellAtacInferDir 
rnaInLink <- rnaInferInputDir
rnaOutLink <- cellRnaIcnvdir 
geneOrderLink <- paste0(base_data_dir, '/gencode_v19_gene_pos.txt')

atac_infercnv_plan <- drake_plan(
    AtacIfcnvOb= target(make_infercnvObj(lib, normal_cells, geneOrderLink, inLink),
                transform = map(lib = !!lbLst,
                            id.vars = !!idLst,
                            .id = id.vars)),
    AtacnIfcnvRes = target(run_infercnv(AtacIfcnvOb, outLink),
                transform = map(AtacIfcnvOb,
                            id.vars = !!idLst,
                            .id = id.vars))
)


rna_infercnv_plan <- drake_plan(
    rnaIfcnvOb= target(make_infercnvObj(lib, normal_cells, geneOrderLink, rnaInLink),
                transform = map(lib = !!lbLst,
                            id.vars = !!idLst,
                            .id = id.vars)),
    rnaIfcnvRes = target(run_infercnv(rnaIfcnvOb, rnaOutLink),
                transform = map(rnaIfcnvOb,
                            id.vars = !!idLst,
                            .id = id.vars))
)

plan <- bind_plans(atac_infercnv_plan, rna_infercnv_plan)

make(plan,lock_envir = TRUE, lock_cache = FALSE, verbose = 0)
vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'infercnv_before.png', font_size = 20 )
# options(clustermq.scheduler = "multicore") # nolint
# make(infercnv_plan, parallelism = "clustermq", jobs = 2, lock_cache = FALSE)

vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_infercnv_after.png', font_size = 20 )