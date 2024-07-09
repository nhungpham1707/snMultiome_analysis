setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

  logistic_atac <- drake_plan(
# atac ---
    dsc_atac = readRDS('output/healthy_data/dsc_atac_hthymrgDim.RDS'),
    train_dsc_atac40k = trainModel(GetAssayData(dsc_atac), classes = dsc_atac$cell_type, maxCells = 40000)

  )

  make(logistic_atac, lock_cache = FALSE,  lock_envir = FALSE)