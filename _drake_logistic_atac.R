setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

loadd(new_atachm_mx)
loadd(atac_hm_tumor_nona)
colnames(new_atachm_mx) <- colnames(atac_hm_tumor_nona)
  logistic_atac <- drake_plan(
# atac ---
    dsc_atac = readRDS('output/healthy_data/dsc_atac_hthymrgDim.RDS'),
    # train_dsc_atac40k = trainModel(GetAssayData(dsc_atac), classes = dsc_atac$cell_type, maxCells = 40000),
    # p_dsc_atac = predictSimilarity(train_dsc_atac40k, new_atachm_mx, 
    #                         classes = atac_hm_tumor_nona$cell_identity,
    #                         minGeneMatch = 0.0,
    #                         logits = F)
    dsc_only_overlap = subset(dsc_atac, features = rownames(new_atachm_mx)),
    train_dsc_atac_onlyoverlap = trainModel(GetAssayData(dsc_only_overlap), classes = sub_dsc$cell_type, maxCells = 40000),
    p_dsc_atac = predictSimilarity(train_dsc_atac_onlyoverlap, new_atachm_mx, 
                            classes = atac_hm_tumor_nona$cell_identity,
                            minGeneMatch = 0.7,
                            logits = F)
  )

  make(logistic_atac, lock_cache = FALSE,  lock_envir = FALSE)