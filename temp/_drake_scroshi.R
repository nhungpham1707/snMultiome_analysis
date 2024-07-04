setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

loadd(mrgGA)
saveRDS(mrgGA, paste0(atcMrgDir, '/mrgGA.RDS'))
hpca.se <- celldex::HumanPrimaryCellAtlasData()
labels = hpca.se$label.main

  scroshi_plan <- drake_plan(
  # run singleR ---
  mrgGA_input = as.matrix(mrgGA[['RNA']]@data),  
  mrgatacGA_singler = SingleR(mrgGA_input,
                ref=hpca.se,
                labels=labels,
                method = NULL,
                clusters = mrgGA$seurat_clusters,
                genes = "de",
                sd.thresh = 1,
                de.method = "classic",
                de.n = NULL,
                de.args = list(),
                aggr.ref = FALSE,
                aggr.args = list(),
                recompute = TRUE,
                restrict = NULL,
                quantile = 0.8,
                fine.tune = TRUE,
                tune.thresh = 0.05,
                prune = TRUE,
                assay.type.test = "logcounts",
                assay.type.ref = "logcounts",
                check.missing = TRUE),
  saveSingleRMrg = saveRDS(mrgatacGA_singler, paste0(atacCellSngRDir, '/sR_mrg_all_cluster.RDS')),
  mrg_scroshi_demo = run_scROSHI_mrg_w_demo_data(sr = mrgGA, cols = my_cols, pt = 1, save_name = 'mrg_all_demo', path_to_save = atacScroshiDir),
  
  mrg_scroshi_atrt = run_scROSHI_mrg_w_atrt_data(sr = mrgGA, cols = my_cols, pt = 1, save_name = 'mrg_all_atrt', path_to_save = atacScroshiDir)
  )

  make(scroshi_plan, lock_envir = FALSE, lock_cache = FALSE)