setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

loadd(mrgGA)
  scroshi_plan <- drake_plan(
  mrg_scroshi_demo = run_scROSHI_mrg_w_demo_data(sr = mrgGA, cols = my_cols, pt = 1, save_name = 'mrg_all_demo', path_to_save = atacScroshiDir),
  
  mrg_scroshi_atrt = run_scROSHI_mrg_w_atrt_data(sr = mrgGA, cols = my_cols, pt = 1, save_name = 'mrg_all_atrt', path_to_save = atacScroshiDir)
  )

  make(scroshi_plan, lock_envir = FALSE, lock_cache = FALSE)