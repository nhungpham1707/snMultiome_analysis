setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()
  
loadd('mrgGA')
saveRDS(mrgGA, paste0(atcMrgDir, '/', 'mrgGA.RDS'))