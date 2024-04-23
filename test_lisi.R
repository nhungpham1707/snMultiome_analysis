## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

loadd(hm_atac_1338e8a3)
lisi_plan <- drake_plan(
lisi_hm_atac_1338e8a3 = calculate_lisi_from_sr(hm_atac_1338e8a3, batch = 'library')
)
make(lisi_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)

print('finished!')