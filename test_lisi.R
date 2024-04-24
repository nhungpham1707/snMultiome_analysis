## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

loadd(hm_atac_00d395cc)
lisi_hm_atac_00d395cc = calculate_lisi_from_sr(hm_atac_00d395cc, batch = 'library')
saveRDS(lisi_hm_atac_00d395cc, 'output/batchEffect/atac/harmony/lisi/hm_atac_00d395cc.RDS')


print('finished!')