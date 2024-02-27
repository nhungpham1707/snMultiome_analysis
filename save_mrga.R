setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
library(drake)

loadd(mrgGA)
saveRDS(mrgGA, '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_atac/merge_all/mrgGA.RDS')

