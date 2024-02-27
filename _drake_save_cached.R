setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()
  
# loadd('mrgatacGA_singler')
# saveRDS(mrgatacGA_singler, paste0(atacCellSngRDir, '/', 'mrgatacGA_singler.RDS'))

# loadd(atac_anchors)
# saveRDS(atac_anchors, '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_atac/merge_all/anchors.RDS')

loadd( "rnaMrgNoNor_sgr")
saveRDS(rnaMrgNoNor_sgr, paste0('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/cell_type/sc_rna/singler/mrgRna_noNOr_singleR.RDS'))