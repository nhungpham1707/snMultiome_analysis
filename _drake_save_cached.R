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

# loadd( "rnaMrgNoNor_sgr")
# saveRDS(rnaMrgNoNor_sgr, paste0('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/cell_type/sc_rna/singler/mrgRna_noNOr_singleR.RDS'))

# loadd(rnaMrgSgr)
# colnames(rnaMrgSgr@meta.data)

# rna <- readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_RNA/merge_all/mrgRna_w_singler.RDS')
# print('read mrrgRna_w_singler is')
# colnames(sr@meta.data)

#

# Convert('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output_before_accident_cleanup/sc_RNA/merge_all/rna.h5ad', dest = 'h5seurat', overwrite = TRUE)
# sr <- LoadH5Seurat('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output_before_accident_cleanup/sc_RNA/merge_all/rna.h5seurat')
# message('save RDS---------------')
# saveRDS(sr, file = '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output_before_accident_cleanup/sc_RNA/merge_all/rna.RDS')

# loadd(gexClusSgr_LX049)
# save_h5ad(gexClusSgr_LX049, 'ikarus/', 'lx049')

loadd(rna_group_sgr)
final_hm_rna <- harmony_n_plot(rna_group_sgr, batch_factor = 'library', theta = 0,
   sigma = 0.1, save_path = "output/batchEffect/rna/harmony/test_final_hm")
p <- DimPlot(final_hm_rna, group.by = 'Subtype',
cols = my_cols)
savePlot('output/batchEffect/rna/harmony/test_final_hm/hm_rna_subtype.png')


message('finished!')