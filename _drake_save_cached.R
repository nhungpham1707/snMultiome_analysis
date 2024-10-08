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

# loadd(rna_group_sgr)
# final_hm_rna <- harmony_n_plot(rna_group_sgr, batch_factor = 'library', theta = 0,
#    sigma = 0.1, save_path = "output/batchEffect/rna/harmony/test_final_hm")
# p <- DimPlot(final_hm_rna, group.by = 'Subtype',
# cols = my_cols)
# savePlot('output/batchEffect/rna/harmony/test_final_hm/hm_rna_subtype.png')

# loadd(rna_w_tumor_label)
# saveRDS(file = 'output/sc_RNA/merge_all/rna_hm.RDS', rna_w_tumor_label)

# loadd(rna_nohm_tumor_label)
# saveRDS(file = 'output/sc_RNA/merge_all/rna.RDS', rna_nohm_tumor_label)

# loadd(atac_group_sgr)
# saveRDS(file = 'output/sc_atac/merge_all/atac.RDS',atac_group_sgr)

# loadd(final_hm_atac_umap)
# saveRDS(file='output/sc_atac/merge_all/atac_hm.RDS', final_hm_atac_umap)

# loadd(hthyDim_adrenal)
# saveRDS(file='output/healthy_data/hthyDim_adrenal', hthyDim_adrenal)
# loadd(atac_hm_w_tumor_label)
# saveRDS(file = 'output/sc_atac/merge_all/atac_hm_w_tumor_label.RDS', atac_hm_w_tumor_label)

# loadd(atac_nohm_tumor_label)
# saveRDS(file = 'output/sc_atac/merge_all/atac_nohm_tumor_label.RDS', atac_nohm_tumor_label)

# loadd(hg38)
# saveRDS(file = 'output/hg38.RDS', hg38)


# loadd(train_dsc_atac40k)
# saveRDS(file='output/logistic_regression/train_dsc_atac_40k.RDS', train_dsc_atac40k )

# loadd(new_atachm_mx)
# saveRDS(file = 'output/logistic_regression/atac_hm_features_above_300cells.RDS', new_atachm_mx )

# message('load data')
# loadd(atac_hthymrgDim)

# message('read dsc_atac')
# dsc_atac <- readRDS('output/logistic_regression/dsc_only_overlap_countMx.RDS')

# message('add cell type')
# dsc_atac$cell_type <- atac_hthymrgDim$cell_type

# message('save rds')
# saveRDS(file = 'output/logistic_regression/dsc_atac_only_overlap.RDS', dsc_atac)

# loadd(dsc_atac_ident)
# saveRDS(file = 'output/logistic_regression/dsc_atac_ident.RDS', dsc_atac_ident)

loadd(rna_hthymrg_clus)
saveRDS(file = 'output/logistic_regression/dscRna.RDS', rna_hthymrg_clus)
message('finished!')

