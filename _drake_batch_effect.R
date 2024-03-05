setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

batch_plan <- drake_plan(
  mrgAtac = readRDS('output/sc_atac/merge_all/mrgDim.RDS'),
    mrgAtacSce = as.SingleCellExperiment(mrgAtac),
    ## visualize batch ----
    # atac_visBatchDate = plotBatchVis(mrgAtacSce, batch = "Date.of.Library", save_path = batchAtacDir, col = my_cols),
    # atac_visBatchLib = plotBatchVis(mrgAtacSce, batch = 'library', save_path = batchAtacDir, col = my_cols),
    # atac_visBatchSample = plotBatchVis(mrgAtacSce, batch = 'sampleID', save_path = batchAtacDir, col = my_cols),
    # atac_visBatchGender = plotBatchVis(mrgAtacSce, batch = 'Gender', save_path = batchAtacDir, col = my_cols),

    ## calculate cms ---
   atac_cms = cms(mrgAtacSce, k =200, group = 'Date.of.Library', res_name = 'dj_date', n_dim = 30, cell_min = 100, dim_red = 'LSI'), 
    atac_cms_50 = cms(atac_cms, k =50, group = 'Date.of.Library', res_name = 'dj_date_k50', n_dim = 30, cell_min = 100, dim_red = 'LSI'), 
   atac_cms_lib = cms(atac_cms_50, k =200, group = 'library', res_name = 'dj_lib', n_dim = 30, cell_min = 100, dim_red = 'LSI'), 
    atac_cms_sid = cms(atac_cms_lib, k =200, group = 'sampleID', res_name = 'dj_sid', n_dim = 30, cell_min = 100, dim_red = 'LSI'), 
    plot_atac_cmsSid = visHist(atac_cms_sid),
    save_cms = saveRDS(atac_cms_sid, file = paste0(batchAtacDir, '/atac_cms.RDS')),
    savePAtacCms = save_plot(paste0(batchAtacDir, '/atac_cms.png'), plot_atac_cmsSid)
)

make(batch_plan, lock_cache = FALSE)
