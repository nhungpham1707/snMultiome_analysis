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
    atac_visBatchDate = plotBatchVis(atac.sce, batch = "Date.of.Library", save_path = batchAtacDir),
    atac_visBatchLib = plotBatchVis(atac.sce, batch = 'library', save_path = batchAtacDir),
    atac_visBatchSample = plotBatchVis(atac.sce, batch = 'sampleID', save_path = batchAtacDir),
    atac_visBatchGender = plotBatchVis(atac.sce, batch = 'Gender', save_path = batchAtacDir),

    ## calculate cms ---
   atac_cms = cms(mrgAtacSce, k =30, group = 'Date.of.Library', res_name = 'dj_date', n_dim = 30, cell_min = 100, dim_red = 'lsi'), 
   atac_cms = cms(atac_cms, k =30, group = 'library', res_name = 'dj_lib', n_dim = 30, cell_min = 100, dim_red = 'lsi'), 
    atac_cms = cms(atac_cms, k =30, group = 'sampleID', res_name = 'dj_sid', n_dim = 30, cell_min = 100, dim_red = 'lsi'), 
    plot_atac_cmsSid = visHist(atac_csm),
    savePAtacCms = save_plot(paste0(batchAtacDir, '/atac_cms.png'), plot_atac_cmsSid)
)

make(batch_plan, lock_cache = FALSE)
