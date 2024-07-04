setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

mrgRna_noNOr_all <- readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_RNA/merge_all/merge_25samples.RDS')

remove_CF_plan <- drake_plan(
genes_to_remove = unique(c(genelists$chr6HLAgenes, genelists$hemo, genelists$stress, genelists$ribo)),
gene_to_retain = setdiff(rownames(mrgRna_noNOr_all), genes_to_remove ),
rna_noCF = subset(mrgRna_noNOr_all, feature = gene_to_retain),
rna_noCF_nor = normalize_dim_plot_sr(rna_noCF, rnaMrgFigDir, lib_name = 'merge_noCF'),
rna_noCF_nor_clu = clustering_rna_data(rna_noCF_nor),
rna_noCF_meta = assign_meta(metadata, gexDemulMeta,
rna_noCF_nor_clu, save_name = paste0(rnaMrgDir, '/mrgRna_noCF_meta.RDS')),
dim_rna_noCF_lib = DimPlot(rna_noCF_meta, group.by = 'library', cols = my_cols, raster = FALSE,pt.size = 1),
save_dim_rna_noCF_lib = savePlot(paste0(rnaMrgFigDir, '/noCF_lib.png'), dim_rna_noCF_lib),
dim_rna_noCF_sub = DimPlot(rna_noCF_meta, group.by = 'Subtype', cols = my_cols, raster = FALSE,pt.size = 1),
save_dim_rna_noCF_sub = savePlot(paste0(rnaMrgFigDir, '/noCF_sub.png'), dim_rna_noCF_sub)
)

make(remove_CF_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
print('finished!')

loadd(rna_noCF_meta)
colnames(rna_noCF_meta@meta.data)