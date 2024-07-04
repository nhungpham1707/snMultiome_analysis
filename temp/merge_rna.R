#!/usr/bin/env Rscript

setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

# filename <- '25012024_all_multiome_lib.csv'
# metadata <- getData(filename, delim = ',')
# lib <- unique(metadata$name)
# path <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/output_19_12_2023/sc_RNA/processing/'
# name <- paste0('_normalize_isoutlier_filter_rna_data_w_features_scale_pca_umap_cluster.RDS')

# srMrg <- readRDS(paste0(path, lib[1],name))

# for (i in 2:length(lib)){
#     message(paste('merge', i, 'samples'))
#     sr <- readRDS(paste0(path, lib[i],name))
#     srMrg <-merge(x = srMrg, y= sr, add.cell.ids = c('',lib[i]))
#     saveRDS(srMrg, file = paste0(rnaMrgDir, '/merge_', i, 'samples.RDS')) 
# }
# message('normalizing merge')

rna_plan <- drake_plan(
srMrg = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_RNA/merge_all/merge_25samples.RDS'),
mrgRnaNor = NormalizeData(srMrg, normalization.method = "LogNormalize"),
mrgRnascale = ScaleData(mrgRnaNor, features = rownames(mrgRnaNor), verbose = FALSE),
mrgRnavar = FindVariableFeatures(mrgRnascale),
mrgRnaSct = SCTransform(mrgRnavar,verbose = FALSE, variable.features.n = 3000, assay = "rna"),
mrgRnapca = RunPCA(mrgRnaSct, verbose = TRUE, npcs = 50),
mrgRnaUmap = RunUMAP(mrgRnapca, dims = 1:30, n.neighbors = 30),
save_norMrgRna = saveRDS(mrgRnaUmap, paste0(rnaMrgDir, '/mrgRna.RDS')),
mrgRnaClu = FindClusters(mrgRnaUmap, algorithm = 3, resolution = 2, verbose = FALSE),
p_rnamrg = DimPlot(mrgRnaClu, raster = FALSE, cols = my_cols),
savep_rnamrg = savePlot(paste0(rnaMrgFigDir, '/mrgRNA_cluster.png'), p_rnamrg)
# mrgRna_wsnn = RunUMAP(mrgRnaSct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_"),
# mrgRna_wsnnClu = FindClusters(mrgRna_wsnn, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE),
# p_mrgRnaWsnn = DimPlot(mrgRna_wsnnClu, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5),
# savep_mrgRnawsnn = savePlot(paste0(rnaMrgFigDir, '/mrgRNA_cluster_wsnn.png'),p_mrgRnaWsnn)
)

make(rna_plan, lock_cache = FALSE)
message('finish mrg and plot rna')