setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# library(SeuratData)
library(SeuratDisk)
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()


# loadd(rna_meta)
# loadd(rnaMrgNoNor_sgr)


# cluster <-rna_meta@active.ident
# png(paste0('output/cell_type/sc_rna/singler/singler_score_heatmap_.png'), width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
# plotScoreHeatmap(
# rnaMrgNoNor_sgr,
# cells.use = NULL,
# labels.use = NULL,
# clusters = cluster, #or NULL
# show.labels = TRUE,
# show.pruned = FALSE,
# max.labels = 40,
# normalize = TRUE,
# cells.order = NULL,
# order.by = "clusters",
# scores.use = NULL,
# calls.use = 0,
# na.color = "gray30",
# cluster_cols = FALSE,
# annotation_col = NULL, #ann_colors gave error
# show_colnames = FALSE,
# color = (grDevices::colorRampPalette(c("#D1147E", "white", "#00A44B")))(100),
# silent = FALSE,
# grid.vars = list()
# )
# dev.off()
# meta <- rnaMrgNoNor_sgr$pruned.labels
# names(meta) <- rownames(rnaMrgNoNor_sgr)
# rnaSr <- AddMetaData(rna_meta, meta, col.name = 'singler.pruned.labels')
# saveRDS(rnaSr, paste0(rnaMrgDir, '/mrgRna_w_singler.RDS'))

# print(colnames(rnaSr@meta.data))
# rnaSr <- readRDS(paste0(rnaMrgDir, '/mrgRna_w_singler.RDS'))
# keep <- "ecMRT" 
# rnaSub <- subset(rnaSr, subset = Subtype %in% keep )
# saveRDS(rnaSub, paste0(rnaMrgDir,'/ecMRT.rds'))
# SaveH5Seurat(rnaSub, filename = paste0(rnaMrgDir, '/ecMRT.h5Seurat'))
# Convert(paste0(rnaMrgDir, '/ecMRT.h5Seurat'), dest = 'h5ad')

# SaveH5Seurat(rnaSr, filename = paste0(rnaMrgDir, '/rna.h5Seurat'))
# Convert(paste0(rnaMrgDir, '/rna.h5Seurat'), dest = 'h5ad')
# print('done without error!')
# p <- DimPlot(rnaSub, group.by = 'Individual.ID', raster = FALSE, cols = my_cols)
# savePlot(paste0(rnaFigDir, '/ecMRT_patient.png'), p )
# p2 <- DimPlot(rnaSub, group.by = 'singler.pruned.labels', raster = FALSE, cols = my_cols)
# savePlot(paste0(rnaFigDir, '/ecMRT_singler.png'), p2 )


# rnaSr <- readRDS(paste0(rnaMrgDir, '/mrgRna_w_singler.RDS'))
# keep <- "ecMRT" 
# rnaSub <- subset(rnaSr, subset = Subtype %in% keep )
# saveRDS(rnaSub, paste0(rnaMrgDir,'/ecMRT.rds'))
# SaveH5Seurat(rnaSub, filename = paste0(rnaMrgDir, '/ecMRT.h5Seurat'))

loadd(mrgAtacDim)
SaveH5Seurat(mrgAtacDim, filename = paste0(atcMrgDir, '/atac.h5Seurat'))
Convert(paste0(atcMrgDir, '/atac.h5Seurat'), dest = 'h5ad')

colnames(mrgAtacDim@meta.data)

print('finished!')