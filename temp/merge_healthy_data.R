#!/usr/bin/env Rscript
library(drake)
library(Seurat)

save_plot <- function(filename,
                      a_plot) {
  
  ggsave(file = filename,
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)
  
}
# normalize -----
normalize_atac <- function(atac_data){
sub_atac_data_outlier_filter_normalize <- RunTFIDF(atac_data) # normalization
sub_atac_data_outlier_filter_normalize <- FindTopFeatures(sub_atac_data_outlier_filter_normalize, min.cutoff = 'q0') # keep feautures in n cells. q0 mean in top 100% cell
sub_atac_data_outlier_filter_normalize <- RunSVD(sub_atac_data_outlier_filter_normalize)
}
# non-linear dimension reduction and clustering ---- 

dimension_reduction <- function(atac_data){
  sub_atac_data_outlier_filter_normalize <- RunUMAP(object= atac_data, reduction = 'lsi', dims = 2:30) # remove the 1st dim because it correlates w seq depth
  sub_atac_data_outlier_filter_normalize <- FindNeighbors(object= sub_atac_data_outlier_filter_normalize, reduction = 'lsi', dims = 2:30)
  sub_atac_data_outlier_filter_normalize <- FindClusters(object = sub_atac_data_outlier_filter_normalize, verbose = FALSE, algorithm = 3)
}

## 
healthyDir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/sc_atac/healthy_data'
dir.create(healthyDir, recursive = TRUE)
dataDir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/'
# test merge without drake 
merge_8_s <- readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/output/sc_atac/merge_all/healthy_merge_8_samples_seurat.RDS')
message('below are tissues that already merged')
unique(merge_8_s$tissue)

tissueL <- c('cerebellum', 'kidney', 'muscle',
'pancreas', 'placenta', 'spleen', 'stomach',
'thymus')
mergSr <- merge_8_s
for (i in 1:length(tissueL)){
    tsr <- readRDS(paste0(dataDir,tissueL[i], '_filtered.seurat.for_website.RDS'))
    mergSr <- merge(x = mergSr, y = tsr, add.cell.id = c('',tissueL[i]))
    n <- 8+i
    saveRDS(mergSr, paste0(healthyDir, '/merge_', n, '_samples.RDS'))
    message(paste('finish merging', n, 'samples'))
}


mrgNor <- normalize_atac(mergSr)
mrgDim <- dimension_reduction(mrgNor)
p = DimPlot(mrgDim, group.by = 'tissue'),
save_plot(paste0(healthyDir, '/all_healthy_data.png'), p)

message('finish merging and dimplot merge all healthy')