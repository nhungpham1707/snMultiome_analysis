library(drake)
library(Seurat)
loadd(rna_nohm_tumor_label)
rna_nohm <- rna_nohm_tumor_label
DimPlot(rna_nohm, group.by = 'cell_identity', cols = my_cols)|
  DimPlot(rna_nohm, group.by = 'RNA_snn_res.0.3', label = T)

rna_nohm$cluster_cellidentity <- rna_nohm$cell_identity

rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'ecMRT' & rna_nohm$RNA_snn_res.0.3 == 7] <- 'ecMRT_1'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'ecMRT' & rna_nohm$RNA_snn_res.0.3 == 11] <- 'ecMRT_2'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'ecMRT' & rna_nohm$RNA_snn_res.0.3 == 13] <- 'ecMRT_3'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'ecMRT' & rna_nohm$RNA_snn_res.0.3 == 14] <- 'ecMRT_4'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'ecMRT' & rna_nohm$RNA_snn_res.0.3 == 16] <- 'ecMRT_5'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'ecMRT' & rna_nohm$RNA_snn_res.0.3 == 20] <- 'ecMRT_6'

rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'FN-eRMS' & rna_nohm$RNA_snn_res.0.3 == 4] <- 'FN-RMS_1'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'FN-eRMS' & rna_nohm$RNA_snn_res.0.3 == 10] <- 'FN-RMS_2'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'FN-eRMS' & rna_nohm$RNA_snn_res.0.3 == 17] <- 'FN-RMS_3'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'FN-eRMS' & rna_nohm$RNA_snn_res.0.3 == 18] <- 'FN-RMS_4'

rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'SySa' & rna_nohm$RNA_snn_res.0.3 == 0] <- 'SySa_1'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'SySa' & rna_nohm$RNA_snn_res.0.3 == 2] <- 'SySa_2'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'SySa' & rna_nohm$RNA_snn_res.0.3 == 15] <- 'SySa_3'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'SySa' & rna_nohm$RNA_snn_res.0.3 == 19] <- 'SySa_4'
rna_nohm$cluster_cellidentity[rna_nohm$cell_identity == 'SySa' & rna_nohm$RNA_snn_res.0.3 == 24] <- 'SySa_5'



rna_nohm$cluster_cellidentity[rna_nohm$cluster_cellidentity == 'SySa' ] <- 'SySa_1'
rna_nohm$cluster_cellidentity[rna_nohm$cluster_cellidentity == 'FN-eRMS' ] <- 'FN-RMS_4'
rna_nohm$cluster_cellidentity[rna_nohm$cluster_cellidentity == 'ecMRT' ] <- 'ecMRT_7'
DimPlot(rna_nohm, group.by = 'cluster_cellidentity', label = T)|
  DimPlot(rna_nohm, group.by = 'cell_identity', label = T)

loadd(dscRnaOverlaptrain80)
p <- predictSimilarity(dscRnaOverlaptrain80, GetAssayData(rna_nohm), classes = rna_nohm$cluster_cellidentity, logits = F)
similarityHeatmap(p)
# sysa_5 cells resemble immune cells (myeloid) from the logistic regression prediction. in rna with harmony, some of these
# cells are closer to macrophage clusters in UMAP. not sure what to do, for now remove them, there 
# are only 23 ccells anyway 
sysa5_cells <- colnames(rna_nohm)[rna_nohm$cluster_cellidentity == 'SySa_5']
length(sysa5_cells)
to_keep <- setdiff(colnames(rna_nohm), sysa5_cells)
length(to_keep)
length(intersect(colnames(rna_nohm), rna_nohm$m_barcode))
rna_nohm$m_barcode <- colnames(rna_nohm)
rna_nohm_nosysa5 <- subset(rna_nohm, m_barcode %in% to_keep )


unique(rna_nohm_nosysa5$Individual.ID[rna_nohm_nosysa5$cluster_cellidentity == 'FP-RMS'])
# "RMS000HVX" PMCID665AAM - from cbioportal it is P3F fusion 
rna_nohm_nosysa5$cluster_cellidentity[rna_nohm_nosysa5$cluster_cellidentity == 'FP-RMS'] <- 'FP-RMS (P3F)'
rna_nohm_nosysa5$cluster_cellidentity[grep('replapse', rna_nohm_nosysa5$cluster_cellidentity)] <- 'ATRT-MYC (relapse)'
rna_nohm_nosysa5$cluster_cellidentity[grep('BrainMet', rna_nohm_nosysa5$cluster_cellidentity)] <- 'ecMRT (metastasis)'

rna_nohm_nosysa5$cluster_cellidentity[grep('ATRT_MYC', rna_nohm_nosysa5$cluster_cellidentity)] <- 'ATRT-MYC'
rna_nohm_nosysa5$cluster_cellidentity[grep('ATRT_TYR', rna_nohm_nosysa5$cluster_cellidentity)] <- 'ATRT-TYR'
rna_nohm_nosysa5$cluster_cellidentity[grep('ATRT_SHH', rna_nohm_nosysa5$cluster_cellidentity)] <- 'ATRT-SHH'
rna_nohm_nosysa5$cluster_cellidentity[grep('T_cells', rna_nohm_nosysa5$cluster_cellidentity)] <- 'T cells'
DimPlot(rna_nohm_nosysa5, group.by = 'cluster_cellidentity', label = T)

p_nosysa5 <- predictSimilarity(dscRnaOverlaptrain80, GetAssayData(rna_nohm_nosysa5), classes = rna_nohm_nosysa5$cluster_cellidentity, logits = F)
similarityHeatmap(p_nosysa5)

saveRDS(file = 'output/sc_RNA/merge_all/rna_nohm_clustertumor_nosysa5.RDS', rna_nohm_nosysa5)

# remove sysa 5 from rna with harmony 
loadd(rna_w_tumor_label_newbc)
rna_hm <- rna_w_tumor_label_newbc
tokeep_hm <- setdiff(colnames(rna_hm), sysa5_cells)
rna_hm$m_barcode <- colnames(rna_hm)
rna_hm_nosysa5 <- subset(rna_hm, subset = m_barcode %in% tokeep_hm)
rna_hm_nosysa5$cell_identity[rna_hm_nosysa5$cell_identity == 'FP-RMS'] <- 'FP-RMS (P3F)'
DimPlot(rna_hm_nosysa5, group.by = 'cell_identity', cols = my_cols)
rna_hm_nosysa5$cell_identity[grep('replapse', rna_hm_nosysa5$cell_identity)] <- 'ATRT-MYC (relapse)'
rna_hm_nosysa5$cell_identity[grep('BrainMet', rna_hm_nosysa5$cell_identity)] <- 'ecMRT (metastasis)'
rna_hm_nosysa5$cell_identity[grep('FN-eRMS', rna_hm_nosysa5$cell_identity)] <- 'FN-RMS'
rna_hm_nosysa5$cell_identity[grep('ATRT_MYC', rna_hm_nosysa5$cell_identity)] <- 'ATRT-MYC'
rna_hm_nosysa5$cell_identity[grep('ATRT_TYR', rna_hm_nosysa5$cell_identity)] <- 'ATRT-TYR'
rna_hm_nosysa5$cell_identity[grep('ATRT_SHH', rna_hm_nosysa5$cell_identity)] <- 'ATRT-SHH'
rna_hm_nosysa5$cell_identity[grep('T_cells', rna_hm_nosysa5$cell_identity)] <- 'T cells'
rna_hm_nosysa5$cell_identity[rna_hm_nosysa5$cell_identity == 'FP-RMS'] <- 'FP-RMS (P3F)'

DimPlot(rna_hm_nosysa5, group.by = 'cell_identity', 
        cols = my_cols, pt.size = 0.8,
        label = T, repel = T)
saveRDS(file = 'output/sc_RNA/merge_all/rna_nohm_clustertumor_nosysa5.RDS', rna_hm_nosysa5)

