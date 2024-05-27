# Some potential tumor clusters contain cells from another tumor type, maybe they are missassigned or doublet. 
# check with scdblFinder did not confirm that they are doublet. 
# hence, in this script, I check each cluster manually, if these contaminated cells have unknown (in gender addmodulescore) or unassigned )in souporcell) from demultiplexing step, remove them. 

# 21 05 2024, Nhung 

# all unknown cells: 
all_unknown <- read.csv('output/cell_type/sc_rna/clean_unknown/all_unknown_bcs_manual.csv') # this file was generated manually 

# for each cluster, if the unknown cells have different tumor type than the majority of the cells in that cluster, note them to remove 
rna <- hmRna_scroshi_atrt
clusters <- unique(rna$RNA_snn_res.0.8)
# exclude healthy clusters those that are classified as immune cells, endothelial cells and liver cells from singleR and scROSHi w immune markers
rna <- add_cancer_score_meta_cols(rna)
rna <- assign_sgr_cancer_score(rna)

DimPlot(rna, group.by = 'sgr_cancer_score')|DimPlot(rna, group.by = 'RNA_snn_res.0.8',label = T) # cluster
healthy_clusters <- c( 8, 18, 20, 21, 13, 29, 28, 26 )
tumor_clusters <- setdiff(clusters, healthy_clusters)

# get the majority tumor type per cluster 
i = 25
cluster_index <- rna$RNA_snn_res.0.8 == tumor_clusters[i]
tumor_type <- data.frame(table(rna$Subtype[cluster_index]))
tumor_type
lib_l <- data.frame(table(rna$library[cluster_index]))
lib_l
main_tumor_type <- tumor_type$Var1[tumor_type$Freq == max(tumor_type$Freq)]
main_tumor_type
cluster_bc <- colnames(rna)[cluster_index]
# unknown cells in the cluster 
unknown_in_cluster <- intersect(all_unknown$x, cluster_bc) 
print(paste('there are', length(unknown_in_cluster), 'unknown cells out of', length(cluster_bc), 'cells in cluster', tumor_clusters[i]))
DimPlot(rna, cells.highlight= unknown_in_cluster)|DimPlot(rna, group.by = 'RNA_snn_res.0.8', cols = my_cols, label = T)
# cluster 9: 783 unknown cells out of 5159 cells. most are atrt, and mrt. remove rms and sysa
sysa_remove <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
rms_remove <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster9_remove <- c(sysa_remove, rms_remove)

# cluster 5: 2404 out of 6423 cells are unknown. most are MRT and MRT_brainMet. remove those from sysa, FP-RMS, FN-RMS 
sysa_remove <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
rms_remove <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster5_remove <- c(sysa_remove, rms_remove)
# cluster1: 1958 out of 7114 unknown cells. most are atrt, mrt. remove sysa and rms 
sysa_remove <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
rms_remove <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster1_remove <- c(sysa_remove, rms_remove)
# cluster 31: only 87 cells, 58 are unknown. it also separate from other --> remove the whole cluster
cluster31_remove <- cluster_bc

# cluster 17: 467 unknown cells out of 1935 cells. keep only FN-eRMS
fn_rms <- cluster_bc[grep('FN-eRMS', rna$Subtype[cluster_index])]
cluster17_remove <- setdiff(cluster_bc, fn_rms)
# cluster 2: 3255 unknown cells out of 6880 cells. most are sysa. there are also many cells from atrt and mrt. for now just remove them 
sysa <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
cluster2_remove <- setdiff(cluster_bc, sysa)
DimPlot(rna, cells.highlight = cluster2_remove)

# cluster 15: 1223 unknown cells out of 2364 cells
sysa_remove <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
rms_remove <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster15_remove <- c(sysa_remove, rms_remove)
DimPlot(rna, cells.highlight = cluster15_remove)
# cluster 7:  713 unknown cells out of 5907 cells. most are atrt, mrt 
sysa_remove <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
rms_remove <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster7_remove <- c(sysa_remove, rms_remove)
DimPlot(rna, cells.highlight = cluster7_remove)
# cluster 12: 341 unknown cells out of 3398 cells. most are ecmrt 
mrt_keep <- cluster_bc[grep('MRT', rna$Subtype[cluster_index])]
cluster12_remove <- setdiff(cluster_bc, mrt_keep)
DimPlot(rna, cells.highlight = cluster12_remove)

# # cluster 22 : 314 sysa, and some cells from atrt-tyr, ecmrt, ecmrt_brainmet, fn-erms, fp-rms. 115 unknown cells, they are in cluster 22 but distribute in many locations in umap, it seems they are in the locaton with clusters of similar tumor type. very weird. for now i just remove them 
cluster22_remove <- cluster_bc

# cluster 16: 1383 unknown cells out of 2165 cells. most are sysa. 
sysa <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
cluster16_remove <- setdiff(cluster_bc, sysa)
# cluster 10: 348 unknown cells out of 4629 cells. most are rms 
rms <-  cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
rms_no_unknown <- setdiff(rms, unknown_in_cluster)
cluster10_remove <- setdiff(cluster_bc, rms_no_unknown)
# cluster 6: 3 unknown cells out of 5948 cells. most are fp-rms 
rms <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster6_remove <- setdiff(cluster_bc, rms)
# cluster 11: 1384 unknown cells out of 4028 cells. most are fn-rms 
fn_rms <- cluster_bc[grep('FN-eRMS', rna$Subtype[cluster_index])]
cluster11_remove <- setdiff(cluster_bc, fn_rms)
# cluster3: 3513 unknown cells out of 6852 cells. most are sysa
sysa <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
cluster3_remove <- setdiff(cluster_bc, sysa)

# 22-05-2024
# cluster19: 659 unknown cells out of 1220 cells. most are sysa 
sysa <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
cluster19_remove <- setdiff(cluster_bc, sysa)
# cluster 23: 307 unknown cells out of 504 cells. most are ATRT-TYR 
atrt_tyr <- cluster_bc[grep('ATRT_TYR', rna$Subtype[cluster_index])]
cluster23_remove <- setdiff(cluster_bc, atrt_tyr)

# cluster4: 3282 unknown cells out of 6691 cells. most are FN-eRMS. a weird cluster, fn-erms also in clusters with mrt and atrt
fn_erms <- cluster_bc[grep('FN-eRMS', rna$Subtype[cluster_index])]
cluster4_remove <- setdiff(cluster_bc, fn_erms)

# cluster 14: 356 unknown cells out of 2574 cells. most are rms 
rms <- cluster_bc[grep('RMS', rna$Subtype[cluster_index])]
cluster14_remove <- setdiff(cluster_bc, rms)
# cluster 24: 163 unknown cells out of 397 cells. only atrt_shh and mrt. no need to remove 
# cluster 0: 4669 unknown cells out of 8541 cells. most are sysa 
sysa <- cluster_bc[grep('SySa', rna$Subtype[cluster_index])]
cluster0_remove <- setdiff(cluster_bc, sysa)

# cluster 30: 92 unknown cells out of 113 cells. only ecMRT. no need to remove 
# cluster 25: 26 unknown cells out of 285 cells. most are mrt
mrt <- cluster_bc[grep('ecMRT', rna$Subtype[cluster_index])]
cluster30_remove <- setdiff(cluster_bc, mrt)

# cluster 27: 114 unknown cells out of 176 cells. sysa. no need to remove 

# cluster 32: 4 unknown cells out of 77 cells. only fn-erms, no need to remove  
combine_remove <- c(cluster19_remove, cluster23_remove, cluster4_remove,cluster14_remove, cluster0_remove, cluster30_remove)



all_remove <- c(cluster12_remove, cluster15_remove, cluster17_remove, cluster1_remove, cluster2_remove, cluster31_remove, cluster5_remove, cluster7_remove, cluster9_remove, cluster22_remove, cluster16_remove, cluster10_remove, cluster6_remove, cluster11_remove, cluster3_remove)
write.csv(file  = 'output/cell_type/sc_rna/clean_unknown/cluster12-17-15-1-2-31-5-7-9-22-16-10-6-11-3-remove.csv', all_remove, row.names = F)


half_remove <- read.csv('output/cell_type/sc_rna/clean_unknown/cluster12-17-15-1-2-31-5-7-9-22-16-10-6-11-3-remove.csv')

all_remove <- c(combine_remove, half_remove$x) # 9301 cells 
write.csv(file  = 'output/cell_type/sc_rna/clean_unknown/all_cluster_remove.csv', all_remove, row.names = F)


to_keep <- setdiff(colnames(rna), all_remove)
sub_rna <- subset(rna, subset = m_barcode %in% to_keep)

DimPlot(sub_rna, group.by = 'Subtype', cols = my_cols)

rna_cluster <- subset(rna, subset = RNA_snn_res.0.8 == tumor_clusters[i])
# # lx101 : ATRT-TYR + Sysa 
# sr <- gexNodb_LX101
# lib <- 'LX101'

# unknown_bc <- plot_unknow_cell(hmRna_scroshi_atrt, sr, lib = lib)

# DimPlot(hmRna_scroshi_atrt, cells.highlight= unknown_bc)|DimPlot(hmRna_scroshi_atrt, group.by = 'Subtype', cols = my_cols)

# DimPlot(hmRna_scroshi_atrt, group.by = 'RNA_snn_res.0.8', cols = my_cols, label = T)|DimPlot(hmRna_scroshi_atrt, group.by = 'Subtype')

# # remove unknown cells with different subtypes in cluster 1,9,7,5 

# index <- colnames(hmRna_scroshi_atrt) %in% unknown_bc

# unknown_df <- data.frame(cluster = hmRna_scroshi_atrt$RNA_snn_res.0.8[index],
# subtype = hmRna_scroshi_atrt$Subtype[index],
# library = hmRna_scroshi_atrt$library[index])

# # extract cluster of interest 
# # cluster 1 is ATRT-shh 
# cluster <- 1
# sub_df <- unknown_df[unknown_df$cluster %in% cluster,]
# table(sub_df$subtype)
# to_keep <- setdiff(colnames(hmRna_scroshi_atrt), rownames(sub_df))
# rna_rm <- subset(hmRna_scroshi_atrt, subset = m_barcode %in% to_keep)

# DimPlot(hmRna_scroshi_atrt, group.by = 'Subtype', cols = my_cols)| DimPlot(rna_rm, group.by = 'Subtype', cols = my_cols)


## after identifying tumor cells (identify_tumor_cells.R), it seems there are still contaminated cells 
# check again. for each cluster, get library list, if the same library has both cells from both tumor type in the same cluster, they are probably doublet 
rna <- rna_nodb_infer

healthy_clusters <- c( 8, 18, 20, 21, 13, 29, 28, 26 )
tumor_clusters <- setdiff(unique(rna$RNA_snn_res.0.8), healthy_clusters)

i = 1
message(paste('investigate cluster', tumor_clusters[i]))
cluster_index <- rna$RNA_snn_res.0.8 == tumor_clusters[i]
tumor_type <- data.frame(table(rna$Subtype[cluster_index]))
tumor_type
lib_l <- data.frame(table(rna$library[cluster_index]))
lib_l

table(rna$library[cluster_index], rna$Subtype[cluster_index])
main_tumor_type <- tumor_type$Var1[tumor_type$Freq == max(tumor_type$Freq)]
main_tumor_type
cluster_bc <- colnames(rna)[cluster_index]
# unknown cells in the cluster 
unknown_in_cluster <- intersect(all_unknown$x, cluster_bc) 
print(paste('there are', length(unknown_in_cluster), 'unknown cells out of', length(cluster_bc), 'cells in cluster', tumor_clusters[i]))
DimPlot(rna, cells.highlight= unknown_in_cluster)|DimPlot(rna, group.by = 'RNA_snn_res.0.8', cols = my_cols, label = T)
