# various strategies were used for tumor cell identification.
# 1. singleR cell annotation. If cells are immune cells or stroma cells, they are not tumor cells
# 2. scROSHI with immune markers. similarly if cells are immune cells, they are not tumor cells. 
# 3. scROSHI with cancer markers. 
# 4. inferCNV, if cells have known abnormal chromosome structure, they are potential tumor 
# 5. for ATRT & MRT: if they do not express SMARCB1, they maybe tumor cells 
# 6. for Sysa and RMS: check expression of marker genes and fusion target genes. 
# for each strategy, cells will get a score of 1 if they are classified as tumor, 0 if otherwise. These score will be summed up. In the end cells with the highest score will be most likely tumor cells. For example, if cells are classified as tumor in both singleR, scROSHi w immune markers, scroshi with cancer markers, infercnv and cancer specific markers, they will have a score of 5. 

# Nhung 14 05 2024

# loadd(hmRna_scroshi_atrt)
loadd(rna_nodb_infer)
rna <- rna_nodb_infer
# add score columns ---
rna <- add_cancer_score_meta_cols(rna)

lib_p <- DimPlot(rna, group.by = 'library', cols = my_cols)
savePlot('output/cell_type/sc_rna/wodb_library.png',lib_p)

sub_p <- DimPlot(rna, group.by = 'Subtype', cols = my_cols)
savePlot('output/cell_type/sc_rna/wodb_subtype.png',sub_p)

clus_p <- DimPlot(rna, group.by = 'RNA_snn_res.0.8', label = T, cols = my_cols)
savePlot('output/cell_type/sc_rna/wodb_cluster.png',clus_p)


# 1. assign score from singleR cell annotation ---
rna <- assign_sgr_cancer_score(rna)
DimPlot(rna, group.by = 'sgr_cancer_score')| DimPlot(rna, group.by = 'Subtype', cols  = my_cols)
p <- DimPlot(rna,group.by = 'sgr_cancer_score', cols = c('grey', 'red') )
p
savePlot('output/cell_type/sc_rna/1.wodb_sgr_score.png',p)

sgr_p <- DimPlot(rna,group.by = 'singleR_labels', cols = my_cols )

savePlot('output/cell_type/sc_rna/wodb_sgr.png',sgr_p)

sgr_gr_p <- DimPlot(rna,group.by = 'group_sgr_labels', cols = c('grey', 'black', 'green', 'yellow', 'blue', 'purple', 'salmon') )

savePlot('output/cell_type/sc_rna/wodb_group_sgr.png',sgr_gr_p)
# 1b. assign score from scroshi ---
# this step is done manually by visualizing the umap and compare immune cells from singleR and scROSHI, if they assign to different cells in the same cluster --> that cluster contain immune cells

healthy_clusters <- c(28,26,8,18,29,13,21,20)
# 2. assign score from scroshi for cancer marker ---
rna <- assign_scroshi_cancer_score(rna)
DimPlot(rna, group.by = 'scroshi_cancer_score')
scroshi_p <- DimPlot(rna, group.by = 'celltype_final', cols = c('blue', 'red', 'black', 'grey'), pt.size = 1)
scroshi_p
savePlot('output/cell_type/sc_rna/2.wodbscroshi_cancer.png',scroshi_p)
scroshi_clusters <- 9 
# 3. assign score from infercnv res --- 
rna <- assign_infercnv_cancer_score(rna)
infer_p <- DimPlot(rna, group.by = 'infer_cancer_score', cols = c('grey', 'red'), pt.size = 1)
savePlot('output/cell_type/sc_rna/3.wodb_infer.png', infer_p)
infer_clusters <-c(4) 

# 4. check SMARCB1 expression for ATRT and MRT---
Idents(rna) <- 'RNA_snn_res.0.8'
smarcb1_p <- DotPlot(rna, feature = 'SMARCB1')
savePlot('output/cell_type/sc_rna/4.wodb_smarcb1.png', smarcb1_p)
# manually infer the expression and assign expression as 0,1 to cluster 
# ATRT clusters: 
# 1,24, 9, 23 : no expression
# MRT clusters: 
# 12,15,5,7 : no expression
#     - 20: more SMARCB1 than other art clusters, more confident that it is healthy liver cells 
smarcb1_clusters <- c(1,24,9,23,12,15,5,7)

rna$smarcb1 <- NA
rna$smarcb1[rna$RNA_snn_res.0.8 %in% c( 1,24,9,23,12,15,5,7)] <- 1 # 1 for tumor cells because did not express smarcb1
smarcb1_dim_p <- DimPlot(rna, group.by = 'smarcb1', pt.size = 1, cols = c('red', 'grey'))
savePlot('output/cell_type/sc_rna/4.wodb_smarcb1_dim.png', smarcb1_dim_p)

# 5. check cancer markers for ATRT, RMS and Sysa ---
rna <- calculate_marker(rna, marker = ATRT_SHH, name = 'atrt_shh')
FeaturePlot(rna, features = 'atrt_shh1')
FeaturePlot(rna, features = ATRT_SHH)
rna <- calculate_marker(rna, marker = ATRT_TYR, name = 'atrt_tyr')
p <- FeaturePlot(rna, features = 'atrt_tyr1', pt.size = 1)
savePlot('output/cell_type/sc_rna/4.wodb_atrt_tyr_markers.png', p)

# strong signal for atrt-tyr. more confident for atrt_tyr cluster 
atrt_tyr_cluster <- 9

rna <- calculate_marker(rna, marker = ATRT_MYC, name = 'atrt_myc')
p <- FeaturePlot(rna, features = 'atrt_myc1', pt.size = 1)
p # no clear signal
savePlot('output/cell_type/sc_rna/4.wodb_atrt_myc_markers.png', p)

rna <- calculate_marker(rna, marker = RMS, name = 'rms')
p <- FeaturePlot(rna, features = 'rms1', pt.size = 1)
p # no clear signal
FeaturePlot(rna, features = RMS)

savePlot('output/cell_type/sc_rna/4.wodb_rms_markers.png', p)
rna <- calculate_marker(rna, marker = Sysa, name = 'sysa')
p <- FeaturePlot(rna, features = 'sysa1',pt.size = 1)
p
savePlot('output/cell_type/sc_rna/4.wodb_sysa_markers.png', p)
p # clear signal for sysa
sysa_clusters <- c(3,16, 0, 27)
# manually inspect the feature plots and identify clusters that express markers for each cancer 
# 6: atrt_tyr, 2,15,0: sysa, 5: rms 
rna$marker_score[rna$RNA_snn_res.0.8 %in% c(sysa_clusters, atrt_tyr_cluster)] <- 1

# check fusion target ---
## sysa fusion --- 
# ref https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8817899/
sysa_fusion_targets <- read.csv('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/sysa_fusion_targets_NIHMS1727124.csv', sep = ';')
sysa_target_up <- sysa_fusion_targets$Direct.targets_up
sysa_target_up <- sysa_target_up[nchar(sysa_target_up) > 0]
rna <- calculate_marker(rna, marker = sysa_target_up, name = 'sysa_target')
p <- FeaturePlot(rna, features = 'sysa_target1', pt.size = 1)
p # not clear signal 


sysa_markers <- read.csv('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/markers/NIHMS1727124-s3a.csv', sep = ';')

sysa_epi <- sysa_markers$Epithelial
sysa_mes <- sysa_markers$Mesenchymal 
rna <- calculate_marker(rna, marker = sysa_epi, name = 'sysa_epi')
p <- FeaturePlot(rna, features = 'sysa_epi1', pt.size = 1)
p 


rna <- calculate_marker(rna, marker = sysa_mes, name = 'sysa_mes')
p_mes <- FeaturePlot(rna, features = 'sysa_mes1', pt.size = 1)
p | p_mes

# rms target ---
rms_target <- read.csv('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/PAX3FOXO_targets.csv', sep = ';', header = F)
rna <- calculate_marker(rna, marker = rms_target[,1], name = 'rms_target')
p <- FeaturePlot(rna, features = 'rms_target1', pt.size = 1)
p


# summary all scores -----
score_df <- data.frame(
    sgr = rna$sgr_cancer_score, 
    infer  = rna$infer_cancer_score, 
    scroshi = rna$scroshi_cancer_score, 
    marker = rna$marker_score,
    smarcb1 = rna$smarcb1)
rna$final_cancer_score <- rowSums(score_df, na.rm = T)
tumor_colors <- c('grey', '#ADE8F4','#48CAE4','#0096C7','#0077B6','#03045E')
 DimPlot(rna, group.by = 'final_cancer_score', cols = tumor_colors, pt =1)



# add score manually ---------------------
rna$manual_tumor <- 'unknown'
tumormap_p <- DimPlot(rna, group.by = 'manual_tumor', cols = tumor_colors, pt.size = 1)
savePlot('output/cell_type/sc_rna/tumor_map.png', tumormap_p)
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% healthy_clusters] <- 'healthy'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/1.nodb_healthy.png',p)

# w scroshi cancer markers 
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% scroshi_clusters] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/2.nobd_w_scroshi.png',p)

# w infercnv 
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% infer_clusters] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/3.w_infercnv.png',p)

# w sysa markers 
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% sysa_clusters] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/4.wodb_w_sysa_markers.png',p)

# w SMARCB1 expression 

rna$manual_tumor[rna$RNA_snn_res.0.8 %in% 
smarcb1_clusters] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/4.nodb_w_no_SMARCB1_expression.png',p)

# w rms target PAX3
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(6)] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'))
p
savePlot('output/cell_type/sc_rna/4.w_rms_markers.png',p)




