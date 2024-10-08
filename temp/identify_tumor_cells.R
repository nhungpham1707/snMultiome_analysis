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


healthy_clusters <- c(28,26,8,18,29,13,21,20)
immune_clusters <- c(8,26,28,18)
t_clusters <- c(18,28, 26 )
macro_clusters <- 8
epi_clusters <- 21
liver_clusters <- 20
endo_clusters <- c(13,29)
fn_rms_clusters <-c(4) 
atrt_shh_clusters <- c(1,24)
atrt_tyr_clusters <- c(23,9)
mrt_clusters <- c(12,15,5,7)
sysa_clusters <- c(3,16, 0, 27)
maybe_rms_p3f_clusters <- c(14, 10)
maybe_rms_p3w_clusters <- 6
maybe_fn_rms_clusters <- 11
maybe_sysa_clusters <- c(2, 19)


rna$cluster_labels <- 'unknown'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% t_clusters] <- 'T_cells'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% macro_clusters] <- 'macrophage/monocyte'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% epi_clusters] <- 'epithelial'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% liver_clusters] <- 'liver'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% endo_clusters] <- 'endothelial'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% fn_rms_clusters] <- 'FN_RMS'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% atrt_shh_clusters] <- 'ATRT_SHH'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% atrt_tyr_clusters] <- 'ATRT_TYR'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% mrt_clusters] <- 'MRT'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% sysa_clusters] <- 'SySa'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% maybe_rms_p3f_clusters] <- 'maybe_P3F_RMS'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% maybe_rms_p3w_clusters] <- 'maybe_P3W_RMS'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% maybe_fn_rms_clusters] <- 'maybe_FN_RMS'
rna$cluster_labels[rna$RNA_snn_res.0.8 %in% maybe_sysa_clusters] <- 'maybe_sysa'

# investigate markers ---
Idents(rna) <- 'cluster_labels'
markers <- c('DES', 'MYOD1', 'MYOG', 'CD3D', 'CD8A', 'GNLY', 'APOE', 'CD14', 'C1QB', 'PAX5', 'MS4A1', 'CD19', 'VWF', 'PECAM1', 'CLDN5', 'SMARCB1')
DotPlot(rna, features= markers) + theme(axis.text.x = element_text(angle = 90))

# look like there are still contaminated cells in clusters

# check for each cluster again 

# RMS markers with only DES, MYOD1, and MYOG 

## maybe_FN_RMS cluster: 
unique(rna$Individual.ID[rna$cluster_labels == 'maybe_FN_RMS'])

## [1] "PMCID545AAO" 38 cells chr8 amp in cbioportal 
# "PMCID464AAA" 55 cells  chr8 amp in cbioportal 
# "PMCID665AAP" 90 cells no
# "PMCID155AAO" 724 yes chr8 
# "PMCID959AAM" 2796 yes chr8
# [6] "PMCID433AAN" 34 cells no chr8

# check infercnv LX187 for PMCID959AAM, heatmap show some amplification for chr8. aneuplidy score plot is weird. maybe cells above 200? not sure. for now don't take this infercnv. there were cells with abnormal infercnv from LX187 in cluster 4
Idents(rna) <- 'cluster_labels'
# genes from cbioportal for PMCID959AAM
DotPlot(rna, features = c('FANCC','CAMTA1'))
DotPlot(rna, features = c('UBR5','BAALC','AGO2','HEY1','NDRG1','PCM1','TONSL','RECQL4'))

# markers from the organoid paper 
rms_markers <- c('MYF6','MYF5', 'MYOD1', 'DES', 'MYOG', 'GPC3', 'FGFR4', 'IGF2', 'IGF1R', 'DUSP6', 'MYCN',  'MYL1', 'MYH3', 'TTN', 'TNNT2', 'TNNT3', 'TNNC2', 'RYR1', 'RYR3', 'VGLL2', 'MYL4', 'DMD', 'PAX7', 'MSTN', 'ERBB3', 'CXCR4', 'CD82', 'NCAM1', 'MEOX2', 'HMGA2', 'NOS1')
DotPlot(rna, features = rms_markers) + theme(axis.text.x = element_text(angle = 90))
# look for unique genes with CNV from cbioportal for FN-RMS compare to the other tumor : too much manual work 

# check other markers for fn-rms 
# HGA2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5657295/
# dotplot show strong expression in maybe-fn-rms but also in mrt and atrt-tyr 

DotPlot(rna, features = c('EPHA2', 'EED', "NELF", "CBS" , 'EPB41L4B'))

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10772578/
DotPlot(rna, features = c('GPC3', 'FGFR4', 'IGF2', 'IGF1R', 'DUSP6', 'MYCN'))
DotPlot(rna, features = c('PAX7', 'MYOD1', 'TTN', 'MYL1',  'MYH3'))
### check tumor by cluster resolution 1
rna <- rna_nodb_infer 
rna <- FindNeighbors(rna, dims = 1:30, assay = 'RNA', reduction='PCA')
rna <- FindClusters(rna, resolution = 1)

immune_res1 <- c(41,39,3,22)
liver_res1 <- c(29)
epi_res1 <- c(35)
endo_res1 <- c(47,45,11,40,46)
rna$tumor_labels <- 'potential_tumor'
rna$tumor_labels[rna$RNA_snn_res.1 %in% immune_res1] <- 'immune'
rna$tumor_labels[rna$RNA_snn_res.1 %in% liver_res1] <- 'liver'
rna$tumor_labels[rna$RNA_snn_res.1 %in% epi_res1] <- 'epithelial'
rna$tumor_labels[rna$RNA_snn_res.1 %in% endo_res1] <- 'endothelial'


# check epithelial cluster near MRT clusters 
rna_tree <- change_tree_label(rna, by = 'tumor_labels', save_name = 'output/cell_type/sc_rna/res1_tumor_tree.png', reduction.method = "PCA", assay.name = 'RNA', dims = 1:30, cluster.col = 'RNA_snn_res.1')

rna_tree <- change_tree_label(rna, by = 'cluster_labels', save_name = 'output/cell_type/sc_rna/res1_cluster_label_tree.png', reduction.method = "PCA", assay.name = 'RNA', dims = 1:30, cluster.col = 'RNA_snn_res.1')

# it seems cluster 27 is closer to healthy and cluster 38 too 
healthy_clusters_res1 <- c(40,11,47,45,46, 35,29,22,3,39,41, 27, 38)
rna$tumor_res1 <- 'unknown'
rna$tumor_res1[rna$RNA_snn_res.1 %in% healthy_clusters_res1] <- 'healthy'
p <- DimPlot(rna, group.by = 'tumor_res1', cols = c('orange', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/1.wo_db_healthy_res1.png',p)

rna$tumor_res1[rna$RNA_snn_res.0.8 %in% scroshi_clusters] <- 'tumor'
p <- DimPlot(rna, group.by = 'tumor_res1', cols = c('orange','blue', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/2.wo_db_scroshi_res1.png',p)

rna$tumor_res1[rna$RNA_snn_res.0.8 %in% infer_clusters] <- 'tumor'
p <- DimPlot(rna, group.by = 'tumor_res1', cols = c('orange','blue', 'grey'), pt.size = 1)
p
savePlot('output/cell_type/sc_rna/3.wo_db_infercnv_res1.png',p)


rna$tumor_labels <- 'unknown'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% t_clusters] <- 'T_cells'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% macro_clusters] <- 'macrophage/monocyte'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% epi_clusters] <- 'epithelial'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% liver_clusters] <- 'liver'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% endo_clusters] <- 'endothelial'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% fn_rms_clusters] <- 'tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% atrt_shh_clusters] <- 'tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% atrt_tyr_clusters] <- 'tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% mrt_clusters] <- 'tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% sysa_clusters] <- 'tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% maybe_rms_p3f_clusters] <- 'maybe_tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% maybe_rms_p3w_clusters] <- 'maybe_tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% maybe_fn_rms_clusters] <- 'maybe_tumor'
rna$tumor_labels[rna$RNA_snn_res.0.8 %in% maybe_sysa_clusters] <- 'maybe_tumor'

# label cell identity
rna$cell_identity <- rna$Subtype
rna$cell_identity[rna$tumor_labels == 'T_cells'] <- 'T_cells'
rna$cell_identity[rna$tumor_labels == 'macrophage/monocyte'] <- 'macrophage/monocyte'
rna$cell_identity[rna$tumor_labels == 'epithelial'] <- 'epithelial'
rna$cell_identity[rna$tumor_labels == 'liver'] <- 'liver'
rna$cell_identity[rna$tumor_labels == 'endothelial'] <- 'endothelial'
rna$cell_identity[rna$tumor_labels == 'maybe_tumor'] <- 'maybe_tumor'


# maybe sysa 
unique(rna$Individual.ID[rna$cluster_labels == 'maybe_sysa'])
"PMCID160AAA" "PMCID467AAP" "PMCID641AAN" "PMCID072AAO" "PMCID338AAA"

# lost of function of oncosuppressor genes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8233868/
sysa_markers <- c('TP53', 'RB1', 'NF1', 'PTEN', 'CDKN2A', 'SMARCB1',  'ATRX')
DotPlot(rna, features = sysa_markers)

# by overexpress MDM2
high_sysa <- c('EGFR', 'SSX', 'ERBB2', 'IGFBP2', 'IGF2' )
DotPlot(rna, features = high_sysa)
# it is not up in my data 

# targets of the fusion 
# https://pubmed.ncbi.nlm.nih.gov/23313505/ 
# gene from figure 3A
sysa_fusion_up_targets <- c('TBX3', 'DLX2', 'PAX3', 'EBF2', 'NKX2_2', 'PRPH', 'FGF9', 'FZD10', 'SOX3', 'AJAP1', 'EBF3', 'SIX1', 'SIX2', 'CABP7', 'CBX4', 'ETV5', 'PDGFRA', 'FOXD1', 'IGF2', 'TWIST1', 'DLX5', 'CRTAC1', 'GFRA2', 'FLRT2', 'HOXD1', 'EPHA4', 'EPHB3', 'EN2', 'MMP2', 'CDH11', 'SIX3', 'DLX1', 'MSX1', 'GBX2', 'PAX7', 'GRIK3', 'FGF19', 'CADM1', 'CACNA1G', 'CRLF1', 'ENC1', 'SHISA2', 'ADAMTS9', 'LRRC4C', 'ELMOD1', 'EEF1A2', 'ALDH1A2', 'NELL1', 'MLLT3') # no, did not show signal 

# FZD10 as target of the fusion https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142991
DotPlot(rna, features = 'FZD10') # no clear signal  

# check subtype of libraries in maybe sysa cluster 

maybe_sysa_lib <- unique(rna$library[rna$cluster_labels == 'maybe_sysa'])
maybe_sysa_lib
metadata[metadata$name %in% maybe_sysa_lib,]
table(rna$library[rna$cluster_labels == 'maybe_sysa'])
# mmost cells are from libraries with only sysa samples --> assign as sysa cluster 

# finalize the tumor identification

# if cells in healthy clusters, label them as healthy. otherwise, they get label as the sample where they are from 

