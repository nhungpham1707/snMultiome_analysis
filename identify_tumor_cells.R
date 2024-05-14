# various strategies were used for tumor cell identification.
# 1. singleR cell annotation. If cells are immune cells or stroma cells, they are not tumor cells
# 2. scROSHI with immune markers. similarly if cells are immune cells, they are not tumor cells. 
# 3. scROSHI with cancer markers.
# 4. inferCNV, if cells have known abnormal chromosome structure, they are potential tumor 
# 5. for ATRT & MRT: if they express SMARCB1, they may not be tumor cells 
# 6. for Sysa and RMS: check expression of marker genes and fusion target genes. 
loadd(hmRna_scroshi_atrt)
# add score columns ---
rna <- add_cancer_score_meta_cols(hmRna_scroshi_atrt)

# 1. assign score from singleR cell annotation ---
rna <- assign_sgr_cancer_score(rna)

# 2. assign score from scroshi ---
# this step is done manually by visualizing the umap and compare immune cells from singleR and scROSHI, if they assign to different cells in the same cluster --> that cluster contain immune cells

# 3. assign score from scroshi for cancer marker ---
rna <- assign_scroshi_cancer_score(rna)

# 4. assign score from infercnv res --- 
rna <- assign_infercnv_cancer_score(rna)

# 5. check SMARCB1 expression ---
Idents(rna) <- 'RNA_snn_res.0.8'
DotPlot(rna, feature = 'SMARCB1')
# manually infer the expression and assign expression as 0,1 to cluster 
# ATRT clusters: 
#     - 4,7: still have a bit of expression 
#     - 26: more SMARCB1 than others —> need to exclude 
#     - 24, 25: no 
# MRT clusters: 
#     - 12,15, 32, 18: a bit expression 
#     - 21: more SMARCB1 than other art clusters, could be healthy liver cells? 
#     - 5, 19, 9, 12: no 

rna$smarcb1 <- NA
rna$smarcb1[rna$RNA_snn_res.0.8 %in% c(4,7,24,25,12,15,32,18,5,19,9,12)] <- 1 # 1 for tumor cells because did not express smarcb1
rna$smarcb1[rna$RNA_snn_res.0.8 %in% c(26, 21)]  <- 0 # 0 for not tumor cell
DimPlot(rna, group.by = 'smarcb1')
# 6. check cancer markers for ATRT, RMS and Sysa ---
rna <- calculate_marker(rna, marker = ATRT_SHH, name = 'atrt_shh')

rna <- calculate_marker(rna, marker = ATRT_TYR, name = 'atrt_tyr')
# strong signal for atrt-tyr
rna <- calculate_marker(rna, marker = ATRT_MYC, name = 'atrt_myc')
p <- FeaturePlot(rna, features = 'myc1', label = T)
p <- FeaturePlot(rna, features = 'atrt_myc1', label = T)
p # no clear signal
savePlot('output/cell_type/sc_rna/markers/atrt_myc.png', p)

rna <- calculate_marker(rna, marker = RMS, name = 'rms')
p <- FeaturePlot(rna, features = 'rms1', label = T)
p # no clear signal
savePlot('output/cell_type/sc_rna/markers/rms.png', p)
rna <- calculate_marker(rna, marker = Sysa, name = 'sysa')
p <- FeaturePlot(rna, features = 'sysa1', label = T)
p
savePlot('output/cell_type/sc_rna/markers/sysa.png', p)
rna <- calculate_marker(rna, marker = sysa_fusion_targets, name = 'sysa_target')
p <- FeaturePlot(rna, features = 'sysa_target1', label = T)
p # clear signal for sysa
rna <- calculate_marker(rna, marker = rms_fusion_targets, name = 'rms_target')
p <- FeaturePlot(rna, features = 'rms_target1', label = T)
p
savePlot('output/cell_type/sc_rna/markers/rms_target.png', p)

# manually inspect the feature plots and identify clusters that express markers for each cancer 
# 6: atrt_tyr, 2,15,0: sysa, 5: rms 
rna$marker_score[rna$RNA_snn_res.0.8 %in% c(6,2,15,0,5)] <- 1

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

# look weird, why cluster 2 with sysa has many cells with 0 score 

# add score manually ---
rna$manual_tumor <- 'unknown'
DimPlot(rna, group.by = 'manual_tumor', cols = tumor_colors)

rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(8,18,13)] <- 'healthy'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'grey'))
p
savePlot('output/cell_type/sc_rna/1.w_sngr.png',p)

# w scroshi cancer markers 
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(7)] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'))
p
savePlot('output/cell_type/sc_rna/2.w_scroshi.png',p)

# w infercnv 
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(2)] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'))
p
savePlot('output/cell_type/sc_rna/3.w_infercnv.png',p)

# w rms target PAX3
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(6)] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'))
p
savePlot('output/cell_type/sc_rna/4.w_rms_markers.png',p)

# w sysa markers 
rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(0,16,1)] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'))
p
savePlot('output/cell_type/sc_rna/5.w_sysa_markers.png',p)

# w SMARCB1 expression 

rna$manual_tumor[rna$RNA_snn_res.0.8 %in% c(4,5,9,12,15,19,21,24,25,27,32)] <- 'tumor'
p <- DimPlot(rna, group.by = 'manual_tumor', cols = c('orange', 'blue', 'grey'))
p
savePlot('output/cell_type/sc_rna/6.w_no_SMARCB1_expression.png',p)
