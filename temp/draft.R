# remove cells with uncertain annotation in DESCARTES 
# (i.e. those with 'unknown', '?', 'xyz positive') and 
# cell types with less than 350 cells (based on previous prediction, 
# these cells do not have good performance)  
# also, cell types are grouped (i.e. lymphoid, myloid cells are grouped to immune cells) 
cleanDscAtac = readRDS('output/healthy_data/groupCelltype_nCleanDscAtac.RDS'), 
cleanDscAtacIdent = change_indent(cleanDscAtac, by = 'group_cell_type'),
cleanDscAtac_markers = FindAllMarkers(cleanDscAtacIdent, only.pos = T, logfc.threshold = 0.20),
topfeatures7CleanAtac = cleanDscAtac_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 7000, 
          wt = avg_log2FC),
  
features_to_keepcleanDsc= topfeatures7CleanAtac$gene,

  train_featurecleanDsc= intersect(features_to_keepcleanDsc, atac_features),
  sub_cleanDscAtac = subset(cleanDscAtacIdent, features = train_featurecleanDsc),
  sub_cleanDscAtac80 = sampling_sr(sub_cleanDscAtac, 80, class_col = 'cell_type', type = 'percent'),
  
  sub_cleanDscAtac_20bc= setdiff(colnames(sub_cleanDscAtac), colnames(sub_cleanDscAta80)),
  sub_cleanDscAtac20 =  subset(sub_cleanDscAtac, subset = cell_bc %in% sub_cleanDscAtac_20bc),
  
  # train ----
  train_cleanDscAtac= trainModel(GetAssayData(sub_cleanDscAtac80), class = sub_cleanDscAtac80$group_cell_type, maxCell = ncol(sub_cleanDscAtac80)),
  # test ----
  p_cleanDscAtac_test20 = predictSimilarity(train_cleanDscAtac, GetAssayData(sub_cleanDscAtac20), 
                                     classes = sub_cleanDscAtac20$group_cell_type, 
                                     logits = F),
  # predict ----
  sub_ataccleanDsc= new_atachmMx_colname[rownames(new_atachm_mx) %in% train_featurecleanDsc,], # 9760 features
  p_cleanDscAtac= predictSimilarity(train_cleanDscAtac, sub_ataccleanDsc, classes = atac_hm_tumor_nona$cell_identity, 
                              logits = F),
