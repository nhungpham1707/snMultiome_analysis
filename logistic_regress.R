# try logistic regression to find cell origin 
# 

loadd(hthyDim_adrenal)
ref_sr <- hthyDim_adrenal
loadd(rna_w_tumor_label)
refMat <- GetAssayData(ref_sr, slot = 'counts')
classes <- ref_sr$cell_type
t_m <- trainModel(refMat,classes)
tgtData <- GetAssayData(rna_w_tumor_label, slot = 'counts')
p_classes <- rna_w_tumor_label$cell_identity
p_m <- predictSimilarity(t_m,tgtData,classes= p_classes,minGeneMatch=0.70)
png(filename = 'output/figures/logistic_regression/trained_kidney.png')
similarityHeatmap(p_m)
dev.off()
# endothelial cells from our scRNAseq are very similar to the vascular endothelial cells in the reference, which is a good positive control 

# test on 75% ref data to train and test on the rest 25% of the ref data 
# for each cell type select 75% cells randomly 
sub_sr <- sampling_sr(ref_sr, 75)
sub_data <- GetAssayData(sub_sr, slot = 'counts')
sub_classes <- sub_sr$cell_type  
t_sub_sr <- trainModel(sub_data, sub_classes)

# test on the remaining 25% cells, if the model was trained correctly, it should predict correctly the label of these cells
ref_sr$all_bc <- colnames(ref_sr)
to_keep <- setdiff(colnames(ref_sr), colnames(sub_sr))
test_sr <- subset(ref_sr, subset = all_bc %in% to_keep)
test_classes <- test_sr$cell_type
test_data <- GetAssayData(test_sr, slot = 'counts')
test_m <- predictSimilarity(t_sub_sr,test_data,classes= test_classes,minGeneMatch=0.70)
png(filename = 'output/figures/logistic_regression/test_adrenal.png')
similarityHeatmap(test_m)
dev.off()
# yes, cell types are correctly correlated to themselve again --> model seems to be fine 

# generate merg reference dataset from descartes
# for each tissue, get a maximum 800 cells per cell type, then merg all the tissue together 

sub_sr <- sampling_sr(sr, percent_to_keep = 800, type = 'number') 

