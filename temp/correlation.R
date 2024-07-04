rna_fea <- findVariableFeatures(rna_w_tumor_label)
rna_g <- rna_fea@assays$RNA@var.features
hthy_fea <- findVariableFeatures(hthyDim_eye)
h_g <- hthy_fea@assays$RNA@var.features
common_g <- intersect(rna_g, h_g)
rna_sub <- rna_w_tumor_label[rownames(rna_w_tumor_label) %in% common_g,]
hty_sub <- hthyDim_eye[rownames(hthyDim_eye) %in% common_g,]