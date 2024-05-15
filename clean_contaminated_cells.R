# lx101 : ATRT-TYR + Sysa 
sr <- gexNodb_LX101
lib <- 'LX101'

unknown_bc <- plot_unknow_cell(hmRna_scroshi_atrt, sr, lib = lib)

DimPlot(hmRna_scroshi_atrt, cells.highlight= unknown_bc)|DimPlot(hmRna_scroshi_atrt, group.by = 'Subtype', cols = my_cols)

DimPlot(hmRna_scroshi_atrt, group.by = 'RNA_snn_res.0.8', cols = my_cols, label = T)|DimPlot(hmRna_scroshi_atrt, group.by = 'Subtype')

# remove unknown cells with different subtypes in cluster 1,9,7,5 

index <- colnames(hmRna_scroshi_atrt) %in% unknown_bc

unknown_df <- data.frame(cluster = hmRna_scroshi_atrt$RNA_snn_res.0.8[index],
subtype = hmRna_scroshi_atrt$Subtype[index],
library = hmRna_scroshi_atrt$library[index])

# extract cluster of interest 
# cluster 1 is ATRT-shh 
cluster <- 1
sub_df <- unknown_df[unknown_df$cluster %in% cluster,]
table(sub_df$subtype)
to_keep <- setdiff(colnames(hmRna_scroshi_atrt), rownames(sub_df))
rna_rm <- subset(hmRna_scroshi_atrt, subset = m_barcode %in% to_keep)

DimPlot(hmRna_scroshi_atrt, group.by = 'Subtype', cols = my_cols)| DimPlot(rna_rm, group.by = 'Subtype', cols = my_cols)