assign_tumor_cells <- function(rna){
    # tumor cells were identified manually by 
  # inspecting singleR, infercnv, scroshi and markers
  # check identify_tumor_cells.R for more details 
t_clusters <- c(18,28, 26 )
macro_clusters <- 8
epi_clusters <- 21
liver_clusters <- 20
endo_clusters <- c(13,29)


rna$cell_identity <- rna$Subtype
rna$cell_identity[rna$RNA_snn_res.0.8 %in% t_clusters] <- 'T_cells'
rna$cell_identity[rna$RNA_snn_res.0.8 %in% macro_clusters] <- 'macrophage/monocyte'
rna$cell_identity[rna$RNA_snn_res.0.8 %in% epi_clusters] <- 'epithelial'
rna$cell_identity[rna$RNA_snn_res.0.8 %in% liver_clusters] <- 'liver'
rna$cell_identity[rna$RNA_snn_res.0.8 %in% endo_clusters] <- 'endothelial'
return(rna)
}
