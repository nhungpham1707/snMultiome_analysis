get_h5_link <- function(lb, metadata){
  link <- unique(metadata$data_link[which(metadata$name == lb)])
  h5Link <- paste0(base_data_dir, '/', link, 
                   '/outs/filtered_feature_bc_matrix.h5')
  h5LinkDf <- data.frame(link_to_h5 = h5Link,
                        lb = lb)# needed for later
}

create_GEX_seurat <- function(h5LinkDf){
  message ('-----------------creating rna seurat----------------')
  
  counts <- Read10X_h5(h5LinkDf$link_to_h5)
  sr <- CreateSeuratObject(counts = counts[["Gene Expression"]], 
                             assay = "RNA")
  sr$barcodes <- colnames(sr)
  sr$library <- h5LinkDf$lb
  return(sr)
}

remove_stress_genes <- function(sr, gene_list){
  message ('-----------------removing stress gene from rna seurat---')
  gene_to_retain <- setdiff(rownames(sr), gene_list)
  message (paste(length(setdiff(rownames(sr), gene_to_retain)), 
                 'stress genes were removed'))
  sr_rm <- subset(sr, feature = gene_to_retain)
  return(sr_rm)
}

check_mito_genes <- function(rna_sr){
  message('---check mt genes ----')
  rna_sr[["percent.mt"]] <- PercentageFeatureSet(rna_sr, pattern = "^MT-")
  return (rna_sr)
}

rna_visualize_metric <- function(rna_sr,  save_path){
    library_name = unique(rna_sr$library)
  VlnPlot(rna_sr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)
  ggsave(file = paste0(save_path, '/', Sys.Date(), "_", 
                       library_name, "scRNA_metric_visualization.png"),
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)
}

filter_rna_sr_w_isoutlier <- function(rna_sr, save_path){
    library_name = unique(rna_sr$library)
  message('--filter rna with isOutlier----')
  rna_sr_sce <- as.SingleCellExperiment(rna_sr)
  # nFeature lower and higher
  nFeature_outlier <- isOutlier(rna_sr$nFeature_RNA, type = 'both')
  nFeature_outlier_threshold <- attr(nFeature_outlier, 'thresholds')
  # plot
  plotColData(rna_sr_sce, y = 'nFeature_RNA',
            colour_by=I(nFeature_outlier)) + 
  theme(text = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 30))

  ggsave(file = paste0(save_path, 
                     '/', library_name, "scRNA_outlier_nFeature_RNA.png"),
       width = 1200 * reso/72, 
       height = 700 * reso/72, units ="px", dpi = reso)
  # percent.mt higher
  percent.mt_outlier <- isOutlier(rna_sr$percent.mt, type = 'higher')
  percent.mt_outlier_threshold <- attr(percent.mt_outlier, 'thresholds')
  if (percent.mt_outlier_threshold[2] == 0){
    message (paste('-----',library_name, 'has mt percent isoutlier threshold 0 so use pre-defined mt_percent_limit as', mt_percent_limit, '-------'))
    percent.mt_outlier_threshold[2] = mt_percent_limit}
  # plot
  plotColData(rna_sr_sce, y = 'percent.mt',
            colour_by=I(percent.mt_outlier)) + 
  theme(text = element_text(size = 30)) + 
  theme(axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 30))

  ggsave(file = paste0(save_path, 
                     '/', library_name, "scRNA_outlier_percent.mt.png"),
       width = 1200 * reso/72, 
       height = 700 * reso/72, units ="px", dpi = reso)
  # filter out 
  subset_rna_sr <- subset(rna_sr, subset = nFeature_RNA > nFeature_outlier_threshold[1] & 
                              nFeature_RNA < nFeature_outlier_threshold[2] & 
                              percent.mt < percent.mt_outlier_threshold[2])
  return (subset_rna_sr)
}
normalize_dim_sr <- function(sr){
  message ('-----------------normalizing rna seurate----------------')
  sr <- NormalizeData(sr, normalization.method = "LogNormalize")
  sr <- ScaleData(sr, features = rownames(sr), verbose = FALSE)
  sr <- FindVariableFeatures(sr)
  # sr <- SCTransform(sr,verbose = FALSE, variable.features.n = 3000)
  sr <- RunPCA(sr, verbose = TRUE, npcs = 50)
  sr <- RunUMAP(sr, dims = 1:30, n.neighbors = 30)
  return(sr)
}

normalize_dim_plot_sr <- function(sr, save_path, lib_name){
  message ('-----------------normalizing rna seurate----------------')
  sr <- NormalizeData(sr, normalization.method = "LogNormalize")
  sr <- ScaleData(sr, features = rownames(sr), verbose = FALSE)
  sr <- FindVariableFeatures(sr)
  # sr <- SCTransform(sr,verbose = FALSE, variable.features.n = 3000)
  sr <- RunPCA(sr, verbose = TRUE, npcs = 50)
    ElbowPlot(sr)
  ggsave(file = paste0(save_path, '/', lib_name, "elbow_plot.png"),
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)
    DimHeatmap(sr, dims = 1, cells = 500, balanced = TRUE)
  
  ggsave(file = paste0(save_path, '/', Sys.Date(), "_", lib_name, "dimheatmap_plot.png"),
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)
  
  DimPlot(sr, reduction = "pca") + NoLegend()
  
  ggsave(file = paste0(save_path, '/', Sys.Date(), "_", lib_name, "_dim_plot.png"),
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)

  sr <- RunUMAP(sr, dims = 1:30, n.neighbors = 30)
  return(sr)
}

clustering_rna_data <- function(sr, dims = 1:15){
  message ('----start clustering rna --------')
  sr <- FindNeighbors(sr, dims = dims)
  sr_w_cluster <- FindClusters(sr, 
                                                          resolution = c(0.1, 0.3, 0.5, 0.7, 1))
  return(sr_w_cluster)

}

plot_cluster <- function(sr, save_path, save_name){
DimPlot(sr, raster = FALSE, label = TRUE)
  ggsave(file = paste0(save_path, '/', Sys.Date(), "_", save_name, 
                       "_cluster.png"),
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)  
  }