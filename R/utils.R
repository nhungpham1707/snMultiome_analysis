# data_file_name is the sample_gender_meta_data 
getData <- function(filename, delim = ','){
  message('--reading metadata file----')
  data <- read.csv(paste0(base_data_dir, '/', filename), sep = delim, row.names = NULL)
  return (data)
}

# generate shorter lib name to be visible in drake graph
splitName <- function(name){
  name <- strsplit(name, split = '_')[[1]][1]
}

# save plot ---
savePlot <- function(filename,
                      a_plot) {
  
  ggsave(file = filename,
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)
  
}

# remove cells that fail from demultiplex 
# these cells cannot be assigned to any sample id or subtype 
remove_na_cells <- function(sr){
  na_cells <- colnames(sr)[is.na(sr$Subtype)] 
  sr$m_barcode <- colnames(sr)
  to_keep <- setdiff(sr$m_barcode, na_cells)
  sr_noNA <- subset(sr, subset = m_barcode %in% to_keep)
}

# make single cell experiment object
make_sce <- function(sr){
sce <- as.SingleCellExperiment(sr)
rowData(sce)$SYMBOL <- rownames(sce) 
colLabels(sce) <- sce$seurat_clusters
return(sce)
}

addLibSubcategory <- function(sr){
  sr$lbsb <- paste0(sr$library, '_', sr$Subtype)
  return (sr)
}

save_h5ad <- function(sr, save_path, save_name){
  save_h5 <- paste0(save_path, '/', save_name, '.h5Seurat')
  SaveH5Seurat(sr, filename = save_h5 )
  Convert(save_h5, dest = 'h5ad')
}

# add extra metadata on the missing patient 
add_missing_patient <- function(sr){
  sr$Individual.ID[sr$Subtype == 'ecMRT_BrainMet'] <- 'PMCID499AAQ'
  return (sr)
}

# check cell cycle genes 
# code copied from 
# https://github.com/scgenomics/scgenomics-public.github.io/blob/main/docs/05-confounders/05-confounders.R 
check_cell_cycle <- function(sr, save_path){
data(cc.genes)

## add cell cycle phase estimates
sr <- CellCycleScoring(object = sr,
        s.features =   cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes,
        assay='RNA') # make sure LogNormalized data is used!

## show the UMAP with type and phase information
p_type <-  DimPlot(sr, pt.size=1, group.by='Subtype', cols=my_cols)
p_phase <- DimPlot(sr, pt.size=1, group.by='Phase')
p_all <- p_type | p_phase
save_plot(paste0(save_path, '/cellcycle_genes.png'), p_all)
return (sr)
}

calculate_lisi_from_sr <- function(sr, batch){
  X <- GetAssayData(sr, slot = 'counts')
  metadata <- sr@meta.data[batch]
  res <- compute_lisi(X, meta_data)
  return (res)
}

change_indent <- function(sr, by){
  Idents(sr) <- by
  return(sr)
}

plot_dataset <- function(metadata, col_to_plot, save_name, cols = 'slateblue1'){
  metadata[,c(col_to_plot)][is.na(metadata[,(col_to_plot)])] <- ''
  stage <- data.frame(table(metadata[,col_to_plot]))
  colnames(stage) <- c('category', 'n')
  ggplot(stage, aes(category, n)) + 
  geom_bar(stat="identity", fill = cols) + coord_flip() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab('') + ylab('') + labs(title = col_to_plot)+ 
  theme(text = element_text(size = 30, hjust = 0))

  ggsave(filename = save_name, 
       width = 1200 * reso/72, height = 700 * reso/72, units ="px", dpi = reso)
}

# transfer labels from rna to atac or between rna before harmony and after harmony
assign_cross_labels <- function(des_sr, source_sr, label_col){
  source_data <- data.frame(category = source_sr@meta.data[,label_col],
    bc = colnames(source_sr))
  des_data <- source_data[source_data$bc %in% colnames(des_sr),]
  metadata <- des_data$category
  names(metadata) <- des_data$bc
  message(paste(length(des_data$bc), 'out of', length(source_data$bc), 'cells found in destination sr!'))
  
  des_sr <- AddMetaData(des_sr, metadata, col.name = label_col)
  return(des_sr)
}