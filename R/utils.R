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