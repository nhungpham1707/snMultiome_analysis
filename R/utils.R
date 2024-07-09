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

set_default_assay <- function(sr, assay ){
  DefaultAssay(sr) <- assay
  return(sr)
}

# some atacseq data such as descartes are seurat object but are not chromatin assay. This prevents the use of built-in functions for chromatin assay such as GetAssayData or granges. this func convert such seurat into seurat with chromatin asay. Annotation file can be generated from 'gethg38annotation' func 
createSrWChromatinAssay <- function(sr, annotation){
  sr_chr <-  CreateChromatinAssay(
    counts = GetAssayData(sr),
    sep = c(":", "-"),
    annotation = annotation,
    min.cells = 10,
    min.features = 200)
  
  sr_chr <- CreateSeuratObject(counts = sr_chr,
                               assay = "peaks" )
  
  metadata <- sr@meta.data
  rownames(metadata) <- colnames(sr)
  sr_chr <- AddMetaData(sr_chr, metadata)
  return(sr_chr)
}

# add cell barcodes after merging as m_barcode in metadata 
add_mbarcode <- function(sr){
  sr$m_barcode <- colnames(sr)
  return(sr)
}

# remove na cells in a column in metadata 
# metadata_colname is string of colname
remove_na_cells_in_metadata <- function(sr, metadata_colname){
  na_cells = colnames(sr)[is.na(sr@meta.data[,metadata_colname])]
  cells_to_keep = setdiff(colnames(sr), na_cells)
  sr_add_mbarcode = add_mbarcode(sr)
  sr_wo_na = subset(sr_add_mbarcode, subset = m_barcode %in% cells_to_keep)
  return (sr_wo_na)
}

# to compare scATACseq from different sources
#  with different peak coordinate,
# we need to convert them into 
# the same peaks. 
# Here we get the overlap peaks, 
# and sum up all small fragments in the
# target data that are overlap 
# with a peak in the reference data
# querry_gr: granges object of the reference data. can be created with granges(sr_w_chromatin_assay)
# hit_gr: granges object of the target data to change peak coordinate
# count_hit: count matrix of the target data (can be created with GetAssayData(sr)) 


# querry_gr, hit_gr, count_hit are data frames

makeCountMx_withSamePeaks <- function(querry_gr, hit_gr, count_hit){
  olap <- findOverlaps(querry_gr,hit_gr)
  olap_df <- as.data.frame(olap) # 568396 peaks 
  unique_olap_querry <- unique(olap_df[,1]) # 262577 peaks 

  new_count_mx <- c()
  names <- c()
  for (i in 1:length(unique_olap_querry)){
    message(paste('running', i, 'peaks out of', length(unique_olap_querry), 'unique overlap peaks'))
    querry_index <- unique_olap_querry[i]
    hit_index <- olap_df[olap_df[,1]==querry_index,2]

    # change coordinate to that in querry (dsc)
    new_coor <- paste0(querry_gr[querry_index]@seqnames,'-', querry_gr[querry_index]@ranges)

    if (length(hit_index) > 1){
      new_value <- colSums(count_hit[hit_index,])
    } else {
      new_value <- count_hit[hit_index,]
    }
    new_count_mx <- rbind(new_count_mx, new_value)
    names <- rbind(names, new_coor)
    rownames(new_count_mx) <- names[,1]
    }
  return (new_count_mx)
}


makeCountMx_withSamePeaks_optimized <- function(querry_gr, hit_gr, count_hit){
  # Find overlaps
  olap <- findOverlaps(querry_gr, hit_gr)
  olap_df <- as.data.frame(olap)
  unique_olap_querry <- unique(olap_df[,1])
  
  # Pre-allocate memory
  n <- length(unique_olap_querry)
  num_cols <- ncol(count_hit)
  new_count_mx <- matrix(0, n, num_cols)
  names <- character(n)
  total_hit <- length(unique_olap_querry)
  # Loop through unique overlap queries
  for (i in seq_along(unique_olap_querry)){
    message(paste('process', i, 'out of', total_hit, 'overlap peaks'))
    querry_index <- unique_olap_querry[i]
    hit_index <- olap_df[olap_df[,1] == querry_index, 2]
    
    # Change coordinate to that in querry (dsc)
    new_coor <- paste0(seqnames(querry_gr[querry_index]), '-', ranges(querry_gr[querry_index]))
    
    # Summarize counts
     if (length(hit_index) > 1){
      new_value <- colSums(count_hit[hit_index,])
    } else {
      new_value <- count_hit[hit_index,]
    }
    
    new_count_mx[i, ] <- new_value
    names[i] <- new_coor
  }
  
  # Convert matrix to data frame
  new_count_mx <- as.data.frame(new_count_mx)
  rownames(new_count_mx) <- names
  return(new_count_mx)
}



makeCountMx_withSamePeaks_optimized2 <- function(querry_gr, hit_gr, count_hit) {
  # Find overlaps
  olap <- findOverlaps(querry_gr, hit_gr)
  olap_df <- as.data.frame(olap)
  
  # Extract unique querry indices
  unique_olap_querry <- unique(olap_df[, 1])
  
  # Pre-allocate memory
  n <- length(unique_olap_querry)
  num_cols <- ncol(count_hit)
  new_count_mx <- matrix(0, n, num_cols)
  names <- character(n)
  
  # Vectorize computation for new_coor
  seqnames_vec <- as.character(seqnames(querry_gr))
  ranges_vec <- as.character(ranges(querry_gr))
  
  # Create a list of indices for each unique querry index
  hit_indices_list <- split(olap_df[, 2], olap_df[, 1])
  
  # Loop through unique overlap queries
  for (i in 10) {
    querry_index <- unique_olap_querry[i]
    hit_index <- hit_indices_list[[as.character(querry_index)]]
    
    # Change coordinate to that in querry (dsc)
    new_coor <- paste0(seqnames_vec[querry_index], '-', ranges_vec[querry_index])
    
    # Summarize counts
    if (length(hit_index) > 1) {
      new_value <- colSums(count_hit[hit_index, ])
    } else {
      new_value <- count_hit[hit_index, ]
    }
    
    new_count_mx[i, ] <- new_value
    names[i] <- new_coor
  }
  
  # Convert matrix to data frame
  new_count_mx <- as.data.frame(new_count_mx)
  rownames(new_count_mx) <- names
  return(new_count_mx)
}



makeCountMx_withSamePeaks_optimized3 <- function(querry_gr, hit_gr, count_hit) {
  # Find overlaps
  olap <- findOverlaps(querry_gr, hit_gr)
  
  # Extract Hits and queryHits directly
  hits <- olap@to
  query_hits <- olap@from
  
  # Extract unique query indices
  unique_query <- unique(query_hits)
  n <- length(unique_query)
  
  # Pre-allocate memory
  num_cols <- ncol(count_hit)
  new_count_mx <- matrix(0, n, num_cols)
  names <- character(n)
  
  # Vectorize computation for new coordinates
  seqnames_vec <- as.character(seqnames(querry_gr))
  ranges_vec <- as.character(ranges(querry_gr))
  
  # Loop through unique query indices
  for (i in seq_len(n)) {
    message(paste('process', i, 'out of', n, 'peaks'))
    query_index <- unique_query[i]
    hit_indices <- hits[query_hits == query_index]
    
    # Change coordinate to that in query (dsc)
    new_coor <- paste0(seqnames_vec[query_index], '-', ranges_vec[query_index])
    
    # Summarize counts using matrix subsetting
    if (length(hit_indices) > 1) {
      new_value <- colSums(count_hit[hit_indices, ])
    } else {
      new_value <- count_hit[hit_indices, ]
    }
    
    new_count_mx[i, ] <- as.numeric(new_value)
    names[i] <- new_coor
  }
  
  # Convert matrix to data frame
  new_count_mx <- as.data.frame(new_count_mx)
  rownames(new_count_mx) <- names
  return(new_count_mx)
}




library(Matrix)
library(foreach)
library(doParallel)

makeCountMx_withSamePeaks_optimized_parallel <- function(querry_gr, hit_gr, count_hit) {
  # Find overlaps
  olap <- findOverlaps(querry_gr, hit_gr)
  
  # Extract Hits and queryHits directly
  hits <- olap@to
  query_hits <- olap@from
  
  # Extract unique query indices
  unique_query <- unique(query_hits)
  n <- length(unique_query)
  
  # Pre-allocate memory
  num_cols <- ncol(count_hit)
  new_count_mx <- matrix(0, n, num_cols)
  names <- character(n)
  
  # Vectorize computation for new coordinates
  seqnames_vec <- as.character(seqnames(querry_gr))
  ranges_vec <- as.character(ranges(querry_gr))
  
  # Set up parallel backend
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Perform parallel computation with progress messages
  new_count_mx <- foreach(i = unique_query, .combine = rbind) %dopar% {
    hit_indices <- hits[query_hits == i]
    
    # Change coordinate to that in query (dsc)
    new_coor <- paste0(seqnames_vec[i], '-', ranges_vec[i])
    
    # Print progress message
    message(sprintf("Processing query index %d/%d", which(unique_query == i), length(unique_query)))
    
    # Summarize counts using matrix subsetting
    if (length(hit_indices) > 1) {
      new_value <- colSums(count_hit[hit_indices,])
    } else {
      new_value <- count_hit[hit_indices, ]
    }
    
    as.numeric(new_value)
  }
  
  # Convert matrix to data frame
  new_count_mx <- as.data.frame(new_count_mx)
  rownames(new_count_mx) <- names
  
  # Deregister parallel backend
  stopCluster(cl)
  registerDoSEQ()  # Switch back to sequential processing
  
  return(new_count_mx)
}


merge_sparseMx <- function(list_of_mx){
  mx1 <- list_of_mx[[1]]
  other_mx <- list_of_mx[[2:length(list_of_mx)]]
  mrgMx <- RowMergeSparseMatrices(mx1, other_mx)
  return(mrgMx)
}