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


# querry_gr: data frame of grange object of querry (reference)
# hit_gr: data frame of grange object of target 
# count_hit: data frames, count matrix of target (output from GetAssayData)

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

assign_colname_newMx <- function(new_count_mx, sr){
  colnames(new_count_mx) <- colnames(sr)
  return(new_count_mx)
}


extract_atac_w_n_features <- function(n, atac_markers,atac_mx, atac_gr, dsc_gr, atac_sr){
topfeatures <- atac_markers %>% 
  group_by(cluster) %>% 
  top_n(n = n, 
        wt = avg_log2FC)

features_to_keep <- topfeatures$gene
sub_atac <- atac_mx[rownames(atac_mx) %in% features_to_keep,]

atac_gr_df = as.data.frame(atac_gr)
index_tokeep = paste0(atac_gr_df$seqnames, '-', atac_gr_df$start, '-', atac_gr_df$end) %in% features_to_keep
new_atac_gr = atac_gr[index_tokeep]
new_atachm_mx = makeCountMx_withSamePeaks_optimized3(dsc_gr,new_atac_gr, sub_atac)
colnames(new_atachm_mx) <- colnames(atac_sr)
return(new_atachm_mx)
}

extract_seurat_w_n_features <- function(n, markers, seurat){
  topfeatures <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = n, 
          wt = avg_log2FC)
  
  features_to_keep <- topfeatures$gene
  message(paste('there are', length(features_to_keep), 'features to extract'))
  sub_sr <-  subset(seurat, features = features_to_keep)
  return(sub_sr)
}

add_barcode_metadata <- function(sr){
  sr$cell_bc <- colnames(sr)
  return(sr)
}

createColPalete <- function(){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  Glasbey = glasbey.colors(32)
  all_col <- unique(c(col_vector, Glasbey))
  return(all_col)
}

# count feature per category, i.e. per cell type 
countFeatures <- function(seuratObj, category = 'cell_identity'){
  cells <- unique(seuratObj@meta.data[,category])
  featureCount <- matrix(data = NA, nrow = length(cells), ncol = 2)
  atacMx <- GetAssayData(seuratObj)

  for (i in 1:length(cells)){
    message(paste('process', i, 'type'))
    subcells <- colnames(seuratObj)[seuratObj@meta.data[,category] == cells[i]]
    subMx <- atacMx[,subcells]
    df_non0 <- subMx[apply(subMx[,-1], 1, function(x) !all(x==0)),]

    featureCount[i,1] <- cells[i]
    featureCount[i,2] <- nrow(df_non0)
  }

  featureCount[,2] <- as.numeric(featureCount[,2])
  colnames(featureCount) <- c('cells', 'count')
  return (featureCount)
}