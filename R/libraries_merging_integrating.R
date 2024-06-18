merge_sr_list <- function(sr_list, savepath){
  mergeSr <- sr_list[[1]]
  message(paste('the first tissue is', unique(sr_list[[1]]$tissue)))
  for (i in 2:length(sr_list)){
    tissue <- unique(sr_list[[i]]$tissue)
    message(paste('merge', i, 'sample', 'tissue', tissue))
    mergeSr <- merge(x = mergeSr, y = sr_list[[i]],
                merge.data = TRUE, 
                add.cell.ids = c('', tissue))
    saveRDS(mergeSr, paste0(savepath, '/merge_', i, '_tissue.RDS'))
  }
  return(mergeSr)
}

# for large dataset, merge required too much memory, to reduce it use integration with anchors 
# https://satijalab.org/seurat/archive/v4.3/integration_large_datasets

integration_w_anchors_subset <- function(sr_list, to_remove){
sub_list <- c()
for (i in 1:length(sr_list)){
  lib <- unique(sr_list[i]$library)
  if (lib %in% to_remove){
    
  }
} 
anchors <- FindIntegrationAnchors(sub_list, reduction = 'rlsi', anchor.features = 2000)
sr.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
}

normalize_anchors <- function(sr.integrated){
  DefaultAssay(sr.integrated) <- "integrated"
  atacSr <- RunTFIDF(sr.integrated) 
  atacSr <- FindTopFeatures(atacSr, min.cutoff = 'q0') 
  # keep feautures in n cells. q0 mean in top 100% cell
  atacSr <- RunSVD(atacSr)
  return(atacSr)
}
# merge pairwise to reduce memory use 
merge_pairwise <- function(sr_list, save_path){
  srMrg <- sr_list[[1]]
  for (i in 2:length(sr_list)){
    message(paste('merging', i, 'samples'))
    sr <- sr_list[[i]]
    lib <-  unique(sr$library)
    srMrg <- merge(x = srMrg, y = sr, add.cell.ids = c('',lib), merge.data = TRUE)
    saveRDS(srMrg, file = paste0(save_path, '/merge_', i, 'samples.RDS'))       
  }
  return(srMrg)
}