merge_sr_list <- function(sr_list, savepath){
  mergeSr <- sr_list[1]
  for (i in 2:length(sr_list)){
    tissue <- unique(sr_list$tissue)
    mergeSr <- merge(x = mergeSr, y = sr_list[i],
                merge.data = TRUE, 
                add.cell.ids = c('', tissue))
    saveRDS(mergeSr, paste0(savepath, '/merge_', i, '_tissue.RDS'))
  }
}

# for large dataset, merge required too much memory, to reduce it use integration with anchors 
# https://satijalab.org/seurat/archive/v4.3/integration_large_datasets

integration_w_anchors <- function(sr_list){
FindIntegrationAnchors(sr_list, reduction = 'rlsi')
sr.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
}
