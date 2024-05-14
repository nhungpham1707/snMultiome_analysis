# by is string that match the colname in metadata that want to visualize. for ex, 'Subtype', 'sampleID'
get_number_of_cell <- function(sr, by){
  count <- as.data.frame(table(sr@meta.data[, c(by)], exclude = NULL))
  colnames(count) <- c(by, 'nCells')
  return (count)
}

change_label <- function(by, count){
  label <- count[,c(by)]
  new_label <- data.frame(old_label = label,
                          new_label = paste0(label, '_', count$nCells))
  return (new_label)
}

generate_new_label_metadata <- function(sr, new_label, by){
  new_label_metadata <- c()
  for (i in 1:length(sr@meta.data[,c(by)])){
    if (!is.na(sr@meta.data[i,c(by)])){
      index <- which(new_label$old_label == sr@meta.data[i,c(by)])
      new_label_metadata[i] <- new_label$new_label[index]  
    } else {
      new_label_metadata[i] <- new_label$new_label[is.na(new_label$old_label)]
    }
  }
  names(new_label_metadata) <- colnames(sr)
  return (new_label_metadata)
}

add_metadata_n_plot_w_new_label <- function(sr, by, new_label_metadata, save_path, col){
  new_sr <- AddMetaData(sr, metadata = new_label_metadata,
                        col.name = paste0(by,'_cell_counts')) 
  p <- DimPlot(new_sr, group.by = paste0(by,'_cell_counts'), raster = FALSE, pt.size = 1, cols = col)
  savePlot(filename = paste0(save_path, '/', by, '_w_cell_counts.png'), p)
  return(new_sr)
}

dimplot_w_nCell_label <- function(sr, by, save_path, col = hue_pal()(100)){
  count <- get_number_of_cell(sr, by)
  new_label <- change_label( by, count)
  new_label_metadata <- generate_new_label_metadata(sr, new_label, by)
  new_sr <- add_metadata_n_plot_w_new_label(sr, by, new_label_metadata, save_path, col = col)
  message (paste('-----finish plotting', by, 'with cell counts -----'))
  return (new_sr)
}

dimplotnSave <- function(sr, save_path,save_name, col = hue_pal()(100)){
  p <- DimPlot(sr, raster = FALSE, pt.size = 1, cols = col)
  savePlot(filename = paste0(save_path,'/', save_name, '.png'), p)
}

dimplotBynSave <- function(sr, by, save_path, save_name, col = hue_pal()(100)){
  p <- DimPlot(sr, group.by = by, raster = FALSE, pt.size = 1, cols = col)
  savePlot(filename = paste0(save_path,'/', save_name, '.png'), p)
}