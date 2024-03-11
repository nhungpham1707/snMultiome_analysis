## check cluster behavior ---- 
### Silhouette widths ------
# ref https://bioconductor.org/books/3.18/OSCA.advanced/clustering-redux.html#quantifying-clustering-behavior
# Cells with large positive silhouette widths are closer to other cells in the same cluster than to cells in different clusters. Thus, clusters with large positive silhouette widths are well-separated from other clusters
calculate_silhouette <- function(sce, reduce_method = 'LSI'){
sil.approx <- approxSilhouette(reducedDim(sce, reduce_method), clusters=colLabels(sce))
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, colLabels(sce), sil.data$other))
sil.data$cluster <- colLabels(sce)
return(sil.data)
}

### cluster purity ------
### The “clustering purity” is defined for each cell as the proportion of neighboring cells that are assigned to the same cluster, after some weighting to adjust for differences in the number of cells between clusteres. Well-separated clusters should exhibit little intermingling and thus high purity values for all member cells
# Some clusters have low purity values that may warrant more careful inspection - these probably represent closely related subpopulations.
calculate_purity <- function(sce, reduce_method = 'LSI'){
pure.sce <- neighborPurity(reducedDim(sce, reduce_method), colLabels(sce))
pure.data <- as.data.frame(pure.sce)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- colLabels(sce)
return(pure.data)
}

