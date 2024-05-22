# function to change the tip label of a hierarchical clustering tree 
# sr: seurat object
# by: which column in the sr@meta.data to 
# label the tree
# save_name: name to save plot with path+name
change_tree_label <- function(sr, by, save_name, reduction.method, assay.name, dims, cluster.col = 'RNA_snn_res.0.8'){
    Idents(sr) <- cluster.col
    sr_tree <- BuildClusterTree(sr, assay = assay.name, dims = dims, reduction = reduction.method)
    cluster_by <- FetchData(sr_tree, c(cluster.col, by))
    distinct_cluster_by <- distinct(cluster_by)
    tree <- Tool(sr_tree, 'BuildClusterTree')
    label <- tree$tip.label 
    for (i in 1:length(label)){
        index <- which(distinct_cluster_by[,cluster.col] == as.numeric(label[i]))
        label[i] <- paste(label[i], paste(distinct_cluster_by[index,2], collapse = "|"))
    }
    tree$tip.label <- label
    reso <- 600
    png(save_name,  width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
    plot(tree)
    dev.off()
    return(tree)
}