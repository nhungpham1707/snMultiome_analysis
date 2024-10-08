# slot_name is the field name in the sr. e.g. "RNA", 'Gene Expression', 'peaks'. For sc_atac with gene activity it is 'RNA'. for sc_rna, it is 'rna'
run_singleR <- function(sr, slot_name = 'RNA'){
    input <- as.matrix(sr[[slot_name]]@data)
    hpca.se <- celldex::HumanPrimaryCellAtlasData()
    labels = hpca.se$label.main
    singler <- SingleR(input,
                ref=hpca.se,
                labels=labels,
                method = NULL,
                clusters = NULL,
                genes = "de",
                sd.thresh = 1,
                de.method = "classic",
                de.n = NULL,
                de.args = list(),
                aggr.ref = FALSE,
                aggr.args = list(),
                recompute = TRUE,
                restrict = NULL,
                quantile = 0.8,
                fine.tune = TRUE,
                tune.thresh = 0.05,
                prune = TRUE,
                assay.type.test = "logcounts",
                assay.type.ref = "logcounts",
                check.missing = TRUE)
    return(singler)
}

plot_singler <- function(singler, sr, save_path){
    cluster <- sr@active.ident
    lib <- unique(sr$library)
    png(paste0(save_path, '/singler_score_heatmap_',lib, '.png'), width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
    plotScoreHeatmap(
    singler,
    cells.use = NULL,
    labels.use = NULL,
    clusters = cluster, #or NULL
    show.labels = TRUE,
    show.pruned = FALSE,
    max.labels = 40,
    normalize = TRUE,
    cells.order = NULL,
    order.by = "clusters",
    scores.use = NULL,
    calls.use = 0,
    na.color = "gray30",
    cluster_cols = FALSE,
    annotation_col = NULL, #ann_colors gave error
    show_colnames = FALSE,
    color = (grDevices::colorRampPalette(c("#D1147E", "white", "#00A44B")))(100),
    silent = FALSE,
    grid.vars = list()
    )
    dev.off()
}
# sR singler res
# sr seurat after normalize and dim reduc
get_sgR_label <- function(sR,sr){
    labels <- sR$labels
    names(labels) <- rownames(sR)
    # index <- which(colnames(sr@meta.data) == 'singleR_labels')
    # if (length(index) == 0 ){
    #     colname <- 'singleR_labels'
    # } else {
    #     colname <- paste0('singleR_labels_', save_name)
    # }
    sr_sngr <- AddMetaData(sr, metadata=labels, col.name = 'singleR_labels')
}

plot_mrg_singler <- function(singler, sr, save_path, save_name){
    cluster <- sr@active.ident
    lib <- save_name
    png(paste0(save_path, '/singler_score_heatmap_',lib, '.png'), width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
    plotScoreHeatmap(
    singler,
    cells.use = NULL,
    labels.use = NULL,
    clusters = cluster, #or NULL
    show.labels = TRUE,
    show.pruned = FALSE,
    max.labels = 40,
    normalize = TRUE,
    cells.order = NULL,
    order.by = "clusters",
    scores.use = NULL,
    calls.use = 0,
    na.color = "gray30",
    cluster_cols = FALSE,
    annotation_col = NULL, #ann_colors gave error
    show_colnames = FALSE,
    color = (grDevices::colorRampPalette(c("#D1147E", "white", "#00A44B")))(100),
    silent = FALSE,
    grid.vars = list()
    )
    dev.off()
}

# # group_singleR_labels <- function(sr){
#     malignant_cells <- setdiff(sr$singleR_labels, nonMalignant_cells)
#     sr$clean_sgr_labels <- sr$singleR_labels
#     sr$clean_sgr_labels[sr$clean_sgr_labels %in% malignant_cells] <- 'maybe_tumors'
#     sr$clean_sgr_labels[sr$clean_sgr_labels %in% immune_cells] <- 'immune_cells'
#     sr$clean_sgr_labels[sr$clean_sgr_labels]
#     sr$clean_sgr_labels[sr$clean_sgr_labels %in% nonMalignant_cells] <- 'others'
#     return(sr)
# }

group_singleR_labels <- function(sr){
    sr$group_sgr_labels <- sr$singleR_labels
    sr$group_sgr_labels[sr$group_sgr_labels %in% brain_cells] <- 'brain_cells'
    sr$group_sgr_labels[sr$group_sgr_labels %in% immune_cells] <- 'immune_cells'
    sr$group_sgr_labels[sr$group_sgr_labels %in% stem_cells] <- 'stem_cells'
    sr$group_sgr_labels[sr$group_sgr_labels %in% bone_cells] <- 'bone_cells'
    sr$group_sgr_labels[sr$group_sgr_labels %in% stroma_cells] <- 'stroma_cells'
    return(sr)
}