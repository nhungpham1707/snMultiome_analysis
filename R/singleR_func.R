# slot_name is the field name in the sr. e.g. "RNA", 'Gene Expression', 'peaks'. For sc_atac with gene activity it is 'RNA'. for sc_rna, it is 'rna'
run_singleR <- function(sr, save_path, slot_name = 'RNA'){
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
    lib <- unique(sr$library)
    # unaligned.idents <- sr@active.ident
    saveRDS(singler, file = paste0(save_path, '/singleR_', lib, '.RDS'))
    res <- c()
    res[['sr']] <- sr
    res[['singler']] <- singler
    res[['cluster']] <- sr@active.ident
    res[['lb']] <- lib
    return(res)
}

plot_singler <- function(res, save_path){
    lib <- res[['lb']]
    singler <- res[['singler']]
    cluster <- res[['cluster']]
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