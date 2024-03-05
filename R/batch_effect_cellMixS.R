# visualize batch effect
# sce : single cell experiment object (convert from sr with as.singlecellexperiment(sr))
# batch: colname of what consider batch effect.
# by default this will be set as date of sequencing
plotBatchVis <- function(sce, batch = "Date.of.Library", save_path, col ){
p <- visGroup(sce, group = batch) + scale_fill_manual(values=col)
save_plot(paste0(save_path, '/batch_vis_by_',
batch, '.png'), p)
}