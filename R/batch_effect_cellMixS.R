# visualize batch effect
# sce : single cell experiment object (convert from sr with as.singlecellexperiment(sr))
# batch: colname of what consider batch effect.
# by default this will be set as date of sequencing
plotBatchVis <- function(sce, batch = "Date.of.Library", save_path, col ){
p <- visGroup(sce, group = batch) + scale_colour_manual(values = col)
save_plot(paste0(save_path, '/batch_vis_by_',
batch, '.png'), p)
}

# dim_met = 'LSI' for atac, 'PCA' for rna
# save_name can be batch correction method to compare before vs 
# after or between different methods
# neighbors can be calculated from the average number of cells per cluster
calculate_cms <- function(sce, neighbors = 200, save_name, dim_met){
    cms <- cms(sce, k =neighbors, group = 'Date.of.Library', res_name = paste0(save_name, '_date_k', neighbors), n_dim = 30, cell_min = 100, dim_red = dim_met)
    cms <- cms(cms, k =neighbors, group = 'library', res_name = paste0(save_name,'_lib_k',neighbors), n_dim = 30, cell_min = 100, dim_red = dim_met)
    cms <- cms(cms, k =neighbors, group = 'sampleID', res_name = paste0(save_name,'_sample_k',neighbors), n_dim = 30, cell_min = 100, dim_red = dim_met)
    cms <- cms(cms, k =neighbors, group = 'Individual.ID', res_name = paste0(save_name,'_patient_k',neighbors), n_dim = 30, cell_min = 100, dim_red = dim_met)
    return (cms)
}

calculate_n_plot_cms <- function(sce, save_path, save_name, neighbors, reduce_method){
    cms = calculate_cms(sce, neighbors = neighbors, save_name = save_name, reduce_method)
    plot_cms = visHist(cms)
    saveRDS(cms, file = paste0(save_path, '/', save_name, '_cms.RDS'))
    save_plot(paste0(save_path, '/', save_name, '_cms.png'), plot_cms)
    return(cms)
}