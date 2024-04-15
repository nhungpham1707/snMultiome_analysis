# ?RunHarmony

# theta Diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.

#   sigma Width of soft kmeans clusters. Default sigma=0.1.Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering.
# lambda Ridge regression penalty. Default lambda=1. 
# Bigger values protect against over correction. If 
# several covariates are specified, then lambda can also 
# be a vector which needs to be equal length with the 
# number of variables to be corrected. In this scenario, 
# each covariate level group will be assigned the scalars
# specified by the user. If set to NULL, harmony will start lambda estimation mode to determine lambdas automatically and try to minimize overcorrection (Use with caution still in beta testing).

harmony_n_plot <- function(sr, batch_factor, reduction = 'pca', assay = 'RNA', theta = 0, sigma = 0.1, lambda = 1,tau = 0, save_path){
    save_name <- paste0('assay_', assay, 'batch_factor_', batch_factor, 'theta_', theta, '_sigma_', sigma, '_tau_', tau, '_lambda_', lambda)
    message(paste('run harmony for', save_name))
    message('run harmony -------')
    hm = RunHarmony(sr, group.by.vars = batch_factor, reduction.use = reduction, project.dim = FALSE, theta = theta, sigma = sigma, assay.use = assay, lambda = lambda, .options = harmony_options(tau = tau))
    message('run find neighbor--------')
    hm = FindNeighbors(object = hm, reduction = "harmony", k.param = 30)
    message('run find cluster-------')
    hm = FindClusters(hm, resolution = c(0.2,0.4,0.6, 0.8,1))
    message('run umap----------')
    hm = RunUMAP(hm, dims = 1:30, reduction = 'harmony')
    message('dimplot---------')
    hm_p = DimPlot(hm, group.by = 'library', raster = FALSE, pt.size = 0.1, cols = my_cols)

    savePlot(paste0(save_path,'/', save_name, '_lib.png'), hm_p)

    hm_type_p = DimPlot(hm, group.by = 'Subtype', raster = FALSE, pt.size = 0.1, cols = my_cols)

    savePlot(paste0(save_path,'/', save_name,'_type.png'), hm_type_p)

    hm_sgr_p = DimPlot(hm, group.by = 'singleR_labels', raster = FALSE, pt.size = 0.1, cols = my_cols)

    savePlot(paste0(save_path,'/', save_name,'_sgr.png'), hm_sgr_p)

    hm_treatment_p = DimPlot(hm, group.by = 'treatment', raster = FALSE, pt.size = 0.1)

    savePlot(paste0(save_path,'/', save_name,'_treatment.png'), hm_treatment_p)

    return(hm)
 }

 