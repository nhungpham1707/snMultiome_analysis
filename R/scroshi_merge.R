# annotate cell types using modified scROSHI functions (see more detail in scROSHI_modification.R)
# input: 
# - either rna or atac with gene activity RDS
# - marker list file 
# - config file to indicate major celltype and subtype, comma separated
# - example data: library(scROSHI)
# data('config')
# data('test_sce_data')
# data('marker_list')
# Output: singlecellExperiment object with extra columns for cell type in metadata

annotate_mrg_w_scROSHI <- function(sr, marker_list, config, pt = 1, cols, save_name, path_to_save){
    message(paste('read sr to annotate with scROSHI for', save_name))
    # sr <- readRDS(srat_filename)
    sr.sce <- as.SingleCellExperiment(sr)
    names(sr.sce@assays)
    rowData(sr.sce)$SYMBOL <- rownames(sr.sce) 
    sce_data = sr.sce
    count_data = "counts"
    message(paste('-----------start scROSHI for', save_name, '-----------'))
    results <- scROSHI_nhung(sce_data = sce_data,
                         celltype_lists = marker_list,
                         type_config = config,
                         count_data = 'counts')
    results_sr <- as.Seurat(results)
    # scroshi_dir <- dir.create(paste0(cell_type_dir/'scROSHI'), recursive = TRUE)
    saveRDS(results_sr, paste0(path_to_save, '/scROSHI_cell_type_annotation_', save_name,'.RDS'))
    # path_to_save <- dir.create(paste0(cell_type_fig_dir/'scROSHI'), recursive = TRUE)

    DimPlot(results_sr, group.by = 'celltype_final', cols = cols, pt.size = pt, raster = FALSE, raster.dpi = c(512, 512))
    ggsave(file = paste0(path_to_save, '/', save_name,"scROSHI_cell_type_annotation.png"), 
        width = 1200 * reso/72, height = 700 * reso/72, units ="px", dpi = reso)
    summary_res <- table(results_sr@meta.data$celltype_final)
    write.csv(summary_res, paste0(path_to_save, '/scROSHI_cell_annotation_', save_name, '.csv'))
    return(results_sr)
}

run_scROSHI_mrg_w_demo_data <- function(sr, cols, pt = 2, save_name, path_to_save){
    data('config')
    data('marker_list')
    results_sr <- annotate_mrg_w_scROSHI(sr, marker_list, config, pt = pt, cols, save_name, path_to_save)
    return(results_sr)
    message('finish scROSHI with demo data no error!')
}

run_scROSHI_mrg_w_c8_data <- function(sr, cols, pt = 1, save_name, path_to_save){
    marker_list <- readRDS(paste0(base_data_dir, '/c8_marker_list.RDS'))
    # marker_list <- marker_list[grep('FETAL', names(marker_list))] 
    marker_list <- marker_list[grep('DESCARTES_FETAL_ADRENAL', names(marker_list))] 
    # config_list <- readRDS(paste0(base_data_dir, '/c8_config.RDS'))
    cell_type <- names(marker_list)
    config_list <- data.frame(Major = cell_type,
                        Subtype = c(rep('none', time = length(cell_type))))
    annotate_mrg_w_scROSHI(sr, marker_list, config, pt = 1, cols, save_name, path_to_save)
    message('finish scROSHI w c8 with no error')
}


run_scROSHI_mrg_w_atrt_data <- function(sr, cols, pt = 1, save_name, path_to_save){
    ATRT_TYR <- c('MITF', 'OTX2', 'TYR', 'PDGFRB', 'JAK1', 'BMP4')
    ATRT_SHH <- c('NOTCH1', 'GLI2', 'MYCN', 'ASCL1', 'HES1', 'DTX1', 'PTCH1', 'BOC')
    ATRT_MYC <- c('HOXC10', 'CCND3', 'MYC')

    ATRT_list <- c()
    ATRT_list[['ATRTTYR']] <- ATRT_TYR
    ATRT_list[['ATRTSHH']] <- ATRT_SHH
    ATRT_list[['ATRTMYC']] <- ATRT_MYC
    ATRT_config <- data.frame(Major = names(ATRT_list),
                            Subtype = c(rep('none', time = length(names(ATRT_list)))))
    marker_list <- ATRT_list
    config <- ATRT_config
    results_sr <- annotate_mrg_w_scROSHI(sr, marker_list, config, pt = 1, cols, save_name, path_to_save)
    message('finish scROSHI w atrt with no error!')
    return(results_sr)
}


run_scROSHI_for_all_samples <- function(metadata, data_type = 'sc_atac'){
 message('---start running annotation w scROSHI------')
    library <- unique(metadata$name)
    if (data_type == 'sc_atac'){
        path_to_save <- cell_type_atac_scroshi_dir
    } else {
        path_to_save <- cell_type_rna_scroshi_dir
    }
    for (i in 1: length(library)) {
        message( paste('----running scROSHI annotation for', i, 'out of', length(library), 'samples -----'))
        library_name <- library[i]
        save_name <- paste0(data_type, '_',library_name)
        res_name <- paste0(path_to_save, '/scROSHI_cell_type_annotation_', save_name,'.RDS')
        if (file.exists(res_name)){
        message (paste(library_name, data_type,'already ran scroshi! skip!'))
        } else {
            if (data_type == 'sc_atac'){
                if (metadata$souporcell_link[which(metadata$name == library_name)[1]] == "") {
                        extension_name <- '_scATAC_isoutlier_normalize_gene_activity.RDS'
                } else {
                        extension_name <- '_scATAC_after_demultiplex_isoutlier_normalize_gene_activity.RDS' }
                srat_filename <- paste0(sc_atac_processing_dir, '/', library_name,extension_name )
                # slot_name = "RNA"
                
            } else {
                srat_filename <- paste0(sc_rna_processing_dir, '/',library_name,'_normalize_isoutlier_filter_rna_data_w_features_scale_pca_umap_cluster.RDS' )
                # slot_name = 'rna'
                }
            sr <- readRDS(srat_filename)

            run_scROSHI_w_demo_data(sr, cols, pt = 2, save_name = save_name , path_to_save)
    }
    message('----finish annotation w scroshi for all samples----')
    }
}
