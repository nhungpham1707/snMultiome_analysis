# sR output from singleR
# sr is seurat atac with gene activity 
make_anno_count_mx <- function(sR, sr, save_path){
    # sr <- sR$sr
    lib <- unique(sr$library)
    # singleR_res = readRDS(paste0(sngR_path, '/singleR_', lib, '.RDS'))
    labels <- sR$labels
    names(labels) <- rownames(sR)
    sr_sngr <- AddMetaData(sr, metadata=labels, col.name = 'singleR_labels')
    count_matrix <- GetAssayData(sr, slot = 'counts')
    saveRDS(count_matrix,paste0(save_path,"/", lib, "_sr_count_matrix.RDS") )
    df <- data.frame(cell = colnames(sr_sngr), type = sr_sngr@meta.data$singleR_labels)
    write.table(df, file = paste0(save_path,"/", lib, "_sr_cell_annotation.txt"), sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE)
}

make_infercnvObj <- function(lib, normal_cells,geneOderLink, inLink){
    cellAnnoLink <- paste0(inLink,"/", lib, "_sr_cell_annotation.txt")
    countMx <- readRDS(paste0(inLink,"/", lib, "_sr_count_matrix.RDS"))

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix= countMx,
    annotations_file= cellAnnoLink,
    gene_order_file= geneOderLink,
    ref_group_names= normal_cells)
    res <- c()
    res[['lib']] <- lib
    res[['infObj']] <- infercnv_obj
    return(res)
}

run_infercnv <- function(res, outLink){
    lib <- res$lib
    outLbLink <- paste0(outLink,'/', lib)
    dir.create(outLbLink, recursive = TRUE)
    infercnv_obj <- res$infObj
    infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.3, 
                             out_dir= outLbLink, 
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=F,
                             analysis_mode = "samples",
                             num_threads = 3,
                             output_format = "pdf",
                             window_length = 101,
                             save_final_rds = T,
                             plot_steps=F,
                             sd_amplifier = 2)

message('finish infercnv without error')
}