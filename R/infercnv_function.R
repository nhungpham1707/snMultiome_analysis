# sR output from singleR
# sr is seurat atac with gene activity 
# make_anno_count_mx <- function(sR, sr, save_path){
#     # sr <- sR$sr
#     lib <- unique(sr$library)
#     # singleR_res = readRDS(paste0(sngR_path, '/singleR_', lib, '.RDS'))
#     labels <- sR$labels
#     names(labels) <- rownames(sR)
#     sr_sngr <- AddMetaData(sr, metadata=labels, col.name = 'singleR_labels')
#     count_matrix <- GetAssayData(sr, slot = 'counts')
#     saveRDS(count_matrix,paste0(save_path,"/", lib, "_sr_count_matrix.RDS") )
#     df <- data.frame(cell = colnames(sr_sngr), type = sr_sngr@meta.data$singleR_labels)
#     write.table(df, file = paste0(save_path,"/", lib, "_sr_cell_annotation.txt"), sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE)
# }

# mrgsr is merge seurat after gene activity
make_anno_count_Mrgmx <- function(mrgsr, save_path){
    count_matrix <- GetAssayData(mrgsr, slot = 'counts')
    saveRDS(count_matrix,paste0(save_path,"/merge_sr_count_matrix.RDS") )
    df <- data.frame(cell = colnames(mrgsr), type = mrgsr@meta.data$singleR_labels)
    write.table(df, file = paste0(save_path,"/merge_sr_cell_annotation.txt"), sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE)
}

# sr after singler annotation. for atac, sr with gene activity after singler annotation
make_anno_count_mx <- function(sr, save_path){
    lib <- unique(sr$library)
    count_matrix <- GetAssayData(sr, slot = 'counts')
    saveRDS(count_matrix,paste0(save_path,"/", lib, "_sr_count_matrix.RDS") )

    # get singlr label
    df <- data.frame(cell = colnames(sr), type = sr$singleR_labels)
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

# analyze infercnv output 
# from scFacility https://github.com/scgenomics/scgenomics-public.github.io/blob/main/docs/16-infercnv/16-infercnv.R
## trick by Jurrian de Kanter
.count_genes_in_runs <- function (x, length=50, SDs=2) {
  ## count total number of genes in stretches longer than LENGTH 
  ## exceeding SDs std.devs from mean. 
  runs <- rle(x > SDs)
  up <- sum(runs$lengths[ runs$lengths>length & runs$values])
  runs <- rle(x <   -SDs)
  down <- sum(runs$lengths[ runs$lengths>length & runs$values])
  ## (don't use abs(x) >  because pos and neg stretches might overlap, so partly cancel)
  sum(up,down)
} # .count_genes_in_runs

aneuploidy_score <- function(adj_expr,
                     length = 70,
                     SDs = 1.5) {
  apply(adj_expr, 2, .count_genes_in_runs, length=length, SDs=SDs)
}  ## aneuploidy_score

calculate_aneuploidy_ref_score <- function(infercnv_obj){
    refexp <- infercnv_obj@expr.data[ , unname(unlist(infercnv_obj@reference_grouped_cell_indices)) ]
    ref.scaled <- scale(refexp)
    refscore <-  aneuploidy_score(ref.scaled, SDs=SDs)

}
calculate_aneuploidy_obs_score <- function(infercnv_obj){
    obsexp <- infercnv_obj@expr.data[ , unname(unlist(infercnv_obj@observation_grouped_cell_indices)) ]

    ## now scale them so that SD is the 'unit' of Modif. Expression:
    obs.scaled <- scale(obsexp)
    obsscore <-  aneuploidy_score(obs.scaled, length=len)
}
plot_aneuploidy_score <- function(infercnv_obj){
    refscore <- calculate_aneuploidy_ref_score(infercnv_obj)
    obsscore <- calculate_aneuploidy_obs_score(infercnv_obj)
    ttest <- t.test(obsscore, refscore)
    plot(density(obsscore), 
      main=sprintf("Aneuploidy score length=%d SDs=%.1f",len, SDs),
      xlab=sprintf("t=%.2f", ttest$statistic))
    rug(obsscore, ticksize=0.01)
    rug(refscore, ticksize= -0.01, col='red')
    lines(density(refscore), col='red')
    legend(x="topright", legend=c('tumor', 'reference'), col=c('black','red'),
       pch=NULL, lty=1, bty='n')
}

add_aneuploidy_score_to_seurat <- function(infercnv_obj, srat, cutoff){
    refscore <- calculate_aneuploidy_ref_score(infercnv_obj)
    obsscore <- calculate_aneuploidy_obs_score(infercnv_obj)
    allscore <- c(obsscore,refscore)
    allscore <- allscore[ colnames(srat) ] # put in order of srat object.
    srat <- AddMetaData(srat, col.name='aneuploidy_score', allscore)

    aneuploid <- ifelse( srat@meta.data$aneuploidy_score < cutoff, 'no', 'yes')
    srat <- AddMetaData(srat, col.name='is_aneuploid', aneuploid) 
}






# from infercnv github
# https://github.com/broadinstitute/infercnv/blob/master/scripts/plot_tumor_vs_normal_chr_densities.R
plot_tumor_vs_normal_chr_densities <- function(infercnv_obj, save_path){
ref_group_cell_indices = infercnv:::get_reference_grouped_cell_indices(infercnv_obj)

chrs = unique(infercnv_obj@gene_order$chr)


for (chr in chrs) {
        
    gene_idx = which(infercnv_obj@gene_order$chr == chr)
    
    ref_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx,ref_group_cell_indices])
    
    df = data.frame(class='normal', vals=ref_data_pts)
    
    for (tumor in names(infercnv_obj@observation_grouped_cell_indices) ) {
        
        tumor_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ tumor ]]
        tumor_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx, tumor_cell_idx])
        
        df = rbind(df, data.frame(class=tumor, vals=tumor_data_pts))
    }

    flog.info(sprintf("Plotting data for chr: %s", chr))
    
    p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + ggtitle(chr) # + scale_y_continuous(trans='log10', limits=c(1,NA))
savePlot(paste0(save_path, '/tumor_vs_normal_chr_', chr, '.png'), p)
    } 
}


# https://github.com/broadinstitute/infercnv/blob/master/scripts/plot_tumor_vs_normal_chr_densities.i3.R
plot_tumor_vs_normal_chr_densities_sd <- function(infercnv_obj, save_path){
ref_group_cell_indices = infercnv:::get_reference_grouped_cell_indices(infercnv_obj)

normal_sd_trend = infercnv:::.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj)

mu = normal_sd_trend$mu
sigma = normal_sd_trend$sigma

chrs = unique(infercnv_obj@gene_order$chr)

delta = infercnv:::get_HoneyBADGER_setGexpDev(gexp.sd=sigma, alpha=0.05, k_cells=7)

for (chr in chrs) {
        
    gene_idx = which(infercnv_obj@gene_order$chr == chr)
    
    ref_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx,ref_group_cell_indices])
    
    df = data.frame(class='normal', vals=ref_data_pts)
    
    for (tumor in names(infercnv_obj@observation_grouped_cell_indices) ) {
        
        tumor_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ tumor ]]
        tumor_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx, tumor_cell_idx])
        
        df = rbind(df, data.frame(class=tumor, vals=tumor_data_pts))
    }

    flog.info(sprintf("Plotting data for chr: %s", chr))
    
    p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + ggtitle(chr) # + scale_y_continuous(trans='log10', limits=c(1,NA))
    
    
    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=mu,'sd'=sigma)) +
        stat_function(fun=dnorm, color='blue', args=list('mean'=mu-delta,'sd'=sigma)) +
        stat_function(fun=dnorm, color='blue', args=list('mean'=mu+delta,'sd'=sigma)) 
    
    savePlot(paste0(save_path, '/tumor_vs_normal_sd_chr_', chr, '.png'), p)

    }
}
