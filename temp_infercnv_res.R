
refexp <- infercnv_obj@expr.data[ , unname(unlist(infercnv_obj@reference_grouped_cell_indices)) ]

obsexp <- infercnv_obj@expr.data[ , unname(unlist(infercnv_obj@observation_grouped_cell_indices)) ]

## now scale them so that SD is the 'unit' of Modif. Expression:
obs.scaled <- scale(obsexp)
ref.scaled <- scale(refexp)

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




library(infercnv)
library(ggplot2)
library(futile.logger)
library(dplyr)


infercnv_obj = readRDS(infercnv_obj_file)

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