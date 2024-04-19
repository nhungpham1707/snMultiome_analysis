
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







infercnv_obj = readRDS(infercnv_obj_file)

