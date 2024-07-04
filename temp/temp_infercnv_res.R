
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

len <- 50 # also try e.g. 30 or 70
SDs <- 1.3 # also try e.g. 1.5 or 2
obsscore <-  aneuploidy_score(obs.scaled, length=len)
refscore <-  aneuploidy_score(ref.scaled, SDs=SDs)

ttest <- t.test(obsscore, refscore)

p <- plot(density(obsscore), 
      main=sprintf("Aneuploidy score length=%d SDs=%.1f",len, SDs),
      xlab=sprintf("t=%.2f", ttest$statistic)) 
rug(obsscore, ticksize=0.01) 
rug(refscore, ticksize= -0.01, col='red') 
lines(density(refscore), col='red') 
legend(x="topright", legend=c('tumor', 'reference'), col=c('black','red'),
       pch=NULL, lty=1, bty='n')

savePlot('lx49_aneuploid_score.png',p)
