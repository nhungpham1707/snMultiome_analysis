# ref 
# https://github.com/plger/scDblFinder
# https://plger.github.io/scDblFinder/articles/scATAC.html

library(scDblFinder)
library(drake)
# loadd(atac.sce)
loadd(gexMito_all_LX049)
lx49 <- make_sce(gexMito_all_LX049)
soc <- read.csv()
known_db <- soc$status == 'doublet'
sce = scDblFinder(lx49, knownDoublets = known_db)
# sce <- scDblFinder(atac.sce, samples = "sampleID", 
#         BPPARAM=MulticoreParam(2))