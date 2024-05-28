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

db <- scDblFinder(atac.sce, samples = "library")

# 
rna <- rna_nodb_infer
rna.sce <- make.sce(rna)
db <- scDblFinder(rna.sce, samples = 'library') # found 8229 cells as doublet, remove but still not satisfy 


# run again with increase expected db rate 

db_0.5 <- scDblFinder(rna.sce, samples = 'library', dbr = 0.5)
# find 21957 cells as doublets 

# check if suspected doublet in cluster 6 res 0.8, p3W rms is in doubletfinder 
suspect_index <- rna$library == 'LX379_LX380_an_596' & rna$RNA_snn_res.0.8 == 6 & rna$Subtype == 'FP-RMS (P3F)'

suspect_cells <- colnames(rna)[suspect_index]
db_res <- colnames(db_0.5_sr)[db_0.5_sr$scDblFinder.class == 'doublet']

# 285 cells are found as doublet from scDblFinder dr 0.5 
length(intersect(db_res, suspect_cells))

DimPlot(db_0.5_sr, cells.highlight= intersect(db_res, suspect_cells))

db_0.2 <- scDblFinder(rna.sce, samples = 'library', dbr = 0.2)