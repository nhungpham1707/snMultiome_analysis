# https://stuartlab.org/signac/articles/motif_vignette
atac_hm = readRDS('output/sc_atac/merge_all/atac_hm_tumor_nona.RDS')

loadd(hg38)
atacchr = createSrWChromatinAssay(atac_hm, hg38)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


atacchr <- AddMotifs(
  object = atacchr,
  genome = hg38,
  pfm = pfm
)
