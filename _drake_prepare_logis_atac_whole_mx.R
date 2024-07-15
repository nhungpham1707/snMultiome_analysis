setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

atac_mx_plan <- drake_plan(
 atac_hthy = readRDS('output/healthy_data/dsc_atac_hthymrgDim.RDS'),
  atac_hm = readRDS('output/sc_atac/merge_all/atac_hm_tumor_nona.RDS'),
  hg38_annotation = readRDS('output/hg38.RDS'),
 dsc_atacchr = createSrWChromatinAssay(atac_hthy, hg38_annotation),
  atacchr = createSrWChromatinAssay(atac_hm, hg38_annotation),
  atac_gr = granges(atacchr),
  dsc_gr = granges(dsc_atacchr),
  atac_hm_count = GetAssayData(atac_hm),
  # reduce data size to speed up the make matrix func
  # remove features with most count 0 
  freq_0 = rowSums(atac_hm_count !=0),
  # keep only features with at least 100 cells non 0 
  atac_non0 = atac_hm_count[freq_0 > 300,], # reduce to 972263 features from 1760372 in atac_hm 
  feature_tokeep = rownames(atac_hm_count)[freq_0 >300],
  atac_gr_df = as.data.frame(atac_gr),
  index_tokeep = paste0(atac_gr_df$seqnames, '-', atac_gr_df$start, '-', atac_gr_df$end) %in% feature_tokeep,
  new_atac_gr = atac_gr[index_tokeep],
  new_atachm_mx = makeCountMx_withSamePeaks_optimized3(dsc_gr,new_atac_gr, atac_non0)
  )

make(atac_mx_plan)
