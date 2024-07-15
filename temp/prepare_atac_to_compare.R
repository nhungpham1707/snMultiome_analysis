setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/')

atac_hthy = readRDS('output/healthy_data/dsc_atac_hthymrgDim.RDS')
atac_hm = readRDS('output/sc_atac/merge_all/atac_hm_tumor_nona.RDS')
hg38_annotation = readRDS('output/hg38.RDS')
dsc_atacchr = createSrWChromatinAssay(atac_hthy, hg38_annotation)
atacchr = createSrWChromatinAssay(atac_hm, hg38_annotation)
atac_gr = granges(atacchr)
dsc_gr = granges(dsc_atacchr)
atac_hm_count = GetAssayData(atac_hm)
# reduce data size to speed up the make matrix func
# remove features with most count 0 
freq_0 = rowSums(atac_hm_count !=0)
# keep only features with at least 100 cells non 0 
atac_non0 = atac_hm_count[freq_0 > 100,] # reduce to 972263 features from 1760372 in atac_hm 
feature_tokeep = rownames(atac_hm_count)[freq_0 >100]
atac_gr_df = as.data.frame(atac_gr)
index_tokeep = paste0(atac_gr_df$seqnames, '-', atac_gr_df$start, '-', atac_gr_df$end) %in% feature_tokeep
new_atac_gr = atac_gr[index_tokeep]
id = seq(1,length(new_atac_gr), by = 10000)
id2 = id[2:length(id)]


l <- cached()
mx_ls <- l[grep('new_atachm_mx_[0-9]+', l)]
mx_ls
length(mx_ls)

loadd(mx_ls)

ncol(get(mx_ls[1]))

add_colname <- function(mx){
  mx <- get(mx)
  colnames(mx) <- colnames(atac_non0)
  return(mx)
}

mx_ls2 <- lapply(mx_ls, add_colname)
mrg_mx <- merge_sparseMx(mx_ls2)
mx_ls3 <- mx_ls2

mx1 <- mx_ls3[[1]]
mx_ls3[[1]] <- NULL
length(mx_ls3)

mrg_mx <- RowMergeSparseMatrices(mx1, mx_ls3)
dim(mrg_mx)
tail(colnames(mrg_mx))
length(colnames(mrg_mx))
length(intersect(colnames(mx1), colnames(mx_ls3[[1]])))
mrg_mx <- bind_rows(mx_ls2)
mrg_mx[1:10,1:10]
dim(mrg_mx)
dim(new_atachm_mx_100001)
rownames(mrg_mx)[1:10]
rownames(mx_ls2[[1]])[1:10]

mrg_mx3 <- bind_rows(mx_ls2, .id = NULL)

dim(atac_non0)

head(rownames(mrg_mx3))
tail(mx_ls)
f_inmrg <- rownames(mrg_mx3)
length(unique(f_inmrg))
f_inmrg2 <- gsub('[0-9]+', '',f_inmrg)
length(unique(f_inmrg2))
f_inmrg2
