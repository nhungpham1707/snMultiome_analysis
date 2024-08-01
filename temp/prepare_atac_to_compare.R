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

# 22 07 2024 ----
source('_drake.R')
loadd(new_atachm_mx)
loadd(atac_hm_tumor_nona)
loadd(train_dsc_atac40k)
new_atachm_mx[1:10,1:10]
colnames(new_atachm_mx) <- colnames(atac_hm_tumor_nona)
p_atac <- predictSimilarity(train_dsc_atac40k, new_atachm_mx, 
                            classes = atac_hm_tumor_nona$cell_identity,
                            minGeneMatch = 0.0,
                            logits = F)

# new idea - extract only overlap feature in dsc and train model based on that 

loadd(atac_hthymrgDim)

dsc_markers <- FindAllMarkers(atac_hthymrgDim)
atac_markers <- FindAllMarkers(atac_hm_tumor_nona)
# keep only top 5k features for each cell types 
topfeatures_dsc <- dsc_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5000, 
        wt = avg_log2FC)

features_to_keep_dsc <- topfeatures_dsc$gene
sub_dsc <- subset(atac_hthymrgDim, features = features_to_keep_dsc)

loadd(atac_markers)
topfeatures_atac <- atac_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5000, 
        wt = avg_log2FC)
features_to_keep_atac <- topfeatures_atac$gene
sub_atac <- subset(atac_hm_tumor_nona, features = features_to_keep_atac) # error: no RNA assay --> why? it is the same with dsc, but dsc work ????

loadd(hg38)
dsc_atacchr = createSrWChromatinAssay(sub_dsc, hg38)
atacchr = createSrWChromatinAssay(sub_atac, hg38)
atac_gr = granges(atacchr)
dsc_gr = granges(dsc_atacchr)

new_atachm_mx = makeCountMx_withSamePeaks_optimized3(dsc_gr,atac_gr, sub_atac)

loadd(atac_hm_tumor_nona)
loadd(atac_markers)
loadd(atac_gr)
loadd(dsc_gr)
loadd(atac_hthymrgDim)
sub_atac_mx <- extract_atac_w_n_features(n = 20000, atac_markers,GetAssayData(atac_hm_tumor_nona), atac_gr, dsc_gr, atac_hm_tumor_nona)

sub_atac_mx <- extract_atac_w_n_features(n = 2000, atac_markers,GetAssayData(atac_hm_tumor_nona), atac_gr, dsc_gr, atac_hm_tumor_nona)

sub_atac_mx <- extract_atac_w_n_features(n = 100, atac_markers,GetAssayData(atac_hm_tumor_nona), atac_gr, dsc_gr, atac_hm_tumor_nona)

sub_atac_mx <- extract_atac_w_n_features(n = 10000, atac_markers,GetAssayData(atac_hm_tumor_nona), atac_gr, dsc_gr, atac_hm_tumor_nona)
saveRDS(file = 'output/logistic_regression/atac_all_markers.RDS', sub_atac_mx)
sub_dsc <- subset(atac_hthymrgDim, features = rownames(sub_atac_mx))

train_sub_dsc <- trainModel(GetAssayData(sub_dsc), classes = sub_dsc$cell_type, maxCells = 4000)
saveRDS(file = 'output/logistic_regression/train_sub_dsc_w_overlap_atac_markers.RDS', train_sub_dsc)

p_sub_dsc_atac <- predictSimilarity(train_sub_dsc, sub_atac_mx, 
                            classes = atac_hm_tumor_nona$cell_identity,
                            minGeneMatch = 0.9,
                            logits = F)
similarityHeatmap(p_sub_dsc_atac)

# 01 08 2024 ----
# try feature selection to reduce the number of feature 
library(bestglm)
library(caret)
library(pROC)
library(glmnet)
library(InformationValue) # no package to install ??
loadd(new_atachm_mx)
loadd(atac_hthymrgDim)

dsc_only_overlap_countMx <- subset(atac_hthymrgDim, features = rownames(new_atachm_mx))
dsc_overlap_mx <- as.matrix(dsc_only_overlap_countMx)
saveRDS(file = 'output/logistic_regression/dsc_only_overlap_countMx.RDS', dsc_overlap_mx)  # take  a long time 

dsc_atac <- readRDS('output/logistic_regression/dsc_only_overlap_countMx.RDS')

dsc_atac$cell_type <- atac_hthymrgDim$cell_type

