# sc_atacseq from our multiome data and online data (i.e. descartes) often have different peaks coordinates. These differences make it difficult (or impossible) to compare them. In this script, we generate a master peak list which contain peaks from both datasets (with modified coordinate to make them similar)
# Nhung et al 2024 


loadd(atac_hthymrgDim)
dsc_atac <- atac_hthymrgDim

loadd(atac_hm_w_tumor_label)
atac <- atac_hm_w_tumor_label
loadd(hg38)

dsc_atacchr <- createSrWChromatinAssay(dsc_atac, hg38)
atacchr <- createSrWChromatinAssay(atac, hg38)
atac_gr <- granges(atacchr)
dsc_gr <- granges(dsc_atacchr)
gr_atac <- atac_gr
gr_dsc <- dsc_gr

olap <- findOverlaps(gr_dsc, gr_atac)
sub_olap <- subsetByOverlaps(gr_dsc, gr_atac)
olap_count <- countOverlaps(gr_dsc, gr_atac )

# https://www.biostars.org/p/473162/
library(tidyverse)
library(valr)

# transform example data from above to data.frame and add required column `"chrom"`
gr_atac_df <-   as.data.frame(gr_atac) %>% 
  mutate(chrom = "chr1") %>% 
  dplyr::select(chrom, start, end) 
gr_dsc_df <- as.data.frame(gr_dsc) %>% 
  mutate(chrom = "chr1") %>% 
  dplyr::select(chrom, start, end) 

# get intersect e.g overlap  

overlap <- valr::bed_intersect(gr_atac_df, gr_dsc_df, suffix = c("_gr_atac", "_gr_dsc"))
overlap_dsc_atac <- valr::bed_intersect(gr_dsc_df, gr_atac_df, suffix = c("_dsc", "_atac"))
overlap_dsc_atac # somehow the coordinate here does not match in atac 
dim(distinct(overlap_dsc_atac[,1:2]))
a <- distinct(overlap[,1:2])
dim(a)
dim(overlap)
# bed_glyph(bed_intersect(gr_atac[1:100,], gr[1:100,]))

# 
count_mtx_dsc <- GetAssayData(dsc_atac)
rownames(dsc_atac)[1]
chr1 <- count_mtx_dsc[grep('chr1-9992-10688', rownames(count_mtx_dsc)),]
length(chr1)
chr1_no0 <- chr1[chr1>0]
max(chr1)
min(chr1)
max(count_mtx_dsc)
count_atac <- atac@assays$peaks@counts
max(count_atac)
min(count_atac)
mrg <- merge(x= atac, y = dsc_atacchr)


## only work with accessiblepeaks ---
pe <- AccessiblePeaks(atac)




# 

loadd(atac_hthymrgDim)
dsc_atac <- atac_hthymrgDim

loadd(atac_hm_w_tumor_label)
atac <- atac_hm_w_tumor_label
loadd(hg38)

dsc_atacchr <- createSrWChromatinAssay(dsc_atac, hg38)
atacchr <- createSrWChromatinAssay(atac, hg38)
atac_gr <- granges(atacchr)
dsc_gr <- granges(dsc_atacchr)
gr_atac <- atac_gr
gr_dsc <- dsc_gr
count <- GetAssayData(atac)
new_atac_count_mx <- makeCountMx_withSamePeaks(gr_dsc, gr_atac, count)




makeCountMx_withSamePeaks <- function(querry_gr, hit_gr, count_hit){
olap <- findOverlaps(querry_gr,hit_gr)
olap_df <- as.data.frame(olap) # 568396 peaks dsc vs atac_nohm 
unique_olap_querry <- unique(olap_df[,1]) # 262577 peaks 

new_count_mx <- c()
names <- c()
for (i in 1:length(unique_olap_querry)){
message(paste('running', i, 'peaks out of', length(unique_olap_querry), 'unique overlap peaks'))
querry_index <- unique_olap_querry[i]
hit_index <- olap_df[olap_df[,1]==querry_index,2]

# change coordinate to that in querry (dsc)
new_coor <- paste0(querry_gr[querry_index]@seqnames,'-', querry_gr[querry_index]@ranges)

if (length(hit_index) > 1){
  new_value <- colSums(count_hit[hit_index,])
} else {
  new_value <- count_hit[hit_index,]
}
new_count_mx <- rbind(new_count_mx, new_value)
names <- rbind(names, new_coor)
rownames(new_count_mx) <- names[,1]
}
return (new_count_mx)
}

olap_querry <- unique_olap_querry[i]

get_new_coor_value <- function(olap_querry, querry_gr, olap_df, count_hit){
querry_index <- olap_querry
hit_index <- olap_df[olap_df[,1]==querry_index,2]
new_count_mx <- c()
name <- c()
# change coordinate to that in querry (dsc)
new_coor <- paste0(querry_gr[querry_index]@seqnames,'-', querry_gr[querry_index]@ranges)

if (length(hit_index) > 1){
  new_value <- colSums(count_hit[hit_index,])
} else {
  new_value <- count_hit[hit_index,]
}
new_count_mx <- rbind(new_count_mx, new_value)
names <- rbind(names, new_coor)
rownames(new_count_mx) <- names[,1]
}
sapply(unique_olap_querry, get_new_coor_value, querry_gr = dsc_gr, olap_df = olap_df, count_hit = count_hit)



# to speed up the process 
# remvoe all 0 rows in ocuntmatrix (not sure why they are still there)

freq_0 <- rowSums(atac_hm_count !=0)
# keep only features with at least 40 cells non 0 
atac_non0 <- atac_hm_count[freq_0 > 100,] # reduce to 972263 features from 1760372 in atac_hm 
feature_tokeep <- rownames(atac_hm_count)[freq_0 >100]

atac <- readRDS('output/sc_atac/merge_all/atac_hm_tumor_nona.RDS')


atac_sub <- subset(atac, features = feature_tokeep, slot = 'count') # does not work, not sure how to define peaks assay 
# remove features in atac_gr 
atac_gr_df <- as.data.frame(atac_gr)
atac_gr_df$name <- paste0(atac_gr_df$seqnames, '-', atac_gr_df$start, '-', atac_gr_df$end)
index_tokeep <- atac_gr_df$name %in% feature_tokeep
new_atac_gr <- atac_gr[index_tokeep]



id <- seq(1,length(new_atac_gr), by = 1000)
subset_atac <- function(id, new_atac_gr, atac_non0){
  
}