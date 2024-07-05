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

olap <- findOverlaps(gr_dsc, gr_atac)

olap_df <- as.data.frame(olap) # 568396 peaks dsc vs atac_nohm 
unique_olap_dsc <- unique(olap_df[,1]) # 262577 peaks 

new_count_mx <- c()
names <- c()
for (i in 1:2){

querry_index <- unique_olap_dsc[i]
hit_index <- olap_df[olap_df[,1]==querry_index,2]


# change coordinate to that in querry (dsc)
new_coor <- paste0(gr_dsc[querry_index]@seqnames,'-', gr_dsc[querry_index]@ranges)


# rownames(count)[hit_index[1]] <- new_coor 
# atac@assays$peaks@counts@Dimnames[[1]][hit_index[1]] <- new_coor
# atac@assays$peaks@data@Dimnames[[1]][hit_index[1]] <- new_coor
# rownames(atac@assays$peaks@meta.features)[hit_index[1]] <- new_coor
# change count value by summing up all peaks in the region of querry peak
if (length(hit_index) > 1){
  new_value <- colSums(count[hit_index,])
} else {
  new_value <- count[hit_index,]
}
# count[hit_index[1]] <- new_value # this take forever ??? 

new_count_mx <- rbind(new_count_mx, new_value)
names <- rbind(names, new_coor)

}