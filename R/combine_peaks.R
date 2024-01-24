makeGrFile <- function(metadata,lb){
  link <- unique(metadata$data_link[which(metadata$name == lb)])
  bedFile <- read.table(paste0(base_data_dir, '/',link, 
                               "/outs/atac_peaks.bed"), 
                        col.names = c("chr", "start", "end"))
  grFile <- makeGRangesFromDataFrame(bedFile)
}

filterBadPeaks <- function(combined.peaks){
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
  return (combined.peaks)
}

filterChrBlacklist <- function(combined.peaks){
  ##Here filter further to only keep standard chromosomes
  combined.peaks <- keepStandardChromosomes(combined.peaks, 
                                            pruning.mode="coarse")
  ##Here filter to remove blacklisted genes
  combined.peaks <- subsetByOverlaps(x = combined.peaks, 
                                     ranges = blacklist_hg38_unified, 
                                     invert = TRUE)
  return (combined.peaks)
}