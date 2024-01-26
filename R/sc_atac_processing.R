getHg38Annotation <- function(){
  hg38 <- EnsDb.Hsapiens.v86
  seqlevelsStyle(hg38) <- 'UCSC'
  annotation <- GetGRangesFromEnsDb(hg38)
  return(annotation)
}

getFragPath <- function(metadata, lb){
  index <- which(metadata$name == lb)[1]
  fragpath <- paste0(base_data_dir, '/', 
                     metadata$data_link[index], 
                     "/outs/atac_fragments.tsv.gz")
  return(fragpath)
}
create_atacAssay <- function(lb, metadata, annotation){
  index <- which(metadata$name == lb)[1]
  fragpath <- getFragPath(metadata, lb)
  counts <- Read10X_h5(paste0(base_data_dir, '/', 
                             metadata$data_link[index], 
                             '/outs/filtered_feature_bc_matrix.h5'))
  peak_counts <- counts[["Peaks"]]
  
  atacAssay <- CreateChromatinAssay(
    counts = peak_counts,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation,
    min.cells = 10,
    min.features = 200)
  return(atacAssay)
}

create_atacDisjoin <- function(atacAssay, lb, metadata, 
                           combined.peaks, annotation){
  fragpath <- getFragPath(metadata, lb)
  frags <- CreateFragmentObject(path = fragpath)
  CRcounts <- FeatureMatrix(fragments = Fragments(atacAssay),
                            features = combined.peaks, 
                            cells = colnames(atacAssay))
  peakAssay <- CreateChromatinAssay(
    counts = CRcounts,
    fragments = fragpath,
    annotation = annotation)
  atacSr <- CreateSeuratObject(counts = peakAssay,
                                assay = "peaks" )
  atacSr$library <- lb
  atacSr$barcodes <- colnames(atacSr)
  message(paste('--------finish creating seurat for disjoin', 
                lb, '-------'))
  return(atacSr)
}

create_atacSr_w_disjoin <- function(lb, metadata, 
                                     combined.peaks, annotation){
  atacAssay <- create_atacAssay(lb, metadata, annotation)
  atacSr <- create_atacDisjoin(atacAssay, lb, metadata, 
                                   combined.peaks, annotation)
  return(atacSr)
}

calculate_fragment <- function(atacSr, metadata){
  barcodes <- colnames(atacSr)
  lb <- unique(atacSr@meta.data$library)
  message(paste('---calculate total fragments for', lb, '----'))
  fragpath <- getFragPath(metadata, lb)
  totalFrag <- CountFragments(fragpath, cells = barcodes) 
  rownames(totalFrag) <- totalFrag$CB
  atacSr$fragments <- totalFrag[colnames(atacSr), "frequency_count"]
  atacSr$mononucleosomal <- totalFrag[colnames(atacSr), "mononucleosomal"]
  atacSr$nucleosome_free <- totalFrag[colnames(atacSr), "nucleosome_free"]
  atacSr$reads_count <- totalFrag[colnames(atacSr), "reads_count"]
  return(atacSr)
}
calculate_frip <- function(atacSr){
  message('----calculate FRiP ------')
  atacSr <- FRiP(
    object = atacSr,
    assay = 'peaks',
    total.fragments = 'fragments'
  )
  return(atacSr)
}
calculate_blacklist <- function(atacSr){
  message('----calculate blacklist ------')
  atacSr$blacklist_fraction <- FractionCountsInRegion(
    object = atacSr, 
    assay = 'peaks',
    regions = blacklist_hg38)
  return(atacSr)
}

calculate_nucleosome_signal <- function(atacSr){
  message('-----calculate nucleosomesignal -----')
  atacSr <- NucleosomeSignal(object = atacSr)   
  # #fragment ratio 147-294: <147 
  return(atacSr)
}
calculate_TSS_enrichment <- function(atacSr){
  message('----calculate TSS enrichment ------')
  atacSr <- TSSEnrichment(object = atacSr, fast = FALSE)
  return(atacSr)
}
calculate_metrics <- function(atacSr, metadata){
  atacSr <- calculate_fragment(atacSr, metadata)
  atacSr <- calculate_frip(atacSr)
  atacSr <- calculate_blacklist(atacSr)
  atacSr <- calculate_nucleosome_signal(atacSr)
  atacSr <- calculate_TSS_enrichment(atacSr)
  lb <- unique(atacSr$library)
  message (paste('---- finish calculate metric for', lb, '-----'))
  return(atacSr)
}

plot_isoutlier <- function(atacSce, yName, toColor, saveName, savePath){
  p <- plotColData(atacSce, y=yName,
                   colour_by=I(toColor)) + 
    theme(text = element_text(size = 30)) + 
    theme(axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 30))
  savePlot(p, filename = paste0(savePath, '/', saveName))
}
# atacSce is single cell experiment object
# atacSce = as.SingleCellExperiment(atacSr)
# nCount_peaks (low and high)
get_nCount_outlier <- function(atacSr, atacSce, figSavePath){
  ncountOutlier <- isOutlier(atacSr$nCount_peaks, type="both")
  ncountThred <- attr(ncountOutlier, "thresholds")
  lb <- unique(atacSr$library)
  plot_isoutlier(atacSce, yName='nCount_peaks',
                 toColor = ncountOutlier, 
                 saveName = paste0('scATAC_outlier_nucleosome_signal_', 
                                   lb, '.png'),
                 savePath = figSavePath)
  return(ncountThred)
}

# blacklist_fraction (high)
get_blacklist_outlier <- function(atacSr, atacSce, figSavePath){
blacklistOutlier <- isOutlier(atacSr$blacklist_fraction, type="higher")
blacklistThred <- attr(blacklistOutlier, "thresholds")
lb <- unique(atacSr$library)
plot_isoutlier(atacSce, yName="blacklist_fraction",
               toColor = blacklistOutlier, 
               saveName = paste0('scATAC_outlier_blacklist_fraction_', 
                                 lb, '.png'),
               savePath = figSavePath)
return(blacklistThred)
}
                                      

# nucleosome_signal (high)
get_nucleosome_outlier <- function(atacSr, atacSce, figSavePath){
  nucleOutlier <- isOutlier(atacSr$nucleosome_signal, type="higher")
  nucleThred <- attr(nucleOutlier, "thresholds")
  lb <- unique(atacSr$library)
  plot_isoutlier(atacSce, yName="nucleosome_signal",
                 toColor = nucleOutlier, 
                 saveName = paste0('scATAC_outlier_nucleosome_signal_', 
                                   lb, '.png'),
                 savePath = figSavePath)
  return(nucleThred)
}

get_TSS_outlier <- function(atacSr, atacSce, figSavePath){
  TSSOutlier <- isOutlier(atacSr$TSS.enrichment, type="lower")
  TSSThred <- attr(TSSOutlier, "thresholds")
  lb <- unique(atacSr$library)
  plot_isoutlier(atacSce, yName="TSS.enrichment",
                 toColor = TSSOutlier, 
                 saveName = paste0('scATAC_outlier_TSS_enrichment_', 
                                   lb, '.png'),
                 savePath = figSavePath)
  return(TSSThred)
}

sc_atac_filter_outliers <- function(atacSr, figSavePath){
  atacSce <- as.SingleCellExperiment(atacSr)
  ncountThred <- get_nCount_outlier(atacSr, atacSce, figSavePath)
  blacklistThred <- get_blacklist_outlier(atacSr, atacSce, figSavePath)
  nucleThred <- get_nucleosome_outlier(atacSr, atacSce, figSavePath)
  TSSThred <- get_TSS_outlier(atacSr, atacSce, figSavePath)
  atacSr <- subset(x = atacSr, 
                   subset = nCount_peaks < ncountThred[2] &             
                     nCount_peaks > ncountThred[1] &
                  blacklist_fraction < blacklistThred[2] & 
                  nucleosome_signal < nucleThred[2] &
                  TSS.enrichment > TSSThred[1])
  return(atacSr)
}   
  
sc_atac_normalize <- function(atacSr){
  atacSr <- RunTFIDF(atacSr) 
  atacSr <- FindTopFeatures(atacSr, min.cutoff = 'q0') 
  # keep feautures in n cells. q0 mean in top 100% cell
  atacSr <- RunSVD(atacSr)
  return(atacSr)
}

sc_atac_dim_redu <- function(atacSr){
  atacSr <- RunUMAP(object= atacSr, reduction = 'lsi', dims = 2:30) 
  # remove the 1st dim because it correlates w seq depth
  atacSr <- FindNeighbors(object= atacSr, reduction = 'lsi', dims = 2:30)
  atacSr <- FindClusters(object = atacSr, verbose = FALSE, algorithm = 3)
  return(atacSr)
}

get_gene_activity <- function(sr){
  gene.activities <- GeneActivity(sr)
  sr[['RNA']] <- CreateAssayObject(counts = gene.activities)
  sr <- NormalizeData(
  object = sr,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sr$nCount_RNA))
  DefaultAssay(sr) <- 'RNA'
  return(sr)
}


