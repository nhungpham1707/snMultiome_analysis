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
  # frags <- CreateFragmentObject(path = fragpath)
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
                 saveName = paste0(lb,'_scATAC_outlier_nucleosome_signal.png'),
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
               saveName = paste0(lb,'_scATAC_outlier_blacklist_fraction.png'),
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
                 saveName = paste0(lb,'_scATAC_outlier_nucleosome_signal.png'),
                 savePath = figSavePath)
  return(nucleThred)
}

get_TSS_outlier <- function(atacSr, atacSce, figSavePath){
  TSSOutlier <- isOutlier(atacSr$TSS.enrichment, type="lower")
  TSSThred <- attr(TSSOutlier, "thresholds")
  lb <- unique(atacSr$library)
  plot_isoutlier(atacSce, yName="TSS.enrichment",
                 toColor = TSSOutlier, 
                 saveName = paste0(lb,'_scATAC_outlier_TSS_enrichment.png'),
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

filter_outliers_healthyAtac <- function(atacSr, figSavePath){
  atacSce <- as.SingleCellExperiment(atacSr)
  ncountOutlier <- isOutlier(atacSr$nCount_peaks, type="both")
  ncountThred <- attr(ncountOutlier, "thresholds")
  lb <- unique(atacSr$tissue)
  plot_isoutlier(atacSce, yName='nCount_peaks',
                 toColor = ncountOutlier, 
                 saveName = paste0(lb,'_scATAC_outlier_nucleosome_signal.png'),
                 savePath = figSavePath)
  
  blacklistOutlier <- isOutlier(atacSr$blacklist_fraction, type="higher")
  blacklistThred <- attr(blacklistOutlier, "thresholds")
  plot_isoutlier(atacSce, yName="blacklist_fraction",
               toColor = blacklistOutlier, 
               saveName = paste0(lb,'_scATAC_outlier_blacklist_fraction.png'),
               savePath = figSavePath)
  # write report
  nOutlier_count <- length(which(ncountOutlier == TRUE))
  nOutlier_blk <- length(which(blacklistOutlier == TRUE))
  total_cells <- ncol(atacSr) 
  report_df <- data.frame(total_cells = total_cells,
              nOutlier_count = nOutlier_count,
              nOutlier_blk = nOutlier_blk)
  write.csv(report_df, paste0(report_dir, '/', lb, 'outlier_statistic.csv'))
  atacSr <- subset(x = atacSr, 
                   subset = nCount_peaks < ncountThred[2] &             
                     nCount_peaks > ncountThred[1] &
                  blacklist_fraction < blacklistThred[2])
  return(atacSr)
}   
  
sc_atac_normalize <- function(atacSr){
  atacSr <- RunTFIDF(atacSr) 
  atacSr <- FindTopFeatures(atacSr, min.cutoff = 'q0') 
  # keep feautures in n cells. q0 mean in top 100% cell
  atacSr <- RunSVD(atacSr)
  return(atacSr)
}

# The combined steps of TF-IDF followed by SVD are known as latent semantic indexing (LSI), and were first introduced for the analysis of scATAC-seq data by Cusanovich et al. 2015
sc_atac_dim_redu <- function(atacSr){
  atacSr <- RunUMAP(object= atacSr, reduction = 'lsi', dims = 2:30) 
  # remove the 1st dim because it correlates w seq depth
  atacSr <- FindNeighbors(object= atacSr, reduction = 'lsi', dims = 2:30)
  atacSr <- FindClusters(object = atacSr, verbose = FALSE, algorithm = 3)
  return(atacSr)
}

# cell cyle-, stress-, hemoglobin-, ribosome-, genes
# sometime affect the clustering. It is best to remove
# them from the variable gene lists before running dim reduc
# for  
remove_confounders <- function(sr, confounder_genes){
  srat_clean <- sr # copy so we can compare

VariableFeatures(srat_clean) <- setdiff( VariableFeatures(srat_clean), remove)
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


# get gene activity when fragment files are not available 
# ref https://biostars.org/p/9525109/
# annotation file is annotation output from getHG38annotation 
# sr is seurat object 
CreateGeneActivityMatrix <- function(
  sr,
  annotation.file,
  seq.levels = c(1:22, "X", "Y"),
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE
) {
    # SeuratObject::plan(multisession, workers = 8)
    # convert peak matrix to GRanges object
    peak.matrix <- as.matrix(GetAssayData(sr))
    peak.matrix <- GetAssayData(sr)
    peak.df <- rownames(x = peak.matrix)
    peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
    peak.df <- as.data.frame(x = peak.df)
    colnames(x = peak.df) <- c("chromosome", 'start', 'end')
    peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
    # if any peaks start at 0, change to 1
    # otherwise GenomicRanges::distanceToNearest will not work
    BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
    # get annotation file, select genes
    #   anno <- rtracklayer::import(con = annotation.file)
    # anno <- GenomeInfoDb::keepSeqlevels(x = anno, value = seq.levels, pruning.mode = 'coarse')

    # gtf <- rtracklayer::import(con = annotation.file)
    # gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')
    gtf <- annotation.file
    # change seqlevelsStyle if not the same
    # if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    #     GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
    # }
    # gtf.genes <- gtf[gtf$type == 'gene']
    # gtf.genes <- gtf[gtf$type == 'cds']
      gtf.genes <- gtf 
    # 
    # Extend definition up/downstream
    if (include.body) {
        gtf.body_prom <- Signac::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
    } else {
        gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
    }
    gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
    keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
    peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
    gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]

    # Some gtf rows will not have gene_name attribute
    # Replace it by gene_id attribute
    # 
    gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]

    peak.ids$gene.name <- gene.ids$gene_name
    peak.ids <- as.data.frame(x = peak.ids)
    peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
    annotations <- peak.ids[, c('peak', 'gene.name')]
    colnames(x = annotations) <- c('feature', 'new_feature')

    # collapse into expression matrix
    peak.matrix <- as(object = peak.matrix, Class = 'matrix')
    all.features <- unique(x = annotations$new_feature)

    if (future::nbrOfWorkers() > 1) {
        mysapply <- future.apply::future_sapply
    } else {
        mysapply <- ifelse(test = verbose, yes = pbapply::pbsapply, no = sapply)
    }
    newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){

    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature

    submat <- peak.matrix[features.use, ]

    if (length(x = features.use) > 1) {
        return(Matrix::colSums(x = submat))
    } else {
        return(submat)
    }
    })

    newmat <- t(x = newmat)

    rownames(x = newmat) <- all.features

    colnames(x = newmat) <- colnames(x = peak.matrix)

    return(as(object = newmat, Class = 'dgCMatrix'))
}

# annotation: output from getHg38annotation 
get_gene_activity_wo_fragFile <- function(sr, annotation){
  gene.activities <- CreateGeneActivityMatrix(sr, annotation)
  sr[['RNA']] <- CreateAssayObject(counts = gene.activities)
  sr <- NormalizeData(
  object = sr,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sr$nCount_RNA))
  DefaultAssay(sr) <- 'RNA'
  return(sr)
}
