get_met <- function(gexSr){
    metaDf <- data.frame(sampleID = gexSr$sampleID,
                        gender = gexSr$curate_gender)
    rownames(metaDf) <- colnames(gexSr)
    return(metaDf)
}

# db - doublet or cells that fail from sample demultiplex 
removeDbAtac <- function(gexSr,atacSr){
    db <- setdiff(colnames(atacSr), colnames(gexSr))
    to_keep <- setdiff(colnames(atacSr), db)
    atacSr <- subset(x = atacSr, subset = barcodes %in% to_keep)
}

addMetAtac <- function(gexSr,atacSr){
    atacSr <- removeDbAtac(gexSr,atacSr)
    metaDf <- get_met(gexSr)
    atacSr <- AddMetaData(atacSr, metaDf)
    return(atacSr)
}

getMultiplexMet <- function(metadata, sLst, sr){
    allMet <- c()
      for (i in 1:length(sLst)){
            sInd <- which(metadata$sampleID == sLst[i])
            sMet <- metadata[sInd,]
            barcodes <- sr$barcodes[which(sr$sampleID == sLst[i])]
            sMet <- map_dfr(seq_len(length(barcodes)), ~sMet)
            sMet$barcodes <- barcodes
            allMet <- rbind(sMet,allMet)
        }
    rownames(allMet) <- allMet$barcodes
    allMet <- allMet[,-ncol(allMet)]
}

getSinglexMet <- function(metadata, sr){
    lib <- unique(sr$library)
    sMet <- metadata[which(metadata$name == lib),]
    barcodes <- sr$barcodes
    sMet <- map_dfr(seq_len(length(barcodes)), ~sMet)
    rownames(sMet) <- barcodes
    return(sMet)
}
addMetaFromFile <- function(metadata, sr){
    lib <- unique(sr$library)
    index <- which(metadata$name == lib)
    submeta <- metadata[index,]
    sLst <- unique(submeta$sampleID)
    if (length(sLst)>1) { # multiplex sample
        meta <- getMultiplexMet(metadata, sLst, sr)} else {
        meta <- getSinglexMet(metadata, sr)
        }
    sr <- AddMetaData(sr, meta)
}

# for samples that have similar genders in 
# multiplex libraries. metadata will assign base
# on souporcell # sr are gex after demultiplex ,after removing unknown and doublet souporcell
addMetaSoc <- function(metadata, sr){
    lib <- unique(sr$library)
    metadata <- metadata[which(metadata$name == lib),]
    gLst <- unique(sr$genotype)
    allMet <- c()
    for (i in 1:length(gLst)){
        barcodes <- sr$barcodes[which(sr$genotype == gLst[i])]
        sMet <- metadata[i,]
        sMet <- map_dfr(seq_len(length(barcodes)), ~sMet)
        sMet$barcodes <- barcodes
        allMet <- rbind(sMet,allMet)
        }
    rownames(allMet) <- allMet$barcodes
    allMet <- allMet[,-ncol(allMet)]
    sr <- AddMetaData(sr, allMet)
}


get_met_soc <- function(gexSr){
    metaDf <- data.frame(sampleID = gexSr$sampleID,
                        gender = gexSr$Gender)
    rownames(metaDf) <- colnames(gexSr)
    return(metaDf)
}
addMetSocAtac <- function(gexSr,atacSr){
    atacSr <- removeDbAtac(gexSr,atacSr)
    metaDf <- get_met_soc(gexSr)
    atacSr <- AddMetaData(atacSr, metaDf)
    return(atacSr)
}