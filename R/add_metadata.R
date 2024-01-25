get_met <- function(gexSr){
    metaDf <- data.frame(sampleID = gexSr$sampleID,
                        gender = gexSr$curate_gender)
    rownames(metaDf) <- colnames(gexSr)
    return(metaDf)
}

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

