get_met <- function(gexSr){
    metaDf <- data.frame(sampleID = gexSr$sampleID,
                        gender = gexSr$curate_gender)
    rownames(metaDf) <- colnames(gexSr)
    return(metaDf)
}

addMetAtac <- function(gexSr,AtacSr){
    metaDf <- get_met(gexSr)
    AtacSr <- AddMetadata(AtacSr, metaDf)
    return(AtacSr)
}

removeDbAtac <- function(AtacSr){
    
}