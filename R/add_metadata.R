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

# add metadata after sample demultiplex to mrg object 

prepare_sr_meta <- function(sr){
    colnames(sr@meta.data)
    sr_barcodes <- data.frame(libary = sr$library,
                      m_barcodes = colnames(sr),
                      barcodes = sr$barcodes)
    sr_barcodes$lib_bar <- paste0(sr_barcodes$libary, '_', sr_barcodes$barcodes)
    return(sr_barcodes)
}

# not all libraries were multiplex, 
# sample Ids and other metadata can be directly
# assigned to these libraries  
get_singlex_meta <- function(metadata, demul_meta, sr_barcodes){
    splex <- setdiff(metadata$name, demul_meta$library)
    splex_meta <- metadata[metadata$name %in% splex,]
    splex_sr <- sr_barcodes[sr_barcodes$libary %in% splex_meta$name,]
    splex_sr_meta <- merge(splex_sr, splex_meta, by.x = 'libary', by.y = 'name')
    splex_sr_meta <- within(splex_sr_meta, rm('libary', 'barcodes', 'lib_bar'))
    return(splex_sr_meta) # somehow ATAC_lib not there 

}

get_multiplex_meta <- function(metadata, demul_meta, sr_barcodes){
    mplex_meta <- metadata[metadata$sampleID %in% demul_meta$sampleID,]
    # mplex_meta <- mplex_meta[,2:ncol(mplex_meta)]
    # head(mplex_meta)
    mplex_meta <- merge(mplex_meta, demul_meta, by.x = 'sampleID', by.y = 'sampleID')
    mplex_meta$lib_bar <- paste0(mplex_meta$library, '_', mplex_meta$barcodes)
    mplex_sr_meta <- merge(mplex_meta, sr_barcodes, by= 'lib_bar')
    mplex_sr_meta <- within(mplex_sr_meta, rm('lib_bar', 'barcodes.x', 'barcodes.y','name', 'library', 'libary'))
    return(mplex_sr_meta)
}

combine_meta <- function(mplex_sr_meta, splex_sr_meta){
    print ('colnames of mplex_sr_meta are')
    print(colnames(mplex_sr_meta))
    print(ncol(mplex_sr_meta))
    print('colnames of splex_sr_meta are')
    print(colnames(splex_sr_meta))
    print(ncol(splex_sr_meta))
    all_meta <- rbind(mplex_sr_meta, splex_sr_meta)
    all_meta2 <- all_meta[!duplicated(all_meta$m_barcodes),] #somehow LX094 barcodes are duplicated
    rownames(all_meta2) <- all_meta2$m_barcodes
    meta_to_add <- within(all_meta2, rm('m_barcodes'))
    return(meta_to_add)
}

assign_meta <- function(metadata, demul_meta, sr, save_name){
    sr_barcodes <- prepare_sr_meta(sr)
    splex_sr_meta <- get_singlex_meta(metadata, demul_meta, sr_barcodes)
    mplex_sr_meta <- get_multiplex_meta(metadata, demul_meta, sr_barcodes)
    meta_to_add <- combine_meta(splex_sr_meta, mplex_sr_meta)
    sr_meta <- AddMetaData(sr, meta_to_add)
    saveRDS(sr_meta, save_name)
    return(sr_meta)
}

add_treatment_meta <- function(sr, meta){
    columns <- setdiff(colnames(meta), colnames(sr@meta.data))
    submeta <- meta[, c('sampleID',columns)]
    bc_sample <- data.frame(bc = colnames(sr),
                            sampleID = sr$sampleID)
    meta_all <- merge(bc_sample, submeta, by = 'sampleID', all = 'T')
    meta_to_add <- meta_all[, columns]
    rownames(meta_to_add) <- meta_all$bc
    sr_meta <- AddMetaData(sr, meta_to_add)
    return(sr_meta)
}
 
