# for some reason, drake combine gave error
# for now assign infercnv res manually and then add to drake again later 

# get cutoff -----
# the cutoff was retrieved manually from plot_aneuploid_score function 

source('test_plan.R')
  
rna_infercnv_cutoff = data.frame(cut_off = c(25,100,100,110,200,25,120,120,25,120),
    lib = c("LX093_LX094_an_163", 
            "LX290_LX291_an_423",  
            "LX065_LX066_an_155",
            "LX103_LX104_an_168", 
            "LX078_LX079_an_161", 
            "LX097_LX098_an_165", 
            "LX095_LX096_an_164", 
            "LX051_LX052_an_128", 
            "LX071_LX072_an_132",
            "LX080_LX081_an_162"))


# assign the cutoff to srat -----
## rna ----
message('check infer res for rna----')
rnaID <- rna_infercnv_cutoff$lib %>% map(splitName)
rna_infer_res <- c()
for (i in 1:length(rnaID)){
    message(paste0('check infer for rna', rnaID[i]))
    name <- paste0('gexClusSgr_', rnaID[i])
    loadd(name)
    sr <- get(name)
    sr_infer = analyze_infercnv_res(sr, rna_infercnv_cutoff, 
    infercnvlink = cellRnaIcnvdir, 
    lib_to_check = rna_infercnv_cutoff$lib[i])

    rna_infer_res <- data.frame(lib = sr_infer$library,
    barcode = sr_infer$barcodes,
    is_aneuploid = sr_infer$is_aneuploid,
    aneuploidy_score = sr_infer$aneuploidy_score) %>% rbind(rna_infer_res,.)
}

message('write rna_infer_res csv----')
write.csv(rna_infer_res, file = paste0(cellRnaIcnvdir, '/rna_infer_res.csv'), row.names = F)

## atac ----

atac_infercnv_cutoff = data.frame(cut_off = c(500, 300),
                                  lib = c('LX078_LX079_an_161', 'LX093_LX094_an_163'))

message('load atac srat')
loadd(atacMeta_LX078)
loadd(atacMeta_specialLib)

atac_infer_res <- c()
message('calculate infer res for lx078------')
lx078_infer = analyze_infercnv_res(atacMeta_LX078, atac_infercnv_cutoff, 
infercnvlink = cellAtacInferDir, 
lib_to_check = 'LX078_LX079_an_161')

atac_infer_res <- data.frame(lib = lx078_infer$library,
barcode = lx078_infer$barcodes,
is_aneuploid = lx078_infer$is_aneuploid,
aneuploidy_score = lx078_infer$aneuploidy_score) %>% rbind(atac_infer_res,.)

message('calculate infer res for lx093----')
lx093_infer = analyze_infercnv_res(atacMeta_specialLib, atac_infercnv_cutoff, 
infercnvlink = cellAtacInferDir, 
lib_to_check = 'LX093_LX094_an_163')

atac_infer_res <- data.frame(lib = lx093_infer$library,
barcode = lx093_infer$barcodes,
is_aneuploid = lx093_infer$is_aneuploid,
aneuploidy_score = lx093_infer$aneuploidy_score) %>% rbind(atac_infer_res,.)

message('write atac infer res to csv -------')
write.csv(atac_infer_res, file = paste0(cellAtacInferDir, '/atac_infer_res.csv'), row.names = F)

message('finish!')

# add to merge sr 
rna_infer_res <- read.csv('output/')

assign_infer_res_to_sr <- function(infer_res, sr){
    infer_res$lib_bc <- paste0(infer_res$lib, '_', infer_res$barcode)
    sr_bc <- data.frame(lib_bc = paste0(sr$library, '_', sr$barcodes))
    meta <- merge(infer_res, sr_bc, by= 'lib_bc')
    meta_to_add <- meta[,c('is_aneuploid', 'aneuploidy_score')]
    rownames(meta_to_add) <- meta$lib_bc
    sr <- AddMetaData(sr, meta_to_add)
}
rna_infer_res$lib_bc <- paste0(rna_infer_res$lib, '_', rna_infer_res$barcode)
mrg_bc <- paste0(mrg_rna$library, "_", mrg_rna$barcodes)
meta <- merge(rna_infer_res, m_bc, by.x = 'lib_bc')
rownames(meta) <- colnames(mrg_rna)
sr <- AddMetaData(mrg_rna, meta)