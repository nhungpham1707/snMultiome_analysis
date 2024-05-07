# calculate scores for cancer cells by adding results from singleR, inferCNV, scROSHI and cancer markers 
add_cancer_score_meta_cols <- function(sr){
    sr$sgr_cancer_score <- NA
    sr$infer_cancer_score <- NA
    sr$scroshi_cancer_score <- NA
    sr$marker_score <- NA
    return (sr)
}
assign_sgr_cancer_score <- function(sr){
    sr$sgr_cancer_score <- 1
    sr$sgr_cancer_score[sr$group_sgr_labels == 'immune_cells'] <- 0
    sr$sgr_cancer_score[sr$group_sgr_labels == 'stroma_cells'] <- 0
    return(sr)
}

assign_infercnv_cancer_score <- function(sr){
    sr$infer_cancer_score <- 0
    sr$infer_cancer_score[sr$is_aneuploid == 'yes'] <- 1
    return(sr)
}

assign_scroshi_cancer_score <- function(sr){
    sr$scroshi_cancer_score <- 0
    sr$scroshi_cancer_score[sr$celltype_final %in% c('ATRTTYR', 'ATRTSHH', 'ATRTMYC', 'RMS', 'Sysa')] <- 1
    return(sr)
}

assign_cancer_marker_score <- function(sr){

}

finalize_cancer_score <- function(sr){
    sr <- add_cancer_score_meta_cols(sr)
    sr <- assign_sgr_cancer_score(sr)
    sr <- assign_infercnv_cancer_score(sr)
    # sr <- assign_scroshi_cancer_score(sr)
    # sr <- assign_cancer_marker_score(sr)
    score_df <- data.frame(
        sgr = sr$sgr_cancer_score, 
        infer  = sr$infer_cancer_score, 
        scroshi = sr$scroshi_cancer_score, 
        marker = sr$marker_score)
    sr$final_cancer_score <- rowSums(score_df, na.rm = T)
    return (sr)
}

calculate_marker <- function(sr, marker, name){
    markers <- marker[marker %in% rownames(sr)]
    sr <- AddModuleScore(sr, features = markers, name = name)
    FeaturePlot(sr, features = paste0(name, '1'), label = T)
    
    return(sr)
}