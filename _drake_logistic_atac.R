setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

new_atachm_mx = readRDS('output/logistic_regression/atac_hm_features_above_300cells.RDS')
loadd(atac_hm_tumor_nona)
colnames(new_atachm_mx) <- colnames(atac_hm_tumor_nona)

logistic_atac <- drake_plan(
# atac ---
    dsc_atac = readRDS('output/healthy_data/dsc_atac_hthymrgDim.RDS'),
    dsc_ident = change_indent(dsc_atac, by ='cell_type'),
    dsc_markers = FindAllMarkers(object = dsc_ident, only.pos = T, logfc.threshold = 0.25),

    n = 7000,
    topfeatures = dsc_markers %>% 
      group_by(cluster) %>% 
      top_n(n = n, 
            wt = avg_log2FC),

    features_to_keep = topfeatures$gene,
    atac_features = rownames(new_atachm_mx),
    train_feature = intersect(features_to_keep, atac_features),
      
    sub_dsc_7k = subset(dsc_atac, features = train_feature),
    train_dsc_7k = trainModel(GetAssayData(sub_dsc_7k), class = sub_dsc_7k$cell_type, maxCell = 82300),
    sub_atac_7k = new_atachm_mx[rownames(new_atachm_mx) %in% train_feature,],
    p_dsc_7k = predictSimilarity(train_dsc_7k, sub_atac_7k, classes = atac_hm$cell_identity, 
                              logits = F, minGeneMatch = 0.7)


  )

  make(logistic_atac, lock_cache = FALSE,  lock_envir = FALSE)