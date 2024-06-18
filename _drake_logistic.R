setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()


loadd(rna_w_tumor_label)
hthy = readRDS('output/healthy_data/merge_15_tissue.RDS')
hthy$all_bc <- colnames(hthy)

logistic_plan <- drake_plan(
    hthy_rna = set_default_assay(hthy, assay = 'RNA'),
    hthy_nor = normalize_dim_plot_sr(hthy_rna, 
                        save_path = healthyDir,
                        lib_name = 'hthy_subset'),
    hthy_clus = clustering_rna_data(hthy_nor),
    hthy_75 = sampling_sr(hthy_clus, 75),
    to_keep = setdiff(colnames(hthy), colnames(hthy_75)),
    hthy_25 =  subset(hthy, subset = all_bc %in% to_keep),
    train_rna_40000 = trainModel(GetAssayData(hthy_75), classes = hthy_75$cell_type, maxCells = 40000),

    predict_rna = predictSimilarity(train_rna_40000, GetAssayData(rna_w_tumor_label), 
    classes = rna_w_tumor_label$cell_identity, 
    minGeneMatch = 0.7),

    predict_25hthy = predictSimilarity(train_rna_40000, 
        GetAssayData(hthy_25),
        classes = hthy_25$cell_type, 
        minGeneMatch = 0.7)
)

make(logistic_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
