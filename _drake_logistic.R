setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()


loadd(rna_w_tumor_label)
loadd(rna_group_sgr)
hthy = readRDS('output/healthy_data/merge_15_tissue.RDS')
hthy$all_bc <- colnames(hthy)

logistic_plan <- drake_plan(
    # test on descartes data ----
    hthy_rna = set_default_assay(hthy, assay = 'RNA'),
    hthy_nor = normalize_dim_plot_sr(hthy_rna, 
                        save_path = healthyDir,
                        lib_name = 'hthy_subset'),
    hthy_clus = clustering_rna_data(hthy_nor),
    hthy_75 = sampling_sr(hthy_clus, 75),
    to_keep = setdiff(colnames(hthy), colnames(hthy_75)),
    hthy_25 =  subset(hthy, subset = all_bc %in% to_keep),
    train_desc_rna_40000 = trainModel(GetAssayData(hthy_75), classes = hthy_75$cell_type, maxCells = 40000),

    predict_rna_desc = predictSimilarity(train_desc_rna_40000, 
        GetAssayData(rna_w_tumor_label), 
        classes = rna_w_tumor_label$cell_identity,
        minGeneMatch = 0.7),

    predict_desc_rna25hthy = predictSimilarity(train_desc_rna_40000, 
        GetAssayData(hthy_25),
        classes = hthy_25$cell_type, 
        minGeneMatch = 0.7),
    
    # test other reference data ----
    xu_atlas = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Terezinha_reference/xu_atlas_2023.RDS'),
    train_xu_atlas = trainModel(GetAssayData(xu_atlas), 
            classes = xu_atlas$final_annotation, 
            maxCells = 40000),
    predict_rnahm_xu_atlas = predictSimilarity(train_xu_atlas, 
        GetAssayData(rna_w_tumor_label), 
        classes = rna_w_tumor_label$cell_identity,
        minGeneMatch = 0.7),
    
    # test on rna without removing patient effect ---
    predict_rna_nohm_descartes = predictSimilarity(train_desc_rna_40000, 
        GetAssayData(rna_group_sgr), 
        classes = rna_group_sgr$RNA_snn_res.0.5,
        minGeneMatch = 0.7),

    predict_rna_nohm_xu = predictSimilarity(train_xu_atlas, 
        GetAssayData(rna_group_sgr),
        classes = rna_group_sgr$RNA_snn_res.0.5,
        minGeneMatch = 0.7)
)

make(logistic_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
