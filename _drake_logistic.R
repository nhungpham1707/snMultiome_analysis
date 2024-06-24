setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()


# loadd(rna_hm)
# loadd(rna)
hthy = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/all_celltypes.downsampled.filtered.RDS')
hthy$all_bc <- colnames(hthy)


logistic_plan <- drake_plan(
    rna = readRDS('output/sc_RNA/merge_all/rna.RDS'),
    rna_hm = readRDS('output/sc_RNA/merge_all/rna_hm.RDS'),
    # test on descartes data ----
    hthy_rna = set_default_assay(hthy, assay = 'RNA'),
    hthy_nor = normalize_dim_plot_sr(hthy_rna, 
                        save_path = healthyDir,
                        lib_name = 'hthy_subset'),
    hthy_clus = clustering_rna_data(hthy_nor),
    hthy_75 = sampling_sr(hthy_clus, 75, class_col = 'cell_type'),
    to_keep = setdiff(colnames(hthy_rna), colnames(hthy_75)),
    hthy_25 =  subset(hthy_rna, subset = all_bc %in% to_keep),
    train_desc_rna_40k = trainModel(GetAssayData(hthy_75), classes = hthy_75$cell_type, maxCells = 40000),

    p_rna_desc = predictSimilarity(train_desc_rna_40k, 
        GetAssayData(rna_hm), 
        classes = rna_hm$cell_identity,
        minGeneMatch = 0.7, logits = F),

    p_desc_rna25hthy = predictSimilarity(train_desc_rna_40k, 
        GetAssayData(hthy_25),
        classes = hthy_25$cell_type, 
        minGeneMatch = 0.7, logits = F),
    
     train_desc_rna_80k = trainModel(GetAssayData(hthy_75), classes = hthy_75$cell_type, maxCells = 80000),

    p_rna_desc_80k = predictSimilarity(train_desc_rna_80k, 
        GetAssayData(rna_hm), 
        classes = rna_hm$cell_identity,
        minGeneMatch = 0.7, logits = F),

    p_desc_rna25hthy_80k = predictSimilarity(train_desc_rna_80k, 
        GetAssayData(hthy_25),
        classes = hthy_25$cell_type, 
        minGeneMatch = 0.7, logits = F),
    # test other reference data ----
    ## xu at alas ----
    xu_atlas = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Terezinha_reference/xu_atlas_2023.RDS'),
    # train_xu_atlas = trainModel(GetAssayData(xu_atlas), 
    #         classes = xu_atlas$final_annotation, 
    #         maxCells = 40000),
    # predict_rnahm_xu_atlas = predictSimilarity(train_xu_atlas, 
    #     GetAssayData(rna_hm), 
    #     classes = rna_hm$cell_identity,
    #     minGeneMatch = 0.7),
    sub_xu = sampling_sr(xu_atlas, 100, type = 'number', class_col = 'final_annotation'),
    train_xu_sub = trainModel(GetAssayData(sub_xu), classes = sub_xu$final_annotation, minCells = 3),
    ## kameneva ---
    ka = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Terezinha_reference/adrenal.human.Kameneva_seurat.scrublet.rds'
),
    train_ka_10k = trainModel(
                    GetAssayData(ka),
                    classes = ka$fate,
                    maxCells = 10000),
    p_ka_10k = predictSimilarity(train_ka_10k, 
            GetAssayData(rna_hm),
            classes = rna_hm$cell_identity,
            minGeneMatch = 0.7, logits = F),
    train_ka_40k = trainModel(
            GetAssayData(ka),
            classes = ka$fate,
            maxCells = 40000),
    p_ka_40k = predictSimilarity(train_ka_40k,
        GetAssayData(rna_hm),
        classes = rna_hm$cell_identity,
        minGeneMatch = 0.7,
        logits = F),
    
    ## xi 2020 ---
    xi = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Jeff_rf/xi_2020.rds'),
    xi_nor = normalize_dim_plot_sr(xi, 
                        save_path = healthyDir,
                        lib_name = 'xi_2020'),
    xi_clus = clustering_rna_data(xi_nor),
    train_xi_40k = trainModel(
                    GetAssayData(xi_clus),
                    classes = xi_clus$cell_type,
                    maxCells = 40000),
    p_xi_40k = predictSimilarity(
                train_xi_40k,
                GetAssayData(rna_hm),
                classes = rna_hm$cell_identity,
                minGeneMatch = 0.7,
                logits = F),
    xi_nojuv = subset(xi, subset = stage_type != 'Juv_SC'),
    train_xi_nojuv = trainModel(
                    GetAssayData(xi_nojuv),
                    classes = xi_nojuv$cell_type,
                    maxCells = 10000),
    p_xi_nojuv = predictSimilarity(
                    train_xi_nojuv,
                    GetAssayData(rna_hm),
                    classes = rna_hm$cell_identity,
                    minGeneMatch = 0.7,
                    logits = F),
    # test on rna without removing patient effect ---
    # p_rna_nohm_descartes_40k = predictSimilarity(train_desc_rna_40k, 
    #     GetAssayData(rna), 
    #     classes = rna$RNA_snn_res.0.5,
    #     minGeneMatch = 0.7, logits = F),
    p_rna_no_hm_xi_40k = predictSimilarity(train_xi_40k,
            GetAssayData(rna),
            classes = rna$RNA_snn_res.0.5,
            minGeneMatch = 0.7,
            logits = F),
    p_rna_nohm_ka_10k = predictSimilarity(train_ka_10k, 
        GetAssayData(rna),
        classes = rna$RNA_snn_res.0.5,
        minGeneMatch = 0.7,
        logits = F),
    p_rna_nohm_ka_40k = predictSimilarity(train_ka_40k,
        GetAssayData(rna),
        classes = rna$cell_identity,
        minGeneMatch = 0.7,
        logits = F)

    # predict_rna_nohm_xu = predictSimilarity(train_xu_atlas, 
    #     GetAssayData(rna),
    #     classes = rna$RNA_snn_res.0.5,
    #     minGeneMatch = 0.7)
)

make(logistic_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
