setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()


# loadd(rna_hm)
# loadd(rna)

loadd(rna_w_tumor_label_newbc)
loadd(rna_hthymrg_clus)
dsc_rna <- rna_hthymrg_clus
dsc_rna$all_bc <- colnames(dsc_rna)
rna <- rna_w_tumor_label_newbc

loadd(atac_hthymrgDim)
 dsc_atac = atac_hthymrgDim
 dsc_atac$cell_bc = colnames(dsc_atac)

new_atachm_mx = readRDS('output/logistic_regression/atac_hm_features_above_300cells.RDS')
loadd(atac_hm_tumor_nona)
colnames(new_atachm_mx) <- colnames(atac_hm_tumor_nona)

logistic_plan <- drake_plan(
    
    # test on descartes data ----
    hthyrna_75 = sampling_sr(dsc_rna, 75, type = 'percent', class_col = 'cell_type'),

    to_keep_rna = setdiff(colnames(hthy_rna), colnames(hthyrna_75)),

    hthyrna_25 =  subset(hthy_rna, subset = all_bc %in% to_keep_rna),

    train_desc_75rna_allCells = trainModel(GetAssayData(hthyrna_75), classes = hthyrna_75$cell_type, maxCells = 82300),

    p_rna_75desc_allCells = predictSimilarity(train_desc_75rna_allCells, 
        GetAssayData(rna_hm), 
        classes = rna_hm$cell_identity,
        minGeneMatch = 0.7, logits = F),

    p_desc_75rnaOn25 = predictSimilarity(train_desc_75rna_allCells, 
        GetAssayData(hthyrna_25),
        classes = hthy_25$cell_type, 
        minGeneMatch = 0.7, logits = F),

    # atac ---
    dsc_atac_only_overlap = subset(dsc_atac, features = rownames(new_atachm_mx)),

    dscOnlyOverlap_markers = FindAllMarkers(object = dsc_atac_only_overlap, only.pos = T, logfc.threshold = 0.25),

    saveMarkers = saveRDS('output/logistic_regression/dscAtacOnlyOverlap_markers.RDS',dscOnlyOverlap_markers ),

    sub_dsc = extract_seurat_w_n_features(n = 3000, dscOnlyOverlap_markers, dsc_atac_only_overlap), 
     
    hthyatac_75 = sampling_sr(sub_dsc, 75, class_col = 'cell_type'),

    to_keep_atac = setdiff(colnames(sub_dsc), colnames(hthyatac_75)),

    hthyatac_25 =  subset(sub_dsc, subset = cell_bc %in% to_keep_atac),

    train_desc_75atac_allCells = trainModel(GetAssayData(hthyatac_75), classes = hthyatac_75$cell_type, maxCells = 82300),

    p_desc_75atacOn25 = predictSimilarity(train_desc_75atac_allCells, 
        GetAssayData(hthyatac_25),
        classes = hthyatac_25$cell_type, 
        minGeneMatch = 0.7, logits = F),

    sub_atac = new_atachm_mx[rownames(new_atachm_mx) %in% rownames(sub_dsc),],

    p_dsc_75atac = predictSimilarity(train_desc_75atac_allCells, 
   sub_atac, classes = atac_hm_tumor_nona$cell_identity, minGeneMatch = 0.7, logits = F)
    # test other reference data ----
    ## xu at alas ----
    # xu_atlas = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Terezinha_reference/xu_atlas_2023.RDS'),
    # train_xu_atlas = trainModel(GetAssayData(xu_atlas), 
    #         classes = xu_atlas$final_annotation, 
    #         maxCells = 40000),
    # predict_rnahm_xu_atlas = predictSimilarity(train_xu_atlas, 
    #     GetAssayData(rna_hm), 
    #     classes = rna_hm$cell_identity,
    #     minGeneMatch = 0.7),
    # sub_xu = sampling_sr(xu_atlas, 100, type = 'number', class_col = 'final_annotation'),
    # train_xu_sub = trainModel(GetAssayData(sub_xu), classes = sub_xu$final_annotation, minCells = 3),
    ## kameneva ---
#     ka = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Terezinha_reference/adrenal.human.Kameneva_seurat.scrublet.rds'
# ),
    # train_ka_10k = trainModel(
    #                 GetAssayData(ka),
    #                 classes = ka$fate,
    #                 maxCells = 10000),
    # p_ka_10k = predictSimilarity(train_ka_10k, 
    #         GetAssayData(rna_hm),
    #         classes = rna_hm$cell_identity,
    #         minGeneMatch = 0.7, logits = F),
    # train_ka_40k = trainModel(
    #         GetAssayData(ka),
    #         classes = ka$fate,
    #         maxCells = 40000),
    # p_ka_40k = predictSimilarity(train_ka_40k,
    #     GetAssayData(rna_hm),
    #     classes = rna_hm$cell_identity,
    #     minGeneMatch = 0.7,
    #     logits = F),
    
    # ## xi 2020 ---
    # xi = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Jeff_rf/xi_2020.rds'),
    # xi_nor = normalize_dim_plot_sr(xi, 
    #                     save_path = healthyDir,
    #                     lib_name = 'xi_2020'),
    # xi_clus = clustering_rna_data(xi_nor),
    # train_xi_40k = trainModel(
    #                 GetAssayData(xi_clus),
    #                 classes = xi_clus$cell_type,
    #                 maxCells = 40000),
    # p_xi_40k = predictSimilarity(
    #             train_xi_40k,
    #             GetAssayData(rna_hm),
    #             classes = rna_hm$cell_identity,
    #             minGeneMatch = 0.7,
    #             logits = F),
    # xi_nojuv = subset(xi, subset = stage_type != 'Juv_SC'),
    # train_xi_nojuv = trainModel(
    #                 GetAssayData(xi_nojuv),
    #                 classes = xi_nojuv$cell_type,
    #                 maxCells = 10000),
    # p_xi_nojuv = predictSimilarity(
    #                 train_xi_nojuv,
    #                 GetAssayData(rna_hm),
    #                 classes = rna_hm$cell_identity,
    #                 minGeneMatch = 0.7,
    #                 logits = F),
    # test on rna without removing patient effect ---
    # p_rna_nohm_descartes_40k = predictSimilarity(train_desc_rna_40k, 
    #     GetAssayData(rna), 
    #     classes = rna$RNA_snn_res.0.5,
    #     minGeneMatch = 0.7, logits = F),
    # p_rna_no_hm_xi_40k = predictSimilarity(train_xi_40k,
    #         GetAssayData(rna),
    #         classes = rna$RNA_snn_res.0.5,
    #         minGeneMatch = 0.7,
    #         logits = F),
    # p_rna_nohm_ka_10k = predictSimilarity(train_ka_10k, 
    #     GetAssayData(rna),
    #     classes = rna$RNA_snn_res.0.5,
    #     minGeneMatch = 0.7,
    #     logits = F),
    # p_rna_nohm_ka_40k = predictSimilarity(train_ka_40k,
    #     GetAssayData(rna),
    #     classes = rna$cell_identity,
    #     minGeneMatch = 0.7,
    #     logits = F),

    # # predict_rna_nohm_xu = predictSimilarity(train_xu_atlas, 
    # #     GetAssayData(rna),
    # #     classes = rna$RNA_snn_res.0.5,
    # #     minGeneMatch = 0.7)

    # # atac ---
    # mrg_dsc_atac = readRDS('output/logistic_regression/mrg_descartes_atac.RDS'),
    # dsc_atac_data = subset(mrg_dsc_atac, subset = source == 'descartes'),
    # atac_data = subset(mrg_dsc_atac, subset = source == 'nhung_etal'),
    # train_dsc_atac = trainModel(GetAssayData(dsc_atac_data), classes = dsc_atac_data$cell_type),

    # atacdata_newbc = paste0(atac_data$barcodes, '_', atac_data$library),
    # atac_data_newbc = RenameCells(atac_data, new.names = atacdata_newbc),
    # atac_data_w_tumor_label = assign_cross_labels(des_sr = atac_data_newbc, source_sr = rna_w_tumor_label_newbc, 
    # label_col = 'cell_identity'),
    # p_dsc_atac = predictSimilarity(train_dsc_atac, 
    #     GetAssayData(atac_data_w_tumor_label),
    #     classes = atac_data_w_tumor_label$cell_identity,
    #     minGeneMatch = 0.7, 
    #     logits = F)

)



make(logistic_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
