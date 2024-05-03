message('---------- Start running _drake.R ------------')
setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
# Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

# read metadata ----
# filename <- '25012024_all_multiome_lib.csv'
filename <- '15042024_add_treatment_metadata.csv'
metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
# get all library ID
alLib <- unique(metadata$name)
alID <- alLib %>% map(splitName)
# get multiplex library that has the same gender
specialLib <- c("LX093_LX094_an_163")
specialLibInd <- grep(specialLib, metadata$name)
nospecialMet <- metadata[-specialLibInd,]
# get libraries that only demultiplex 
# with souporcell
soclId <- specialLib %>% map(splitName)
# get multiplex libraries list 
mulLib <- unique(nospecialMet$name[nchar(nospecialMet$souporcell_link) > 0])
# mulLib <- mulLib[1]
mulId <- mulLib %>% map(splitName)
# get single libraries list 
sngLib <- unique(nospecialMet$name[nchar(nospecialMet$souporcell_link) == 0])
# sngLib <- sngLib[1]
snglId <- sngLib %>% map(splitName)
# get all samples to make combine peaks
lbLst <- unique(c(specialLib, mulLib, sngLib)) 
idLst <- lbLst %>% map(splitName)


# define drake plan ----
# ----------------------------------------------------------
# plan to generate a common peak file to merge all atac samples with disjoin
# ----------------------------------------------------------

combine_peak_plan <- drake_plan(
  grFile = target(makeGrFile(metadata, lb),
                  transform = map(lb = !!lbLst, 
                                  id_var = !!idLst,
                                  .id = id_var)),
  allpeaks = target(disjoin(c(grFile)),
                    transform = combine(grFile,
                                        id.var = !!idLst,
                                        .id = id_var)),
  allpeaksFil = filterBadPeaks(allpeaks),
  allpeaksFilChr = filterChrBlacklist(allpeaksFil), 
)

# ----------------------------------------------------------
# plan to process library that pool together 2 samples of
# the same gender and the same tumor type
# steps include:
# - make sr with common peaks file
# - filter low quality cells w isoutlier
# - normalize, dim reduc
# - get gene activity to run single r and infercnv
# - run singler
# - prep for infercnv: generate gene matrix and annotation
#    file for infercnv 
# - sample demultiplex using souporcell
# - assign metadata
# ----------------------------------------------------------

process_special_lib_plan <- drake_plan(
  ## process ---
  atacSr_specialLib = create_atacSr_w_disjoin(specialLib, metadata, allpeaksFilChr, hg38),
  atacSrMe_specialLib = calculate_metrics(atacSr_specialLib, metadata),
  atacSrFil_specialLib = sc_atac_filter_outliers(atacSrMe_specialLib, figSavePath = atacProcessFigDir),
  atacSrNor_specialLib = sc_atac_normalize(atacSrFil_specialLib),
  atacSrDim_specialLib = sc_atac_dim_redu(atacSrNor_specialLib),

  atacSrGeA_specialLib = get_gene_activity(atacSrDim_specialLib),

  ## run singler ---
  sR_specialLib = run_singleR(atacSrGeA_specialLib),
  psR_specialLib = plot_singler(sR_specialLib, atacSrGeA_specialLib, save_path=atacCellSngRFigDir),
  atacGASgR_specialLib = get_sgR_label(sR_specialLib, atacSrGeA_specialLib),
  atacSgR_specialLib = get_sgR_label(sR_specialLib, atacSrDim_specialLib),

  ## prep for infercnv ---
  preInfer_specialLib = make_anno_count_mx(atacGASgR_specialLib, save_path=AtacInferInputDir),

  ## samples demultiplex with only souporcell ---
  h5Link_specialLib = get_h5_link(specialLib, metadata),
  gexSr_specialLib = create_GEX_seurat(h5Link_specialLib),
  gexSrstress_specialLib = remove_stress_genes(gexSr_specialLib, stress_gene_list),
  gexNor_specialLib = normalize_dim_sr(gexSrstress_specialLib),
  gexSoc_specialLib = assign_metadata_from_souporcell(gexNor_specialLib, metadata),
  vis_demul_specialLib = visualize_soc_demulti(gexSoc_specialLib, atacProcessFigDir),
  gexNoDb_specialLib = remove_soc_db_unknown(gexSoc_specialLib),
  gexSID_specialLib = addMetaSoc(metadata, gexNoDb_specialLib),
  atacDemul_specialLib = addMetSocAtac(gexSID_specialLib,atacSgR_specialLib),
  atacMeta_specialLib = addMetaFromFile(metadata, atacDemul_specialLib)

  # scroshi ----
  # atacspeciallib_gScroshi_demo = target(run_scROSHI_w_demo_data(sr = atacMeta_specialLib, 
  #       cols = my_cols, pt = 1,  save_name = 'w_demo_marker', save_path = atacScroshiDir)),
  
  # atacspeciallib_Scroshi_atrt = target(run_scROSHI_w_cancer_marker(sr = atacMeta_specialLib, 
  #   cols = my_cols, pt = 1, save_name = 'w_cancer_marker', save_path = atacScroshiDir))
  )

# ----------------------------------------------------------
# process all libraries. multiplex scATCseq libraries 
# have additional steps to demultiplex using 
# souporcell and gender from scRNAseq data
# ----------------------------------------------------------

  process_plan <- drake_plan(
  hg38 = getHg38Annotation(),
  # multiplex libraries -------------
  ## process ---
  atacSr = target(create_atacSr_w_disjoin(lb, metadata, allpeaksFilChr, hg38),
                  transform = map(lb = !!mulLib,
                                  id_var = !!mulId,
                                  .id = id_var)),
  atacSrMe = target(calculate_metrics(atacSr, metadata),
                    transform = map(atacSr,
                                    id_var = !!mulId,
                                    .id = id_var)),
  atacSrFil = target(sc_atac_filter_outliers(atacSrMe, 
                     figSavePath = atacProcessFigDir),
                     transform = map(atacSrMe,
                                     id_var = !!mulId,
                                     .id = id_var)),
  atacSrNor = target(sc_atac_normalize(atacSrFil),
                     transform = map(atacSrFil,
                                     id.var = !!mulId,
                                     .id = id_var)),
  atacSrDim = target(sc_atac_dim_redu(atacSrNor),
                     transform = map(atacSrNor,
                                     id.var = !!mulId,
                                     .id = id_var)),
  atacSrGeA = target(get_gene_activity(atacSrDim),
                     transform = map(atacSrDim,
                                     id.var = !!mulId,
                                     .id = id_var)),
  ## run singler ---
  sR = target(run_singleR(atacSrGeA),
              transform = map(atacSrGeA,
                              id_var = !!mulId,
                              .id = id_var)),
  psR = target(plot_singler(sR, atacSrGeA, save_path=atacCellSngRFigDir),
               transform = map(sR,atacSrGeA,
                               id_var = !!mulId,
                               .id = id_var)),
  atacSgR = target(get_sgR_label(sR, atacSrDim),
                   transform = map(sR,atacSrDim,
                                   id_var = !!mulId,
                                   .id = id_var)),
  atacGASgR = target(get_sgR_label(sR, atacSrGeA),
                   transform = map(sR,atacSrGeA,
                                   id_var = !!mulId,
                                  .id = id_var)),                          
  ## prep for infercnv ---
  preInfer = target(make_anno_count_mx(atacGASgR, save_path=AtacInferInputDir),
                    transform = map(atacGASgR,
                                    id_var = !!mulId,
                                    .id = id_var)),
  ## demultiplex -- 
  h5Link = target(get_h5_link(lb, metadata),
                  transform = map(lb = !!mulLib,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexSr = target(create_GEX_seurat(h5Link),
                 transform = map(h5Link,
                                 id.var = !!mulId,
                                 .id = id.var)),
  gexSrstress = target(remove_stress_genes(gexSr, stress_gene_list),
                       transform = map(gexSr,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexNor = target(normalize_dim_sr(gexSrstress),
                  transform = map(gexSrstress,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexSex = target(calculate_sex_genes_module(gexNor),
                  transform = map(gexNor,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexSexMeta = target(add_sex_metadata(gexSex),
                      transform = map(gexSex,
                                      id.var = !!mulId,
                                      .id = id.var)),
  gexSoc = target(assign_metadata_from_souporcell(gexSexMeta, metadata),
                  transform = map(gexSexMeta,
                                  id.var = !!mulId,
                                  .id = id.var)),
  vis_demul = target(visualize_soc_gender_demulti(gexSoc, 
                                                  atacProcessFigDir),
                     transform = map(gexSoc,id.var = !!mulId,
                                     .id = id.var)),
  gexNoDb = target(remove_doublet_n_unknown(gexSoc),
                   transform = map(gexSoc,
                                   id.var = !!mulId,
                                   .id = id.var)),
  gexNoUnknown = target(fix_unknown_gender_or_souporcell(gexNoDb),
                        transform = map(gexNoDb,
                                        id.var = !!mulId,
                                        .id = id.var)),
  gexNoContra = target(remove_contradict_cells(gexNoUnknown),
                       transform = map(gexNoUnknown,
                                       id.var = !!mulId,
                                       .id = id.var)),
  gexSID = target(assign_sample_name_and_tumor_type(metadata, gexNoContra),
                  transform = map(gexNoContra,
                                  id.var = !!mulId,
                                  .id = id.var)),
  # prepare demultiplex metadata 
  gexDemulMeta = target(generate_demultiplex_metadata(c(gexSID)),
                  transform = combine(gexSID,
                  id.var = !!mulId,
                  .id = id.var)),
  atacDemul = target(addMetAtac(gexSID,atacSgR),
                     transform = map(gexSID,atacSgR,
                                     id.var = !!mulId,
                                     .id = id.var)),
  atacMeta = target(addMetaFromFile(metadata, atacDemul),
                    transform = map(atacDemul,
                                    id.var = !!mulId,
                                    .id = id.var)),
  # scroshi ---
  # atacScroshi_demo = target(run_scROSHI_w_demo_data(sr = atacMeta, cols = my_cols, pt = 1, 
  #                     save_name = 'w_demo_marker', save_path = atacScroshiDir),
  #                     transform = map(atacMeta,
  #                       id.var = !!mulId,
  #                       .id = id.var)),
  
  # atacScroshi_atrt = target(run_scROSHI_w_cancer_marker(sr = atacMeta, cols = my_cols, pt = 1, 
  #                     save_name = 'w_cancer_marker', save_path = atacScroshiDir),
  #                     transform = map(atacMeta,
  #                       id.var = !!mulId,
  #                       .id = id.var)),
  # single libraries ----
  ## process ---
  atacSrsg = target(create_atacSr_w_disjoin(lb, metadata, allpeaksFilChr, hg38),
                    transform = map(lb = !!sngLib,
                                    id_var = !!snglId,
                                    .id = id_var)),
  atacSrMesg = target(calculate_metrics(atacSrsg, metadata),
                      transform = map(atacSrsg,
                                      id_var = !!snglId,
                                      .id = id_var)),
  atacSrFilsg = target(sc_atac_filter_outliers(atacSrMesg, 
                   figSavePath = atacProcessFigDir),
                       transform = map(atacSrMesg,
                                       id_var = !!snglId,
                                       .id = id_var)),
  atacSrNorsg = target(sc_atac_normalize(atacSrFilsg),
                       transform = map(atacSrFilsg,
                                       id.var = !!snglId,
                                       .id = id_var)),
  atacSrDimsg = target(sc_atac_dim_redu(atacSrNorsg),
                       transform = map(atacSrNorsg,
                                       id.var = !!snglId,
                                       .id = id_var)),
  atacSrGeAsg = target(get_gene_activity(atacSrDimsg),
                       transform = map(atacSrDimsg,
                                       id.var = !!snglId,
                                       .id = id_var)),
  ## run singler ---
  atac_sRsg = target(run_singleR(atacSrGeAsg),
                transform = map(atacSrGeAsg,
                                id_var = !!snglId,
                                .id = id_var)),
  p_atac_sRsg = target(plot_singler(atac_sRsg, atacSrGeAsg, save_path=atacCellSngRFigDir),
                 transform = map(atac_sRsg,atacSrGeAsg,
                                 id_var = !!snglId,
                                 .id = id_var)),
  atacsgSgR = target(get_sgR_label(atac_sRsg, atacSrDimsg),
                     transform = map(atac_sRsg,atacSrDimsg,
                                     id_var = !!snglId,
                                     .id = id_var)),
  atacsgGASgR = target(get_sgR_label(atac_sRsg, atacSrGeAsg),
                    transform = map(atac_sRsg,atacSrGeAsg,
                              id_var = !!snglId,
                              .id = id_var)),
  ## prep for infercnv ---
  preInfersg = target(make_anno_count_mx(atacsgGASgR, save_path=AtacInferInputDir),
                      transform = map(atacsgGASgR,
                                      id_var = !!snglId,
                                      .id = id_var)),
  ## add metadata ----
  atacMetasg = target(addMetaFromFile(metadata, atacsgSgR),
                      transform = map(atacsgSgR,
                                      id.var = !!snglId,
                                      .id = id.var)),
  # atacSgScroshi_demo = target(run_scROSHI_w_demo_data(sr = atacMetasg, cols = my_cols, pt = 1, 
  #                     save_name = 'w_demo_marker', save_path = atacScroshiDir),
  #                     transform = map(atacMetasg,
  #                       id.var = !!snglId,
  #                       .id = id.var)),
  
  # atacSgScroshi_atrt = target(run_scROSHI_w_cancer_marker(sr = atacMetasg, cols = my_cols, pt = 1, 
  #                     save_name = 'w_cancer_marker', save_path = atacScroshiDir),
  #                     transform = map(atacMetasg,
  #                       id.var = !!snglId,
  #                       .id = id.var)),
  ## merge atac-----
  mrgAtac = target(merge_pairwise(c(atacMeta_specialLib, atacMeta, atacMetasg),atcMrgDir),
            transform = combine(atacMeta,atacMetasg,
                    id.var = !!c(mulLib, sngLib),
                    .id = id.var)),
  mrgAtacNor = target(sc_atac_normalize(mrgAtac)),
  mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor)),
  mrgPtype = dimplot_w_nCell_label(mrgAtacDim, by = 'Subtype',atacMrgFigDir , col = my_cols2),
  mrgPsID = dimplot_w_nCell_label(mrgAtacDim, by = 'sampleID',atacMrgFigDir , col = my_cols2),

  atac_noNA = remove_na_cells(mrgAtacDim),

  # group singleR cell types ---
  atac_group_sgr = group_singleR_labels(atac_noNA),

  # prep mrg atac for infercnv
  mrgGA = get_gene_activity(atac_group_sgr),

  ## process rna ------------------------------------------
  h5Link_all = target(get_h5_link(lb, metadata),
                  transform = map(lb = !!alLib,
                                  id.var = !!alID,
                                  .id = id.var)),
  gexSr_all = target(create_GEX_seurat(h5Link_all),
                  transform = map(h5Link_all,
                              id.var = !!alID,
                              .id = id.var)),
  gexSrstress_all = target(remove_stress_genes(gexSr_all,
                          stress_gene_list),
                          transform = map(gexSr_all,
                              id.var = !!alID,
                              .id = id.var)),
  gexMito_all = target(check_mito_genes(gexSrstress_all),
                    transform = map(gexSrstress_all,
                    id.var = !!alID,
                    .id = id.var)),
  p_rnametric = target(rna_visualize_metric(gexMito_all,rnaFigDir),
                    transform = map(gexMito_all,
                    id.var = !!alID,
                    .id = id.var)),
  gexFilt = target(filter_rna_sr_w_isoutlier(gexMito_all, rnaFigDir),
                  transform = map(gexMito_all,
                  id.var = !!alID,
                  .id = id.var)),
  gexNor_all = target(normalize_dim_plot_sr(gexFilt, rnaFigDir, alLib),
                          transform = map(gexFilt,!!alLib,
                          id.var = !!alID,
                              .id = id.var)),
  gexClus = target(clustering_rna_data(gexNor_all),
                  transform = map(gexNor_all,
                  id.var = !!alID,
                              .id = id.var)),
  p_rnaClus = target(plot_cluster(gexClus, rnaFigDir, alLib), 
              transform = map(gexClus, !!alLib,
                  id.var = !!alID,
                              .id = id.var)),
  
   # run singleR rna ----
  rna_sgr = target(run_singleR(gexClus),
                transform = map(gexClus,
                                id_var = !!alID,
                                .id = id_var)),
  prna_sgr = target(plot_singler(rna_sgr, gexClus, save_path=rnaCellSngRFigDir),
                 transform = map(rna_sgr,gexClus,
                                 id_var = !!alID,
                                 .id = id_var)),
  gexClusSgr = target(get_sgR_label(rna_sgr, gexClus),
                transform = map(rna_sgr, gexClus,
                            id_var = !!alID,
                            .id = id_var)),
  # merge rna ---
   mrgRna_all = target(merge_pairwise(c(gexClusSgr), rnaMrgDir),
            transform = combine(gexClusSgr,
            id.var = !!alID,
            .id = id.var)),
  mrgRnaNor_all = normalize_dim_plot_sr(mrgRna_all, rnaMrgFigDir, lib_name = 'merge'),
  mrgRnaClu = clustering_rna_data(mrgRnaNor_all),
  rna_meta = assign_meta(metadata, gexDemulMeta,
  mrgRnaClu, save_name = paste0(rnaMrgDir, '/mrgRna_meta.RDS')),
  rna_noNA = remove_na_cells(rna_meta),
  rna_fix = fix_special_lib_rna(gexSID_specialLib, rna_noNA, gexNoDb_specialLib),
  rna_group_sgr = group_singleR_labels(rna_fix),
  # prep rna for infercnv -----
  # prepare count matrix 
  preInferRna = target(make_anno_count_mx(gexClusSgr, save_path = rnaInferInputDir ),
                    transform = map(gexClusSgr,
                  id.var = !!alID,
                  .id = id.var))
  # scroshi rna ---
  # rnaScroshi_demo = target(run_scROSHI_w_demo_data(sr = gexClusSgr, cols = my_cols, pt = 1, 
  #                     save_name = 'w_demo_marker', save_path = CellRnaScroshiDir),
  #                     transform = map(gexClusSgr,
  #                       id.var = !!alID,
  #                       .id = id.var)),
  
  # rnaScroshi_atrt = target(run_scROSHI_w_cancer_marker(sr = gexClusSgr, cols = my_cols, pt = 1, 
  #                     save_name = 'w_cancer_marker', save_path = CellRnaScroshiDir),
  #                     transform = map(gexClusSgr,
  #                       id.var = !!alID,
  #                       .id = id.var))

  
)

cluster_behavior_plan <- drake_plan( 
  # make atac sce ----
  atac.sce = make_sce(atac_group_sgr),
  # make rna sce -----
  rna.sce = make_sce(rna_group_sgr),
  # check cluster behavior ---- 
  ## Silhouette widths ------
  ### atac ----
  sil.atac = calculate_silhouette(atac.sce, reduce_method = 'LSI'),
  sil.atac.p =  ggplot(sil.atac, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") 
  + scale_colour_manual(values = my_cols) + 
    theme( panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size =20), 
          axis.title.y = element_text(size = 20)), 
  save_silAtac_p = savePlot('output/batchEffect/atac_cluster_behavior.png', sil.atac.p),
  silAtac_tab = table(Cluster=colLabels(atac.sce), sil.atac$closest),  # cluster 0, 1 and 4 have many cells that can easily mix with other clusters
  ### rna -----
  sil.rna = calculate_silhouette(rna.sce, reduce_method = 'PCA'),
  sil.rna.p = ggplot(sil.rna, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size =20), 
          axis.title.y = element_text(size = 20)),
  save_sil.rna.p = savePlot('output/batchEffect/rna_silhouette_cluster_behavior.png', sil.rna.p),
  table(Cluster=colLabels(rna.sce), sil.rna$closest),
  # cluster 7 & 23 have many cells that can easily mix with other clusters
  ## cluster purity ------
  ### atac ----
  pure.atac = calculate_purity(atac.sce, reduce_method = 'LSI'),
  pure.atac.p = ggplot(pure.atac , aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
  scale_colour_manual(values = my_cols) + 
    theme( panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size =20), 
          axis.title.y = element_text(size = 20)),
  save_pure.atac.p = savePlot('output/batchEffect/atac_cluster_purity.png', pure.atac.p),
  # To determine which clusters contaminate each other, we can identify the cluster with the most neighbors for each cell. In the table below, each row corresponds to one cluster; large off-diagonal counts indicate that its cells are easily confused with those from another cluster.
   ### rna ----
  pure.rna = calculate_purity(rna.sce, 'PCA'),
  pure.rna.p = ggplot(pure.rna, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") + 
    theme( panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          text = element_text(size =20), 
          axis.title.y = element_text(size = 20)),
  save_pure.rna.p = savePlot('output/batchEffect/rna_cluster_purity.png', pure.rna.p)
)

batch_detection_plan <- drake_plan(
  ## visualize batch ----
  ### atac ---
    atac_visBatchDate = plotBatchVis(atac.sce, batch = "Date.of.Library", save_path = batchAtacDir, col = my_cols),
    atac_visBatchLib = plotBatchVis(atac.sce, batch = 'library', save_path = batchAtacDir, col = my_cols),
    atac_visBatchSample = plotBatchVis(atac.sce, batch = 'sampleID', save_path = batchAtacDir, col = my_cols),
    atac_visBatchGender = plotBatchVis(atac.sce, batch = 'Gender', save_path = batchAtacDir, col = my_cols),
  ### rna -----
    rna_VisBatchDate = plotBatchVis(rna.sce, batch = "Date.of.Library", save_path = batchRnaDir, col = my_cols),
    rna_VisBatchLib = plotBatchVis(rna.sce, batch = 'library', save_path = batchRnaDir, col = my_cols),
    rna_VisBatchSample = plotBatchVis(rna.sce, batch = 'sampleID', save_path = batchRnaDir, col = my_cols),
    rna_VisBatchGender = plotBatchVis(rna.sce, batch = 'Gender', save_path = batchRnaDir, col = my_cols),

    ## calculate cms ---
    ### atac ----
    atac_cms = calculate_n_plot_cms(atac.sce, save_path = batchAtacDir, neighbors = 200, save_name = 'dj', 'LSI'),
    ### rna ----
    rna_cms = calculate_n_plot_cms(rna.sce, save_path = batchRnaDir, save_name = 'no_correction', neighbors = 200, 'PCA'), 

    ## check cell cycle ----
    rna_cellCycle = check_cell_cycle(rna_group_sgr, save_path = batchRnaDir)
)

# ------------------------------------------------------
# ------------------------------------------------------
batch_correction_plan <- drake_plan(
  # harmony ----------------------------------
  batch_factors = c('library', 'Individual.ID'),
  theta = seq(0, 0.5, by = 0.1),
  sigma = seq(0.1, 1, by = 0.2),
  ## atac -----------
  atac_lbsb = addLibSubcategory(atac_group_sgr), 

  # hm_atac = target(harmony_n_plot(atac_lbsb, batch_factor = batch_factors,theta = theta, sigma = sigma, save_path = batchAtacHarmonyDir, assay = 'peaks', reduction = 'lsi'), dynamic = cross(batch_factors, theta, sigma)),
  
  ## rna ----
  rna_lbsb = addLibSubcategory(rna_group_sgr),

  # hm_rna = target(harmony_n_plot(rna_group_sgr, batch_factor = batch_factors,theta = theta, sigma = sigma, save_path = batchRnaHarmonyDir), dynamic = cross(batch_factors, theta, sigma)),
  ## final hm ----
  final_theta = 0,
  final_sigma = 0.1, 
  # final_hm_rna = harmony_n_plot(rna_group_sgr, batch_factor = 'library', theta = final_theta,
  #  sigma = final_sigma, save_path = batchRnaHarmonyDir),
  
  final_hm_rna = RunHarmony(rna_group_sgr, group.by.vars = 'library', theta = 0),
  final_hm_rna_nb = FindNeighbors(object = final_hm_rna, reduction = "harmony", k.param = 30, dims = 1:30),
  final_hm_rna_clus = FindClusters(final_hm_rna_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  final_hm_rna_umap = RunUMAP(final_hm_rna_clus, dims = 1:30, reduction = 'harmony'),
  hm_rna_p = DimPlot(final_hm_rna_umap, group.by = 'Subtype', cols = my_cols),
  save_hm_rna_p = savePlot(paste0(batchRnaHarmonyDir, '/final_hm_subtype.png'), hm_rna_p),

  final_hm_atac = RunHarmony(atac_group_sgr, group.by.vars = 'library', theta = 0, reduction.use = 'lsi',  assay.use = 'peaks', project.dim = FALSE),
  final_hm_atac_nb = FindNeighbors(object = final_hm_atac, reduction = "harmony", k.param = 30),
  final_hm_atac_clus = FindClusters(final_hm_atac_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  final_hm_atac_umap = RunUMAP(final_hm_atac_clus, dims = 1:30, reduction = 'harmony'),
  hm_atac_p = DimPlot(final_hm_atac_umap, group.by = 'Subtype', cols = my_colors),
  save_hm_atac_p = savePlot(paste0(batchAtacHarmonyDir, '/final_hm_subtype.png'), hm_atac_p),
  # remove MHC genes and other confounding genes ----
  ## rna ---
  genes_to_remove = unique(c(genelists$chr6HLAgenes, genelists$hemo, genelists$stress, genelists$ribo)), 
  gene_to_retain = setdiff(rownames(rna_lbsb), genes_to_remove ),
  rna_noCF = subset(rna_lbsb, feature = gene_to_retain),
  rna_noCF_nor = normalize_dim_plot_sr(rna_noCF, rnaMrgFigDir, lib_name = 'merge_noCF'),
  rna_noCF_nor_clu = clustering_rna_data(rna_noCF_nor),
  rna_noCF_meta = assign_meta(metadata, gexDemulMeta,
  rna_noCF_nor_clu, save_name = paste0(rnaMrgDir, '/mrgRna_noCF_meta.RDS')),
  dim_rna_noCF_lib = DimPlot(rna_noCF_meta, group.by = 'library', cols = my_cols, raster = FALSE,pt.size = 1),
  save_dim_rna_noCF_lib = savePlot(paste0(rnaMrgFigDir, '/noCF_lib.png'), dim_rna_noCF_lib),
  dim_rna_noCF_sub = DimPlot(rna_noCF_meta, group.by = 'Subtype', cols = my_cols, raster = FALSE,pt.size = 1),
  save_dim_rna_noCF_sub = savePlot(paste0(rnaMrgFigDir, '/noCF_sub.png'), dim_rna_noCF_sub),
  hm_rna_noCF = harmony_n_plot(rna_noCF, batch_factor = 'library', theta = final_theta,
   sigma = final_sigma, save_path = batchRnaHarmonyDir)
 

  # calculate lisi ----
  ## before correction ----
  # # lisi_atac = calculate_lisi_from_sr(atac_lbsb, batch = 'library'),
  # lisi_rna = calculate_lisi_from_sr(rna_lbsb, batch = 'library'),
  # lisi_hm_rna = calculate_lisi_from_sr(final_hm_rna, batch = 'library')

)

cell_annotation_plan <- drake_plan(
  # scroshi merg rna ---
  hmRna_scroshi_demo = run_scROSHI_w_demo_data(sr = final_hm_rna_umap, cols = my_cols, pt = 1, save_name = 'hm_rna_w_demo_marker', save_path = CellRnaScroshiDir),
  
  hmRna_scroshi_atrt = run_scROSHI_w_cancer_marker(sr = final_hm_rna_umap, cols = my_cols, pt = 1, save_name = 'hm_rna_w_atrt', save_path = CellRnaScroshiDir),

  # scroshi merg atac ---
  hmAtac_scroshi_demo = run_scROSHI_w_demo_data(sr = final_hm_atac_umap, cols = my_cols, pt = 1, save_name = 'hm_atac_w_demo_marker', save_path = CellAtacScroshiDir),
  
  # hmAtac_scroshi_atrt = run_scROSHI_w_cancer_marker(sr = final_hm_atac, cols = my_cols, pt = 1, save_name = 'hm_atac_w_atrt', save_path = CellAtacScroshiDir),
  # infercnv res intepretation ----
  ## rna ---
  # did not run for LX069, LX099, LX183, LX189 
  # no clear difference w ref cells: LX 049, LX051, LX053, LX065, LX067, LX074, LX093, LX097, LX290
  # cut off after inspect manually using plot_aneuploidy_score function 
  atac_infercnv_cutoff = data.frame(cut_off = c(500, 300),
                                  lib = c('LX078_LX079_an_161', 'LX093_LX094_an_163')),
  # atac_infer = target(analyze_infercnv_res(c(atacMeta_specialLib, atacMeta, atacMetasg),infercnv_cut_off = atac_infercnv_cut_off, outlink = cellAtacInferDir),
  #           transform = combine(atacMeta,atacMetasg,
  #                   id.var = !!c(mulLib, sngLib),
  #                   .id = id.var)),
  

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
            "LX080_LX081_an_162")),
  # rna_infer_res_LX093 = target(analyze_infercnv_res(srat = c(gexClusSgr), rna_infercnv_cutoff,
  # infercnvlink = cellRnaIcnvdir, lib_to_check = "LX093_LX094_an_163"),
  # transform = combine(gexClusSgr)),
  # rna_infer = target(analyze_infercnv_res(srat = c(gexClusSgr),
  #                                           infercnv_cut_off = rna_infercnv_cutoff, 
  #                                           outlink = cellRnaIcnvdir ),
  #                                           transform = combine(gexClusSgr,
  #                                                       id.var = !!alID,
  #                                                       .id = id.var)),
                                                                                        

  # calculate markers ---
   rna_markers = FindAllMarkers(object = final_hm_rna, only.pos = T, logfc.threshold = 0.25),
   atac_markers = FindAllMarkers(object = final_hm_atac, only.pos = T, logfc.threshold = 0.25),

   ## show it visually:
# DoHeatmap(subset(srat, downsample = 50),
#           features = top10$gene,
#           group.colors=cluster.colors,
#           assay='RNA',
          # slot='scale.data'
          # ) + theme(axis.text.y=element_text(size=6)) # try here https://github.com/scgenomics/scgenomics-public.github.io/blob/main/docs/14-enrich/14-enrich.R
)
 
plan <- bind_plans(combine_peak_plan, process_special_lib_plan, process_plan, cell_annotation_plan, cluster_behavior_plan, batch_detection_plan, batch_correction_plan)

drake_config(plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
