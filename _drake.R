setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

## read metadata ----
filename <- '25012024_all_multiome_lib.csv'
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
# region ---
lbLst <- unique(c(specialLib, mulLib, sngLib)) 
idLst <- lbLst %>% map(splitName)
# endregion
## define plan ----

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


process_special_lib_plan <- drake_plan(
  ## process ---
  atacSr_specialLib = create_atacSr_w_disjoin(specialLib, metadata, allpeaksFilChr, hg38),
  atacSrMe_specialLib = calculate_metrics(atacSr_specialLib, metadata),
  atacSrFil_specialLib = sc_atac_filter_outliers(atacSrMe_specialLib, 
                               figSavePath = atacProcessFigDir),
  atacSrNor_specialLib = sc_atac_normalize(atacSrFil_specialLib),
  atacSrDim_specialLib = sc_atac_dim_redu(atacSrNor_specialLib),
  atacSrGeA_specialLib = get_gene_activity(atacSrDim_specialLib),
  ## run singler ---
  sR_specialLib = run_singleR(atacSrGeA_specialLib),
  psR_specialLib = plot_singler(sR_specialLib, atacSrGeA_specialLib, 
                                save_path=atacCellSngRFigDir),
  atacSgR_specialLib = get_sgR_label(sR_specialLib, atacSrDim_specialLib),
  ## prep for infercnv ---
  preInfer_specialLib = make_anno_count_mx(sR_specialLib,atacSrGeA_specialLib, 
                                           save_path=AtacInferInputDir),
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

)

process_plan <- drake_plan(
  hg38 = getHg38Annotation(),
  # multiplex sample
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
  ## prep for infercnv ---
  preInfer = target(make_anno_count_mx(sR,atacSrGeA, save_path=AtacInferInputDir),
                    transform = map(sR,atacSrGeA,
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
  saveGexDemulMeta = saveRDS(gexDemulMeta, paste0(rnaProcessDir, '/gex_demultiplex_metadata.RDS')),
  atacDemul = target(addMetAtac(gexSID,atacSgR),
                     transform = map(gexSID,atacSgR,
                                     id.var = !!mulId,
                                     .id = id.var)),
  atacMeta = target(addMetaFromFile(metadata, atacDemul),
                    transform = map(atacDemul,
                                    id.var = !!mulId,
                                    .id = id.var)),
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
  sRsg = target(run_singleR(atacSrGeAsg),
                transform = map(atacSrGeAsg,
                                id_var = !!snglId,
                                .id = id_var)),
  psRsg = target(plot_singler(sRsg, atacSrGeAsg, save_path=atacCellSngRFigDir),
                 transform = map(sRsg,atacSrGeAsg,
                                 id_var = !!snglId,
                                 .id = id_var)),
  atacsgSgR = target(get_sgR_label(sRsg, atacSrDimsg),
                     transform = map(sRsg,atacSrDimsg,
                                     id_var = !!snglId,
                                     .id = id_var)),
  ## prep for infercnv ---
  preInfersg = target(make_anno_count_mx(sRsg,atacSrGeAsg, save_path=AtacInferInputDir),
                      transform = map(sRsg,atacSrGeAsg,
                                      id_var = !!snglId,
                                      .id = id_var)),
  ## add metadata ----
  atacMetasg = target(addMetaFromFile(metadata, atacsgSgR),
                      transform = map(atacsgSgR,
                                      id.var = !!snglId,
                                      .id = id.var)),
  ## merge atac-----
  mrgAtac = target(merge_pairwise(c(atacMeta_specialLib, atacMeta, atacMetasg),atcMrgDir),
            transform = combine(atacMeta,atacMetasg,
                    id.var = !!c(mulLib, sngLib),
                    .id = id.var)),
  mrgAtacNor = target(sc_atac_normalize(mrgAtac)),
  mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor)),
  mrgPtype = dimplot_w_nCell_label(mrgAtacDim, by = 'Subtype',atacMrgFigDir , col = my_cols2),
  mrgPsID = dimplot_w_nCell_label(mrgAtacDim, by = 'sampleID',atacMrgFigDir , col = my_cols2),
  # # prep mrg atac for infercnv
  mrgGA = get_gene_activity(mrgAtacDim),
  ## process rna ----
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

  # # merge rna with normalized samples -----
  mrgRna_all = target(merge_pairwise(c(gexClus), rnaMrgDir),
            transform = combine(gexClus,
            id.var = !!alID,
            .id = id.var)),
  mrgRnaNor_all = normalize_dim_plot_sr(mrgRna_all, rnaMrgFigDir, lib_name = 'merge'),
  saverna_all = saveRDS(mrgRnaNor_all, paste0(rnaMrgDir, '/mrgRna.RDS')),
  mrgRnaClu = clustering_rna_data(mrgRnaNor_all),
  p_mrgRna_all = plot_cluster(mrgRnaClu, rnaMrgFigDir, save_name = 'merge' ),
  saveplotmrgrna = savePlot(paste0(rnaMrgFigDir, '/mrgRNA.png'), p_mrgRna_all),

   # merge rna without normalized samples ----
  mrgRna_noNOr_all = target(merge_pairwise(c(gexFilt), rnaMrgDir),
            transform = combine(gexFilt,
            id.var = !!alID,
            .id = id.var)),
  mrgRnaNor_noNOr_all = normalize_dim_plot_sr(mrgRna_noNOr_all, rnaMrgFigDir, lib_name = 'merge'),
  saverna_noNOr_all = saveRDS(mrgRnaNor_noNOr_all, paste0(rnaMrgDir, '/mrgRna_noNOr.RDS')),
  mrgRnaClu_noNOr = clustering_rna_data(mrgRnaNor_noNOr_all),
  saverna_noNOr_clu = saveRDS(mrgRnaClu_noNOr, paste0(rnaMrgDir, '/mrgRna_noNOrClu.RDS')),
  rna_meta = assign_meta(metadata, gexDemulMeta,
  mrgRnaClu_noNOr, save_name = paste0(rnaMrgDir, '/mrgRna_meta.RDS')),
   # run singleR rna ---
  rnaMrgNoNor_sgr = run_singleR(mrgRnaClu_noNOr),
  savernaSgR = saveRDS(rnaMrgNoNor_sgr, paste0(cellRNAsingRdir, '/mrgRna_noNOr_singleR.RDS')),
  rnaMrgSgr = get_sgR_label(rnaMrgNoNor_sgr, rna_meta),
  # prep rna for infercnv -----
  # prepare count matrix 
  preInferRna = target(make_anno_count_rna_mx(rnaMrgSgr, gexClus, save_path = rnaInferInputDir ),
                    transform = map(gexClus,
                  id.var = !!alID,
                  .id = id.var))
)

cell_annotation_plan <- drake_plan(
  # scroshi
   mrgRna_scroshi_demo = run_scROSHI_w_demo_data(sr = rnaMrgSgr, cols = my_cols, pt = 1, save_name = 'merge_all_rna_w_demo_marker', save_path = CellRnaScroshiDir),
  
  mrgRna_scroshi_atrt = run_scROSHI_w_atrt_data(sr = mrgRna_scroshi_demo, cols = my_cols, pt = 1, save_name = 'merge_rna_w_atrt', save_path = CellRnaScroshiDir)
)

cluster_behavior_plan <- drake_plan(
  # remove NA cells in rna ----
  rna_noNA = remove_na_cells(rna_meta),
  # make atac sce ----
  atac.sce = make_sce(mrgAtacDim),
  # make rna sce -----
  rna.sce = make_sce(rna_noNA),
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
)

batch_correction_plan <- drake_plan(
  # try with harmony -----
  ## atac ----
  hm_lib.atac = RunHarmony(mrgAtacDim, group.by.vars = 'library', reduction.use = 'lsi',  assay.use = 'peaks', project.dim = FALSE),
  hm_lib.atac_umap = RunUMAP(hm_lib.atac, dims = 2:30, reduction = 'harmony'),
  ## rna ----
  hm_lib.rna = RunHarmony(rnaMrgSgr, group.by.vars = 'library', reduction.use = 'pca', project.dim = FALSE),
  hm_lib.rna_nb = FindNeighbors(object = hm_lib.rna, reduction = "harmony"),
  hm_lib.rna_clus = FindClusters(hm_lib.rna_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  hm_lib.rna_umap = RunUMAP(hm_lib.rna_clus, dims = 2:30, reduction = 'harmony'),
  hm_rna_singr_p = DimPlot(hm_lib.rna_umap, group.by = 'singleR_labels', raster = FALSE, cols = my_cols),
  svae_hm_rna_singr = savePlot('output/batchEffect/hm_rna_singr.png', hm_rna_singr_p),
  hm_rna_lib_p = DimPlot(hm_lib.rna_umap, group.by = 'library', raster = FALSE, pt.size = 0.1, cols = my_cols),
  save_hmrna_lib = savePlot('output/batchEffect/hm_rna_lib.png', hm_rna_lib_p),
  hm_rna_type_p = DimPlot(hm_lib.rna_umap, group.by = 'Subtype', raster = FALSE, pt.size = 0.1, cols = my_cols),
  save_hm_rna_type = savePlot('output/batchEffect/hm_rna_type.png', hm_rna_type_p),
  
  # harmony on new category (lib+sub) ---
  ## atac ----
  mrgAtacLbSb = addLibSubcategory(mrgAtacDim),
  saveh5_atac = save_h5ad(mrgAtacLbSb, save_path = atcMrgDir, save_name = 'atac_lbsb'),
  hm_lbsb.atac = RunHarmony(mrgAtacLbSb, group.by.vars = 'lbsb', reduction.use = 'lsi',  assay.use = 'peaks', project.dim = FALSE),
  hm_lbsb.atac_nb = FindNeighbors(object = hm_lbsb.atac, reduction = "harmony", k.param = 30),
  hm_lbsb.atac_clus = FindClusters(hm_lbsb.atac_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  hm_lbsb.atac_umap = RunUMAP(hm_lbsb.atac_clus, dims = 2:30, reduction = 'harmony'),
  hm_lbsb_atc_singr_p = DimPlot(hm_lbsb.atac_umap, group.by = 'singleR_labels', raster = FALSE, cols = my_cols),
  save_hm_lbsb_atac_singr = savePlot(paste0(batchAtacHarmonyDir,'/hm_lbsb_atac_singr.png'), hm_lbsb_atc_singr_p),
  hm_atac_lbsb_p = DimPlot(hm_lbsb.atac_umap, group.by = 'library', raster = FALSE, pt.size = 0.1, cols = my_cols),
  save_hmrna_lib_p = savePlot(paste0(batchAtacHarmonyDir,'/hm_rna_lib.png'), hm_atac_lbsb_p),
  hm_atac_lbsb_type_p = DimPlot(hm_lbsb.atac_umap, group.by = 'Subtype', raster = FALSE, pt.size = 0.1, cols = my_cols),
  save_hm_atac_lbsb_type_p = savePlot(paste0(batchAtacHarmonyDir,'/hm_rna_type.png'), hm_atac_lbsb_type_p),
  
  ## rna ----
  mrgRnaLbSb = addLibSubcategory(rnaMrgSgr),
  saveh5_rna = save_h5ad(mrgRnaLbSb, save_path = rnaMrgDir, save_name = 'rna_lbsb'),
  hm_lbsb.rna = RunHarmony(mrgRnaLbSb, group.by.vars = 'lbsb', reduction.use = 'pca', project.dim = FALSE),
  hm_lbsb.rna_nb = FindNeighbors(object = hm_lbsb.rna, reduction = "harmony", k.param = 30),
  hm_lbsb.rna_clus = FindClusters(hm_lbsb.rna_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  hm_lbsb.rna_umap = RunUMAP(hm_lbsb.rna_clus, dims = 1:30, reduction = 'harmony'),
  hm_lbsb_rna_singr_p = DimPlot(hm_lbsb.rna_umap, group.by = 'singleR_labels', raster = FALSE, cols = my_cols),
  save_hm_lbsb_rna_singr = savePlot(paste0(batchRnaHarmonyDir,'/hm_lbsb_rna_singr.png'), hm_lbsb_rna_singr_p),
  hm_rna_lbsb_p = DimPlot(hm_lbsb.rna_umap, group.by = 'library', raster = FALSE, pt.size = 0.1, cols = my_cols),
  save_hmrna_lib_p = savePlot(paste0(batchRnaHarmonyDir,'/hm_rna_lib.png'), hm_rna_lbsb_p),
  hm_rna_lbsb_type_p = DimPlot(hm_lbsb.rna_umap, group.by = 'Subtype', raster = FALSE, pt.size = 0.1, cols = my_cols),
  save_hm_rna_lbsb_type_p = savePlot(paste0(batchRnaHarmonyDir,'/hm_rna_type.png'), hm_rna_lbsb_type_p)

)
# convert seurat to anndata https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html 


plan <- bind_plans(combine_peak_plan,process_special_lib_plan,  process_plan, cell_annotation_plan, cluster_behavior_plan, batch_detection_plan, batch_correction_plan)


# options(clustermq.scheduler = "multicore") # nolint
# make(plan, parallelism = "clustermq", jobs = 1, lock_cache = FALSE)
make(plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
# vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_cleancode_pipeline.png', font_size = 20 )

# plot cms score

# loadd(atac_cms)
# loadd(rna_cms)
# his_p <- visHist(atac_cms)
# metric_p <- visMetric(atac_cms, metric_var = 'cms_smooth.date')
# plot_grip(his_p, metric_p, ncol = 2)

# names(colData(atac_cms))
# cms_df <- data.frame(date = atac_cms$cms.dj_date_k200,
#         smooth_date = atac_cms$cms_smooth.dj_date_k200,
#         lib = atac_cms$cms.dj_lib_k200,
#         smooth_lib = atac_cms$cms_smooth.dj_lib_k200,
#         patient = atac_cms$cms.dj_patient_k200,
#         smooth_patient = atac_cms$cms_smooth.dj_patient_k200 )

# saveRDS(cms_df, paste0(batchAtacDir,'/cms_df.rds'))
# his_p <- visHist(atac_cms, metric = 'cms.dj_date_k200' )
# metric_p <- visMetric(atac_cms, metric_var = 'cms_smooth.dj_date_k200')
# # p <- plot_grip(his_p, metric_p, ncol = 2)
# savePlot(paste0(batchAtacDir, '/dj_date_hist.png'), his_p )
# savePlot(paste0(batchAtacDir, '/dj_date_metric.png'), metric_p )
# print('rna colnames')
# names(colData(rna_cms))

# his_rna <- visHist(rna_cms, metric = "cms.no_correction_patient_k200" )
# savePlot(paste0(batchAtacDir, '/rna_patient_hist.png'), his_rna)
# metric_rna <- visMetric(rna_cms, metric_var = "cms_smooth.no_correction_patient_k200")
# savePlot(paste0(batchAtacDir, '/rna_patient_metric.png'), metric_rna)