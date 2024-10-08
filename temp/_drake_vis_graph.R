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
  # mrgPlibnCell = dimplot_w_nCell_label(mrgAtacDim, by = 'library',atacMrgFigDir , col = my_cols2),
  # mrgP = dimplotnSave(mrgAtacDim, atacMrgFigDir, save_name = 'cluster'),
  # mrgPlib = dimplotBynSave(mrgAtacDim, by = 'library',
  #                          atacMrgFigDir, save_name = 'lib_no_nCell', col = my_cols2),
  # # prep mrg atac for infercnv
  mrgGA = get_gene_activity(mrgAtacDim),
  # preInferMrg = make_anno_count_Mrgmx(mrgGA, save_path=AtacInferInputDir),
  # integration w anchors -----------------
  # atac_anchors = target(FindIntegrationAnchors(c(atacMeta_specialLib, atacMeta, atacMetasg), reduction = 'rlsi', anchor.features = 2000),
  #                   transform = combine(atacMeta,atacMetasg,
  #                   id.var = !!c(mulLib, sngLib),
  #  atac_anchors = target(integration_w_anchors_subset(c(atacMeta_specialLib, atacMeta, atacMetasg), to_remove = 'atac'),
  #                   transform = combine(atacMeta,atacMetasg,
  #                   id.var = !!c(mulLib, sngLib),
  #                   .id = id.var)),
  # atac.integrated = IntegrateData(anchorset = atac_anchors, dims = 1:50),
  # int.AtacNor = normalize_anchors(atac.integrated),
  # int.AtacDim = sc_atac_dim_redu(int.AtacNor),
  # saveInt = saveRDS(int.AtacDim, paste0(atcMrgDir, '/int.AtacDim.RDS')),
  # int.Ptype = dimplotnSave (int.AtacDim, atacMrgFigDir, save_name = 'int.atac.png', col = my_cols2),
  # int.sID = dimplotBynSave(int.AtacDim, by = 'sampleID',atacMrgFigDir, save_name = 'int.atac.sID.png', col = my_cols2),
  # int.subtype = dimplotBynSave(int.AtacDim, by = 'Subtype',atacMrgFigDir, save_name = 'int.atac.subtype.png', col = my_cols2),
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
  p_mrgRna_all_0.5 = plot_cluster(mrgRnaClu, rnaMrgFigDir, group.by = 'rna_snn_res.0.5', save_name = 'merge' ),
   saveplotmrgrna_0.5 = savePlot(paste0(rnaMrgFigDir, '/mrgRNA.png'), p_mrgRna_all_0.5),
   p_mrgRna_all_0.3 = plot_cluster(mrgRnaClu, rnaMrgFigDir, group.by = 'rna_snn_res.0.3', save_name = 'merge' ),
   saveplotmrgrna_0.3 = savePlot(paste0(rnaMrgFigDir, '/mrgRNA.png'), p_mrgRna_all_0.3),
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
  p_mrgRna_noNOr_all_0.5 = plot_cluster(mrgRnaClu_noNOr, rnaMrgFigDir, group.by = 'rna_snn_res.0.5', save_name = 'merge_noNOr' ),
   saveplotmrgrna_noNOr_0.5 = savePlot(paste0(rnaMrgFigDir, '/mrgRNA_noNOr.png'), p_mrgRna_noNOr_all_0.5),
   p_mrgRna_noNOr_all_0.3 = plot_cluster(mrgRnaClu_noNOr, rnaMrgFigDir, group.by = 'rna_snn_res.0.3', save_name = 'merge_noNOr' ),
   saveplotmrgrna_noNOr_0.3 = savePlot(paste0(rnaMrgFigDir, '/mrgRNA_noNOr.png'), p_mrgRna_noNOr_all_0.3),
   p_mrgRna_noNOr_all = plot_cluster(mrgRnaClu_noNOr, rnaMrgFigDir, save_name = 'merge_noNOr' ),
   saveplotmrgrna_noNOr = savePlot(paste0(rnaMrgFigDir, '/mrgRNA_noNOr.png'), p_mrgRna_noNOr_all)

)

plan <- bind_plans(combine_peak_plan,process_special_lib_plan,  process_plan)

vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_cleancode_pipeline.png', font_size = 20 )

# test if res exceeb byte is because of drake
# message('load drake atac')
# loadd("atacSrDim_LX049")                        
# loadd("atacSrDim_LX051")                        
# loadd("atacSrDim_LX065") 

# message('merging atac')
# mergetest <- merge(x = atacSrDim_LX049, y =c(atacSrDim_LX051, atacSrDim_LX065 ), add.cell.ids = lbLst, merge.data = TRUE)
# message('finish merging atac without error')

