setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

## read metadata ----
filename <- '30012024_remove_relapse.csv'
metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
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

## define plan ----

combine_peak_plan <- drake_plan(
  grFilenoR = target(makeGrFile(metadata, lb),
                  transform = map(lb = !!lbLst, 
                                  id_var = !!idLst,
                                  .id = id_var)),
  allpeaksnoR = target(disjoin(c(grFilenoR)),
                    transform = combine(grFilenoR,
                                        id.var = !!idLst,
                                        .id = id_var)),
  allpeaksnoRFil = filterBadPeaks(allpeaksnoR),
  allpeaksnoRFilChr = filterChrBlacklist(allpeaksnoRFil), 
)


process_special_lib_plan <- drake_plan(
  ## process ---
  atacSrnoR_specialLib = create_atacSr_w_disjoin(specialLib, metadata, allpeaksnoRFilChr, hg38),
  atacSrnoRMe_specialLib = calculate_metrics(atacSrnoR_specialLib, metadata),
  atacSrnoRFil_specialLib = sc_atac_filter_outliers(atacSrnoRMe_specialLib, figSavePath = atacnoRProcessFigDir),
  atacSrnoRNor_specialLib = sc_atac_normalize(atacSrnoRFil_specialLib),
  atacSrnoRDim_specialLib = sc_atac_dim_redu(atacSrnoRNor_specialLib),
  atacSrnoRGeA_specialLib = get_gene_activity(atacSrnoRDim_specialLib),
  ## run singler ---
  sRnoR_specialLib = run_singleR(atacSrnoRGeA_specialLib),
  psRnoR_specialLib = plot_singler(sRnoR_specialLib, atacSrnoRGeA_specialLib, save_path=atacnoRCellSngRFigDir),
  atacnoRSgR_specialLib = get_sgR_label(sRnoR_specialLib, atacSrnoRDim_specialLib),
  ## prep for infercnv ---
  preInfer_noRspecialLib = make_anno_count_mx(sRnoR_specialLib,atacSrnoRGeA_specialLib, save_path=AtacnoRInferInputDir),
  ## samples demultiplex with only souporcell ---
  h5Link_noRspecialLib = get_h5_link(specialLib, metadata),
  gexSrnoR_specialLib = create_GEX_seurat(h5Link_noRspecialLib),
  gexSrnoRstress_specialLib = remove_stress_genes(gexSrnoR_specialLib, stress_gene_list),
  gexNornoR_specialLib = normalize_dim_sr(gexSrnoRstress_specialLib),
  gexSocnoR_specialLib = assign_metadata_from_souporcell(gexNornoR_specialLib, metadata),
  vis_demul_noRspecialLib = visualize_soc_demulti(gexSocnoR_specialLib, atacnoRProcessFigDir),
  gexnoRNoDbnoR_specialLib = remove_soc_db_unknown(gexSocnoR_specialLib),
  gexSIDnoR_specialLib = addMetaSoc(metadata, gexnoRNoDbnoR_specialLib),
  atacDemulnoR_specialLib = addMetSocAtac(gexSIDnoR_specialLib,atacnoRSgR_specialLib),
  atacMetanoR_specialLib = addMetaFromFile(metadata, atacDemulnoR_specialLib)
)

process_plan <- drake_plan(
  hg38 = getHg38Annotation(),
  # multiplex sample
  ## process ---
  atacSrnoR = target(create_atacSr_w_disjoin(lb, metadata, allpeaksnoRFilChr, hg38),
                  transform = map(lb = !!mulLib,
                                  id_var = !!mulId,
                                  .id = id_var)),
  atacSrnoRMe = target(calculate_metrics(atacSrnoR, metadata),
                    transform = map(atacSrnoR,
                                    id_var = !!mulId,
                                    .id = id_var)),
  atacSrnoRFil = target(sc_atac_filter_outliers(atacSrnoRMe, 
                                             figSavePath = atacnoRProcessFigDir),
                     transform = map(atacSrnoRMe,
                                     id_var = !!mulId,
                                     .id = id_var)),
  atacSrnoRNor = target(sc_atac_normalize(atacSrnoRFil),
                     transform = map(atacSrnoRFil,
                                     id.var = !!mulId,
                                     .id = id_var)),
  atacSrnoRDim = target(sc_atac_dim_redu(atacSrnoRNor),
                     transform = map(atacSrnoRNor,
                                     id.var = !!mulId,
                                     .id = id_var)),
  atacSrnoRGeA = target(get_gene_activity(atacSrnoRDim),
                     transform = map(atacSrnoRDim,
                                     id.var = !!mulId,
                                     .id = id_var)),
  ## run singler ---
  sRnoR = target(run_singleR(atacSrnoRGeA),
              transform = map(atacSrnoRGeA,
                              id_var = !!mulId,
                              .id = id_var)),
  psRnoR = target(plot_singler(sRnoR, atacSrnoRGeA, save_path=atacnoRCellSngRFigDir),
               transform = map(sRnoR,atacSrnoRGeA,
                               id_var = !!mulId,
                               .id = id_var)),
  atacSgRnoR = target(get_sgR_label(sRnoR, atacSrnoRDim),
                   transform = map(sRnoR,atacSrnoRDim,
                                   id_var = !!mulId,
                                   .id = id_var)),
  ## prep for infercnv ---
  preInfernoR = target(make_anno_count_mx(sRnoR,atacSrnoRGeA, save_path=AtacnoRInferInputDir),
                    transform = map(sRnoR,atacSrnoRGeA,
                                    id_var = !!mulId,
                                    .id = id_var)),
  
  ## demultiplex -- 
  h5LinknoR = target(get_h5_link(lb, metadata),
                  transform = map(lb = !!mulLib,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexSrnoR = target(create_GEX_seurat(h5LinknoR),
                 transform = map(h5LinknoR,
                                 id.var = !!mulId,
                                 .id = id.var)),
  gexSrnoRstress = target(remove_stress_genes(gexSrnoR, stress_gene_list),
                       transform = map(gexSrnoR,
                                       id.var = !!mulId,
                                       .id = id.var)),
  gexNornoR = target(normalize_dim_sr(gexSrnoRstress),
                  transform = map(gexSrnoRstress,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexSexnoR = target(calculate_sex_genes_module(gexNornoR),
                  transform = map(gexNornoR,
                                  id.var = !!mulId,
                                  .id = id.var)),
  gexSexnoRMeta = target(add_sex_metadata(gexSexnoR),
                      transform = map(gexSexnoR,
                                      id.var = !!mulId,
                                      .id = id.var)),
 gexSocnoR = target(assign_metadata_from_souporcell(gexSexnoRMeta, metadata),
                  transform = map(gexSexnoRMeta,
                                  id.var = !!mulId,
                                  .id = id.var)),
  vis_demul = target(visualize_soc_gender_demulti(gexSocnoR, 
                                                  atacnoRProcessFigDir),
                     transform = map(gexSocnoR,id.var = !!mulId,
                                     .id = id.var)),
  gexNoDbnoR = target(remove_doublet_n_unknown(gexSocnoR),
                   transform = map(gexSocnoR,
                                   id.var = !!mulId,
                                   .id = id.var)),
  gexNoUnknownnoR = target(fix_unknown_gender_or_souporcell(gexNoDbnoR),
                        transform = map(gexNoDbnoR,
                                        id.var = !!mulId,
                                        .id = id.var)),
  gexNoContranoR = target(remove_contradict_cells(gexNoUnknownnoR),
                       transform = map(gexNoUnknownnoR,
                                       id.var = !!mulId,
                                       .id = id.var)),
  gexSIDnoR = target(assign_sample_name_and_tumor_type(metadata, gexNoContranoR),
                  transform = map(gexNoContranoR,
                                  id.var = !!mulId,
                                  .id = id.var)),
  
  atacDemulnoR = target(addMetAtac(gexSIDnoR,atacSgRnoR),
                     transform = map(gexSIDnoR,atacSgRnoR,
                                     id.var = !!mulId,
                                     .id = id.var)),
  atacMetanoR = target(addMetaFromFile(metadata, atacDemulnoR),
                    transform = map(atacDemulnoR,
                                    id.var = !!mulId,
                                    .id = id.var)),
  # single libraries ----
  ## process ---
  atacSrnoRsg = target(create_atacSr_w_disjoin(lb, metadata, allpeaksnoRFilChr, hg38),
                    transform = map(lb = !!sngLib,
                                    id_var = !!snglId,
                                    .id = id_var)),
  atacSrnoRMesg = target(calculate_metrics(atacSrnoRsg, metadata),
                      transform = map(atacSrnoRsg,
                                      id_var = !!snglId,
                                      .id = id_var)),
  atacSrnoRFilsg = target(sc_atac_filter_outliers(atacSrnoRMesg, 
                   figSavePath = atacnoRProcessFigDir),
                       transform = map(atacSrnoRMesg,
                                       id_var = !!snglId,
                                       .id = id_var)),
  atacSrnoRNorsg = target(sc_atac_normalize(atacSrnoRFilsg),
                       transform = map(atacSrnoRFilsg,
                                       id.var = !!snglId,
                                       .id = id_var)),
  atacSrnoRDimsg = target(sc_atac_dim_redu(atacSrnoRNorsg),
                       transform = map(atacSrnoRNorsg,
                                       id.var = !!snglId,
                                       .id = id_var)),
  atacSrnoRGeAsg = target(get_gene_activity(atacSrnoRDimsg),
                       transform = map(atacSrnoRDimsg,
                                       id.var = !!snglId,
                                       .id = id_var)),
  ## run singler ---
  sRsgnoR = target(run_singleR(atacSrnoRGeAsg),
                transform = map(atacSrnoRGeAsg,
                                id_var = !!snglId,
                                .id = id_var)),
  psRnoRsg = target(plot_singler(sRsgnoR, atacSrnoRGeAsg, save_path=atacnoRCellSngRFigDir),
                 transform = map(sRsgnoR,atacSrnoRGeAsg,
                                 id_var = !!snglId,
                                 .id = id_var)),
  atacsgSgRnoR = target(get_sgR_label(sRsgnoR, atacSrnoRDimsg),
                     transform = map(sRsgnoR,atacSrnoRDimsg,
                                     id_var = !!snglId,
                                     .id = id_var)),
  ## prep for infercnv ---
  preInfernoRsg = target(make_anno_count_mx(sRsgnoR,atacSrnoRGeAsg, save_path=AtacnoRInferInputDir),
                      transform = map(sRsgnoR,atacSrnoRGeAsg,
                                      id_var = !!snglId,
                                      .id = id_var)),
  ## add metadata ----
  atacMetasgnoR = target(addMetaFromFile(metadata, atacsgSgRnoR),
                      transform = map(atacsgSgRnoR,
                                      id.var = !!snglId,
                                      .id = id.var)),
  ## merge -----
  mrgAtacnoR = target(merge_pairwise(c(atacMetanoR_specialLib, atacMetanoR, atacMetasgnoR),atcnoRMrgDir),
            transform = combine(atacMetanoR,atacMetasgnoR,
                    id.var = !!c(mulLib, sngLib),
                    .id = id.var)),
  mrgAtacnoRNor = target(sc_atac_normalize(mrgAtacnoR)),
  mrgAtacnoRDim = target(sc_atac_dim_redu(mrgAtacnoRNor)),
  mrgnoRPtype = dimplot_w_nCell_label(mrgAtacnoRDim, by = 'Subtype',atacnoRMrgFigDir , col = my_cols2),
  mrgnoRPsID = dimplot_w_nCell_label(mrgAtacnoRDim, by = 'sampleID',atacnoRMrgFigDir , col = my_cols2),
  mrgnoRPlibnCell = dimplot_w_nCell_label(mrgAtacnoRDim, by = 'library',atacnoRMrgFigDir, col = my_cols2),
  mrgnoRP = dimplotnSave(mrgAtacnoRDim, atacnoRMrgFigDir, save_name = 'cluster'),
  mrgnoRPlib = dimplotBynSave(mrgAtacnoRDim, by = 'library',
                           atacnoRMrgFigDir, save_name = 'lib_no_nCell', col = my_cols2),
  # prep for infercnv
  mrgnoRGA = get_gene_activity(mrgAtacnoRDim),
  preInfernoRMrg = make_anno_count_Mrgmx(mrgnoRGA, save_path=AtacnoRInferInputDir),

  # dimplot each sample
  # dimP = target(dimplotnSave(sr, atacnoRMrgFigDir, save_name = 'cluster'),add.cell.ids = ),
  #                  transform = map(atacMeta,atacMetasg,
  #                 id.var = !!c(specialLib,mulLib,sngLib),
  
  # integration w anchors
  # anchors = target(FindIntegrationAnchors(c(atacMeta_specialLib, c(atacMeta,atacMetasg)), 
  #                  reduction = 'rlsi'),
  #                 transform = combine(atacMeta,atacMetasg,
  #                     id.var = !!c(mulLib, sngLib),
  #                     .id = id.)),
  # srInt = IntegrateData(anchorset = anchors, dims = 1:50)
)

plan <- bind_plans(combine_peak_plan,process_special_lib_plan,  process_plan)
# options(clustermq.scheduler = "multicore") # nolint
# make(plan, parallelism = "clustermq", jobs = 1, lock_cache = FALSE)
# make(plan, lock_cache = FALSE)
vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'cleancode_no_relapse_pipeline.png', font_size = 20 )

