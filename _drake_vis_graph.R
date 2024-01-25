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
specialLib <- c("LX093_LX094_an_163")
specialLibInd <- grep(specialLib, metadata$name)
nospecialMet <- metadata[-specialLibInd,]
test_data <- nospecialMet[c(1:3,5:6),] # to test code 
lbLst <- unique(test_data$name)
idLst <- lbLst %>% map(splitName)
# get multiplex libraries list 
mulLib <- unique(test_data$name[nchar(test_data$souporcell_link) > 0])
mulId <- mulLib %>% map(splitName)
# get single libraries list 
sngLib <- unique(test_data$name[nchar(test_data$souporcell_link) == 0])
snglId <- sngLib %>% map(splitName)
# get libraries that only demultiplex 
# with souporcell
soclId <- specialLib %>% map(splitName)
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
  psR_specialLib = plot_singler(sR_specialLib, atacSrGeA_specialLib, save_path=atacCellSngRFigDir),
  ## prep for infercnv ---
  preInfer_specialLib = make_anno_count_mx(sR_specialLib,atacSrGeA_specialLib, save_path=AtacInferInputDir),
  ## samples demultiplex with only souporcell ---
  h5Link_specialLib = get_h5_link(specialLib, metadata),
  gexSr_specialLib = create_GEX_seurat(h5Link_specialLib),
  gexSrstress_specialLib = remove_stress_genes(gexSr_specialLib, stress_gene_list),
  gexNor_specialLib = normalize_dim_sr(gexSrstress_specialLib),
  gexSoc_specialLib = assign_metadata_from_souporcell(gexNor_specialLib, metadata),
  vis_demul_specialLib = visualize_soc_demulti(gexSoc_specialLib, atacProcessFigDir),
  gexNoDb_specialLib = remove_soc_db_unknown(gexSoc_specialLib),
  gexSID_specialLib = addMetaSoc(gexNoDb_specialLib),
  atacDemul_specialLib = addMetSocAtac(gexSID_specialLib,atacSrDim_specialLib),
  atacMeta_specialLib = addMetaFromFile(ori_metadata, atacDemul_specialLib)
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

  atacDemul = target(addMetAtac(gexSID,atacSrDim),
                        transform = map(gexSID,atacSrDim,
                        id.var = !!mulId,
                        .id = id.var)),
  atacMeta = target(addMetaFromFile(ori_metadata, atacDemul),
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
  ## prep for infercnv ---
  preInfersg = target(make_anno_count_mx(sRsg,atacSrGeAsg, save_path=AtacInferInputDir),
                        transform = map(sRsg,atacSrGeAsg,
                        id_var = !!snglId,
                        .id = id_var)),
  atacMetasg = target(addMetaFromFile(ori_metadata, atacSrDimsg),
                    transform = map(atacSrDimsg,
                    id.var = !!snglId,
                    .id = id.var)),
  ## merge -----

  mrgAtac = target(merge(x =  atacDemul_specialLib, y= c(atacMeta,atacMetasg), add.cell.ids = c(specialLib,lbLst, sngLib)),
                   transform = combine(atacMeta,atacMetasg,
                                       id.var = !!idLst,
                                       .id = id.var)),
  mrgAtacNor = target(sc_atac_normalize(mrgAtac)),
  mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor))
  )

plan <- bind_plans(process_special_lib_plan, combine_peak_plan, process_plan)
# options(clustermq.scheduler = "multicore") # nolint
# make(plan, parallelism = "clustermq", jobs = 2, lock_cache = FALSE)
vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'cleancode_pipeline.png', font_size = 20 )

# test if res exceeb byte is because of drake
# message('load drake atac')
# loadd("atacSrDim_LX049")                        
# loadd("atacSrDim_LX051")                        
# loadd("atacSrDim_LX065") 

# message('merging atac')
# mergetest <- merge(x = atacSrDim_LX049, y =c(atacSrDim_LX051, atacSrDim_LX065 ), add.cell.ids = lbLst, merge.data = TRUE)
# message('finish merging atac without error')