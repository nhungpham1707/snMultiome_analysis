setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
lapply(list.files("./R", full.names = TRUE), source)

## read metadata ----
filename <- '30_11_2023_all_multiome_libraries.csv' 
metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
specialLib <- c("LX093_LX094_an_163")
specialLibInd <- grep(specialLib, metadata$name)
metadata <- metadata[-specialLibInd,]
test_data <- metadata[c(1:3,5:6),] # to test code 
lbLst <- unique(test_data$name)
idLst <- lbLst %>% map(splitName)
mulLib <- unique(test_data$name[nchar(test_data$souporcell_link) > 0])
mulId <- mulLib %>% map(splitName)

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

process_plan <- drake_plan(
  hg38 = getHg38Annotation(),
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
  sR = target(run_singleR(atacSrGeA),
                        transform = map(atacSrGeA,
                        id_var = !!mulId,
                        .id = id_var)),
  psR = target(plot_singler(sR, atacSrGeA, save_path=atacCellSngRFigDir),
                        transform = map(sR,atacSrGeA,
                        id_var = !!mulId,
                        .id = id_var)),
  preInfer = target(make_anno_count_mx(sR,atacSrGeA, save_path=AtacInferInputDir),
                        transform = map(sR,atacSrGeA,
                        id_var = !!mulId,
                        .id = id_var)),
  
  # demultiplex ----
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
                        .id = id.var))
  
  # mrgAtac = target(merge(c(atacSrDim), add.cell.ids = lbLst),
  #                  transform = combine(atacSrDim,
  #                                      id.var = !!idLst,
  #                                      .id = id.var)),
  # mrgAtacNor = target(sc_atac_normalize(mrgAtac)),
  # mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor))
  )

# plan <- bind_plans(combine_peak_plan, process_plan, sample_demultiplex_plan)
plan <- bind_plans(combine_peak_plan, process_plan)
# options(clustermq.scheduler = "multicore") # nolint
# make(plan, parallelism = "clustermq", jobs = 2, lock_cache = FALSE)
make(plan,lock_envir = TRUE, lock_cache = FALSE, verbose = 0)

vis_drake_graph(plan, targets_only = TRUE, lock_cache = FALSE, file = 'cleancode_pipeline.png', font_size = 20 )

# test if res exceeb byte is because of drake
# message('load drake atac')
# loadd("atacSrDim_LX049")                        
# loadd("atacSrDim_LX051")                        
# loadd("atacSrDim_LX065") 

# message('merging atac')
# mergetest <- merge(x = atacSrDim_LX049, y =c(atacSrDim_LX051, atacSrDim_LX065 ), add.cell.ids = lbLst, merge.data = TRUE)
# message('finish merging atac without error')