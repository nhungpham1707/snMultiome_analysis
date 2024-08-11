# the project is run with drake, 
# a workflow manager for R. Refer to 
# https://github.com/ropensci/drake for more detail

message('---------- Start running _drake.R ------------')
setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
source("./global_variables.R")
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
# with souporcell (because they contain samples
# with the same gender)
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

# read metadata for healthy data that are used for identifying cell origin ----
hthy_dataDir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes'
healthy_metadata <- read.csv(paste0(hthy_dataDir,'/filtered.cell_metadata.for_website.txt.gz'), sep = '\t')
all_hthytissue_list <- unique(healthy_metadata$tissue)
hthytissue_list <- all_hthytissue_list

# define drake plan ----
# ----------------------------------------------------------
# plan to generate a common peak file to merge all atac samples with disjoin.
# check here for more info https://stuartlab.org/signac/articles/merging
# ----------------------------------------------------------

## generate common peak ----

combine_peak_plan <- drake_plan(
  grFile = target(makeGrFile(metadata, lb),
                  transform = map(lb = !!lbLst, 
                                  id.vars = !!idLst,
                                  .id = id.vars)),
  allpeaks = target(disjoin(c(grFile)),
                    transform = combine(grFile,
                                        id.vars = !!idLst,
                                        .id = id.vars)),
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
# - prep for infercnv: generate gene matrix 
# and annotation file for infercnv 
# - sample demultiplex using souporcell
# - assign metadata
# ----------------------------------------------------------

# process special lib ----
process_special_lib_plan <- drake_plan(
  ## process ----
  atacSr_specialLib = create_atacSr_w_disjoin(specialLib, metadata, 
                                              allpeaksFilChr, hg38),
  atacSrMe_specialLib = calculate_metrics(atacSr_specialLib, metadata),
  atacSrFil_specialLib = sc_atac_filter_outliers(atacSrMe_specialLib, 
                                                 figSavePath = atacProcessFigDir),
  atacSrNor_specialLib = sc_atac_normalize(atacSrFil_specialLib),
  atacSrDim_specialLib = sc_atac_dim_redu(atacSrNor_specialLib),

  atacSrGeA_specialLib = get_gene_activity(atacSrDim_specialLib),

  ## run singler ----
  sR_specialLib = run_singleR(atacSrGeA_specialLib),
  psR_specialLib = plot_singler(sR_specialLib, atacSrGeA_specialLib, 
                                save_path=atacCellSngRFigDir),
  atacGASgR_specialLib = get_sgR_label(sR_specialLib, atacSrGeA_specialLib),
  atacSgR_specialLib = get_sgR_label(sR_specialLib, atacSrDim_specialLib),

  ## prep input data for infercnv ----
  preInfer_specialLib = make_anno_count_mx(atacGASgR_specialLib, 
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

# ----------------------------------------------------------
# process all libraries. multiplex scATCseq libraries 
# have additional steps to demultiplex using 
# souporcell and gender from scRNAseq data
# ----------------------------------------------------------

# process all libraries ----
process_plan <- drake_plan(
  hg38 = getHg38Annotation(),
  ## multiplex atacseq libraries -------------
  ### process ----
  atacSr = target(create_atacSr_w_disjoin(lb, metadata, allpeaksFilChr, hg38),
                  transform = map(lb = !!mulLib,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  atacSrMe = target(calculate_metrics(atacSr, metadata),
                    transform = map(atacSr,
                                    id.vars = !!mulId,
                                    .id = id.vars)),
  atacSrFil = target(sc_atac_filter_outliers(atacSrMe, 
                     figSavePath = atacProcessFigDir),
                     transform = map(atacSrMe,
                                     id.vars = !!mulId,
                                     .id = id.vars)),
  atacSrNor = target(sc_atac_normalize(atacSrFil),
                     transform = map(atacSrFil,
                                     id.vars = !!mulId,
                                     .id = id.vars)),
  atacSrDim = target(sc_atac_dim_redu(atacSrNor),
                     transform = map(atacSrNor,
                                     id.vars = !!mulId,
                                     .id = id.vars)),
  atacSrGeA = target(get_gene_activity(atacSrDim),
                     transform = map(atacSrDim,
                                     id.vars = !!mulId,
                                     .id = id.vars)),
  ### run singler ----
  sR = target(run_singleR(atacSrGeA),
              transform = map(atacSrGeA,
                              id.vars = !!mulId,
                              .id = id.vars)),
  psR = target(plot_singler(sR, atacSrGeA, save_path=atacCellSngRFigDir),
               transform = map(sR,atacSrGeA,
                               id.vars = !!mulId,
                               .id = id.vars)),
  atacSgR = target(get_sgR_label(sR, atacSrDim),
                   transform = map(sR,atacSrDim,
                                   id.vars = !!mulId,
                                   .id = id.vars)),
  atacGASgR = target(get_sgR_label(sR, atacSrGeA),
                   transform = map(sR,atacSrGeA,
                                   id.vars = !!mulId,
                                  .id = id.vars)),                          
  ### prep for infercnv ----
  preInfer = target(make_anno_count_mx(atacGASgR, save_path=AtacInferInputDir),
                    transform = map(atacGASgR,
                                    id.vars = !!mulId,
                                    .id = id.vars)),
  ### demultiplex ----
  h5Link = target(get_h5_link(lb, metadata),
                  transform = map(lb = !!mulLib,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  gexSr = target(create_GEX_seurat(h5Link),
                 transform = map(h5Link,
                                 id.vars = !!mulId,
                                 .id = id.vars)),
  gexSrstress = target(remove_stress_genes(gexSr, stress_gene_list),
                       transform = map(gexSr,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  gexNor = target(normalize_dim_sr(gexSrstress),
                  transform = map(gexSrstress,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  gexSex = target(calculate_sex_genes_module(gexNor),
                  transform = map(gexNor,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  gexSexMeta = target(add_sex_metadata(gexSex),
                      transform = map(gexSex,
                                      id.vars = !!mulId,
                                      .id = id.vars)),
  gexSoc = target(assign_metadata_from_souporcell(gexSexMeta, metadata),
                  transform = map(gexSexMeta,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  vis_demul = target(visualize_soc_gender_demulti(gexSoc, 
                                                  atacProcessFigDir),
                     transform = map(gexSoc,id.vars = !!mulId,
                                     .id = id.vars)),
  gexNoDb = target(remove_doublet_n_unknown(gexSoc),
                   transform = map(gexSoc,
                                   id.vars = !!mulId,
                                   .id = id.vars)),
  gexNoUnknown = target(fix_unknown_gender_or_souporcell(gexNoDb),
                        transform = map(gexNoDb,
                                        id.vars = !!mulId,
                                        .id = id.vars)),
  gexNoContra = target(remove_contradict_cells(gexNoUnknown),
                       transform = map(gexNoUnknown,
                                       id.vars = !!mulId,
                                       .id = id.vars)),
  gexSID = target(assign_sample_name_and_tumor_type(metadata, gexNoContra),
                  transform = map(gexNoContra,
                                  id.vars = !!mulId,
                                  .id = id.vars)),
  ### prepare demultiplex metadata ----
  gexDemulMeta = target(generate_demultiplex_metadata(c(gexSID)),
                  transform = combine(gexSID,
                  id.vars = !!mulId,
                  .id = id.vars)),
  atacDemul = target(addMetAtac(gexSID,atacSgR),
                     transform = map(gexSID,atacSgR,
                                     id.vars = !!mulId,
                                     .id = id.vars)),
  atacMeta = target(addMetaFromFile(metadata, atacDemul),
                    transform = map(atacDemul,
                                    id.vars = !!mulId,
                                    .id = id.vars)),

  ## single atacseq libraries ----
  ### process ----
  atacSrsg = target(create_atacSr_w_disjoin(lb, metadata, allpeaksFilChr, hg38),
                    transform = map(lb = !!sngLib,
                                    id.vars = !!snglId,
                                    .id = id.vars)),
  atacSrMesg = target(calculate_metrics(atacSrsg, metadata),
                      transform = map(atacSrsg,
                                      id.vars = !!snglId,
                                      .id = id.vars)),
  atacSrFilsg = target(sc_atac_filter_outliers(atacSrMesg, 
                   figSavePath = atacProcessFigDir),
                       transform = map(atacSrMesg,
                                       id.vars = !!snglId,
                                       .id = id.vars)),
  atacSrNorsg = target(sc_atac_normalize(atacSrFilsg),
                       transform = map(atacSrFilsg,
                                       id.vars = !!snglId,
                                       .id = id.vars)),
  atacSrDimsg = target(sc_atac_dim_redu(atacSrNorsg),
                       transform = map(atacSrNorsg,
                                       id.vars = !!snglId,
                                       .id = id.vars)),
  atacSrGeAsg = target(get_gene_activity(atacSrDimsg),
                       transform = map(atacSrDimsg,
                                       id.vars = !!snglId,
                                       .id = id.vars)),
  ### run singler ----
  atac_sRsg = target(run_singleR(atacSrGeAsg),
                transform = map(atacSrGeAsg,
                                id.vars = !!snglId,
                                .id = id.vars)),
  p_atac_sRsg = target(plot_singler(atac_sRsg, atacSrGeAsg, save_path=atacCellSngRFigDir),
                 transform = map(atac_sRsg,atacSrGeAsg,
                                 id.vars = !!snglId,
                                 .id = id.vars)),
  atacsgSgR = target(get_sgR_label(atac_sRsg, atacSrDimsg),
                     transform = map(atac_sRsg,atacSrDimsg,
                                     id.vars = !!snglId,
                                     .id = id.vars)),
  atacsgGASgR = target(get_sgR_label(atac_sRsg, atacSrGeAsg),
                    transform = map(atac_sRsg,atacSrGeAsg,
                              id.vars = !!snglId,
                              .id = id.vars)),
  ### prep for infercnv ----
  preInfersg = target(make_anno_count_mx(atacsgGASgR, save_path=AtacInferInputDir),
                      transform = map(atacsgGASgR,
                                      id.vars = !!snglId,
                                      .id = id.vars)),
  ### add metadata ----
  atacMetasg = target(addMetaFromFile(metadata, atacsgSgR),
                      transform = map(atacsgSgR,
                                      id.vars = !!snglId,
                                      .id = id.vars)),
  ## merge all atac-----
  mrgAtac = target(merge_pairwise(c(atacMeta_specialLib, atacMeta, atacMetasg),atcMrgDir),
            transform = combine(atacMeta,atacMetasg,
                    id.vars = !!c(mulLib, sngLib),
                    .id = id.vars)),
  mrgAtacNor = target(sc_atac_normalize(mrgAtac)),
  mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor)),
  mrgPtype = dimplot_w_nCell_label(mrgAtacDim, by = 'Subtype',atacMrgFigDir , col = my_cols2),
  mrgPsID = dimplot_w_nCell_label(mrgAtacDim, by = 'sampleID',atacMrgFigDir , col = my_cols2),

  atac_noNA = remove_na_cells(mrgAtacDim),

  ## group atac singleR cell types ----
  atac_group_sgr = group_singleR_labels(atac_noNA),

  ## prep mrg atac for infercnv ----
  atacGA = get_gene_activity(atac_group_sgr),

  ## process rna ------------------------------------------
  h5Link_all = target(get_h5_link(lb, metadata),
                  transform = map(lb = !!alLib,
                                  id.vars = !!alID,
                                  .id = id.vars)),
  gexSr_all = target(create_GEX_seurat(h5Link_all),
                  transform = map(h5Link_all,
                              id.vars = !!alID,
                              .id = id.vars)),
  gexSrstress_all = target(remove_stress_genes(gexSr_all,
                          stress_gene_list),
                          transform = map(gexSr_all,
                              id.vars = !!alID,
                              .id = id.vars)),
  gexMito_all = target(check_mito_genes(gexSrstress_all),
                    transform = map(gexSrstress_all,
                    id.vars = !!alID,
                    .id = id.vars)),
  p_rnametric = target(rna_visualize_metric(gexMito_all,rnaFigDir),
                    transform = map(gexMito_all,
                    id.vars = !!alID,
                    .id = id.vars)),
  gexFilt = target(filter_rna_sr_w_isoutlier(gexMito_all, rnaFigDir),
                  transform = map(gexMito_all,
                  id.vars = !!alID,
                  .id = id.vars)),
  gexNor_all = target(normalize_dim_plot_sr(gexFilt, rnaFigDir, alLib),
                          transform = map(gexFilt,!!alLib,
                          id.vars = !!alID,
                              .id = id.vars)),
  gexClus = target(clustering_rna_data(gexNor_all),
                  transform = map(gexNor_all,
                  id.vars = !!alID,
                              .id = id.vars)),
  p_rnaClus = target(plot_cluster(gexClus, rnaFigDir, alLib), 
              transform = map(gexClus, !!alLib,
                  id.vars = !!alID,
                              .id = id.vars)),
  
   ## run singleR rna ----
  rna_sgr = target(run_singleR(gexClus),
                transform = map(gexClus,
                                id.vars = !!alID,
                                .id = id.vars)),
  prna_sgr = target(plot_singler(rna_sgr, gexClus, save_path=rnaCellSngRFigDir),
                 transform = map(rna_sgr,gexClus,
                                 id.vars = !!alID,
                                 .id = id.vars)),
  gexClusSgr = target(get_sgR_label(rna_sgr, gexClus),
                transform = map(rna_sgr, gexClus,
                            id.vars = !!alID,
                            .id = id.vars)),
  ## merge rna ----
   mrgRna_all = target(merge_pairwise(c(gexClusSgr), rnaMrgDir),
            transform = combine(gexClusSgr,
            id.vars = !!alID,
            .id = id.vars)),
  mrgRnaNor_all = normalize_dim_plot_sr(mrgRna_all, rnaMrgFigDir, lib_name = 'merge'),
  mrgRnaClu = clustering_rna_data(mrgRnaNor_all),
  rna_meta = assign_meta(metadata, gexDemulMeta,
  mrgRnaClu, save_name = paste0(rnaMrgDir, '/mrgRna_meta.RDS')),
  rna_noNA = remove_na_cells(rna_meta),
  rna_fix = fix_special_lib_rna(gexSID_specialLib, rna_noNA, gexNoDb_specialLib),
  rna_group_sgr = group_singleR_labels(rna_fix),
  ## rename cell barcodes ----
  new_bc = paste0(rna_group_sgr$barcodes, '_', rna_group_sgr$library),
  rna_new_bc = RenameCells(rna_group_sgr, new.names = new_bc),

  # prep rna for infercnv -----
  # prepare count matrix 
  preInferRna = target(make_anno_count_mx(gexClusSgr, save_path = rnaInferInputDir ),
                    transform = map(gexClusSgr,
                  id.vars = !!alID,
                  .id = id.vars)),
   ## remove confounding genes before running infercnv ----
  gex_noCF = target(subset(gexClusSgr, feature = gene_to_retain),
            transform = map(gexClusSgr,
                          id.vars = !!alID,
                          .id = id.vars)),
  preInferRna_noCF = target(make_anno_count_mx(gex_noCF, save_path = rnaInferInputDir ),
                    transform = map(gex_noCF,
                  id.vars = !!alID,
                  .id = id.vars)),
  
  # reanalyze unknown and unassign cells from souporcells and gender demultiplex ----
  # since they look like they could be doublets 
  unknown_cells = target(get_unknown_unassign_cells(gexNoDb),   transform = map(gexNoDb,
            id.vars = !!mulId,
            .id = id.vars)),
  unknown_bcs = target(get_bc_in_mrg(rna_group_sgr, unknown_cells, lib = mulId),
          transform = map(unknown_cells, !!mulId,
          id.vars = !!mulId,
          .id = id.vars)),
  # all_unknown is generated manually by loadding all unknown_bcs and combine in a list 

  
)

# inspect cluster quality -----
cluster_behavior_plan <- drake_plan( 
  # make atac sce ----
  atac.sce = make_sce(atac_group_sgr),
  # # make rna sce -----
  rna.sce = make_sce(rna_new_bc),
  # # check cluster behavior ---- 
  # ## Silhouette widths ------
  # ### atac ----
  sil.atac = calculate_silhouette(atac.sce, reduce_method = 'LSI'),
  sil.atac.p =  ggplot(sil.atac, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley")
  + scale_colour_manual(values = my_cols) +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_silAtac_p = savePlot('output/batchEffect/atac_cluster_behavior.png', sil.atac.p),
  silAtac_tab = table(Cluster=colLabels(atac.sce), sil.atac$closest),  # cluster 0, 1 and 4
  # have many cells that can easily mix with other clusters
  
  # ### rna -----
  sil.rna = calculate_silhouette(rna.sce, reduce_method = 'PCA'),
  sil.rna.p = ggplot(sil.rna, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_sil.rna.p = savePlot('output/batchEffect/rna_silhouette_cluster_behavior.png', sil.rna.p),
  table(Cluster=colLabels(rna.sce), sil.rna$closest),
  # # cluster 7 & 23 have many cells that can easily mix with other clusters
  # ## cluster purity ------
  # ### atac ----
  pure.atac = calculate_purity(atac.sce, reduce_method = 'LSI'),
  pure.atac.p = ggplot(pure.atac , aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
  scale_colour_manual(values = my_cols) +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_pure.atac.p = savePlot('output/batchEffect/atac_cluster_purity.png', pure.atac.p),
  # # To determine which clusters contaminate each other, we can identify the cluster with the most neighbors for each cell. In the table below, each row corresponds to one cluster; large off-diagonal counts indicate that its cells are easily confused with those from another cluster.
  #  ### rna ----
  pure.rna = calculate_purity(rna.sce, 'PCA'),
  pure.rna.p = ggplot(pure.rna, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_pure.rna.p = savePlot('output/batchEffect/rna_cluster_purity.png', pure.rna.p)
)


# from the UMAP, it seems clusters are driven based on other factors (e.g. patient, library)
# rather than biology. Here we detect which factor drive the cluster
# detect batch effect -----
batch_detection_plan <- drake_plan(
  ## visualize batch ----
  ### atac ----
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
    atac_cms = calculate_n_plot_cms(atac.sce, save_path = batchAtacDir, 
                                    neighbors = 200, save_name = 'dj', 'LSI'),
  ### rna ----
    rna_cms = calculate_n_plot_cms(rna.sce, save_path = batchRnaDir, 
                                   save_name = 'no_correction', neighbors = 200, 'PCA'),

  ## check cell cycle ----
    rna_cellCycle = check_cell_cycle(rna_group_sgr, save_path = batchRnaDir)
)

# ------------------------------------------------------
# it seems clusters are driven by library. Many techniques (harmony, seurat cca, mnn, 
# epiConv, BAVARIA, Scanorama, SysVi) have been tested, 
# among them harmony seems to be the best and easiest one to correct the data 
# harmony was run with various theta (0,0.1,0.2,0.5,1) and batch factors (library, patient)
# to select the best parameters
# ------------------------------------------------------
batch_correction_plan <- drake_plan(
  # harmony ----------------------------------
  final_hm_rna = RunHarmony(rna_group_sgr, group.by.vars = 'library', 
                            theta = 0.8),
  final_hm_rna_nb = FindNeighbors(object = final_hm_rna, reduction = "harmony", 
                                  k.param = 30, dims = 1:30),
  final_hm_rna_clus = FindClusters(final_hm_rna_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  final_hm_rna_umap = RunUMAP(final_hm_rna_clus, dims = 1:30, reduction = 'harmony'),
  hm_rna_p = DimPlot(final_hm_rna_umap, group.by = 'Subtype', cols = my_cols),
  save_hm_rna_p = savePlot(paste0(batchRnaHarmonyDir, '/final_hm_subtype.png'), hm_rna_p),

  final_hm_atac = RunHarmony(atac_group_sgr, group.by.vars = 'library', 
                             theta = 0.2, reduction.use = 'lsi',  assay.use = 'peaks', 
                             project.dim = FALSE),
  final_hm_atac_nb = FindNeighbors(object = final_hm_atac, reduction = "harmony", k.param = 30),
  final_hm_atac_clus = FindClusters(final_hm_atac_nb, resolution = c(0.2,0.4,0.6, 0.8,1)),
  final_hm_atac_umap = RunUMAP(final_hm_atac_clus, dims = 1:30, reduction = 'harmony'),
  hm_atac_p = DimPlot(final_hm_atac_umap, group.by = 'Subtype', cols = my_cols),
  save_hm_atac_p = savePlot(paste0(batchAtacHarmonyDir, '/final_hm_subtype.png'), hm_atac_p),
  
  # remove MHC genes and other confounding genes ----
  # conclusion after removing CF genes: it does not improve the clustering 
  ## rna ----
  genes_to_remove = unique(c(genelists$chr6HLAgenes, genelists$hemo, 
                             genelists$stress, genelists$ribo)), 
  gene_to_retain = setdiff(rownames(rna_group_sgr), genes_to_remove ),
  rna_noCF = subset(rna_group_sgr, feature = gene_to_retain),
  rna_noCF_nor = normalize_dim_plot_sr(rna_noCF, rnaMrgFigDir, lib_name = 'merge_noCF'),
  rna_noCF_nor_clu = clustering_rna_data(rna_noCF_nor),
  rna_noCF_meta = assign_meta(metadata, gexDemulMeta,
  rna_noCF_nor_clu, save_name = paste0(rnaMrgDir, '/mrgRna_noCF_meta.RDS')),
  dim_rna_noCF_lib = DimPlot(rna_noCF_meta, group.by = 'library', 
                             cols = my_cols, raster = FALSE,pt.size = 1),
  save_dim_rna_noCF_lib = savePlot(paste0(rnaMrgFigDir, '/noCF_lib.png'), 
                                   dim_rna_noCF_lib),
  dim_rna_noCF_sub = DimPlot(rna_noCF_meta, group.by = 'Subtype', 
                             cols = my_cols, raster = FALSE,pt.size = 1),
  save_dim_rna_noCF_sub = savePlot(paste0(rnaMrgFigDir, '/noCF_sub.png'), 
                                   dim_rna_noCF_sub)
)

# annotate cell and identify tumor -----
cell_annotation_plan <- drake_plan(
  ## from infercnv ----
  # infercnv was run in a different conda env r43_copy, script _drake_infercnv.R
  ### infercnv res intepretation ----
  # infer_res was generated manually from get_infercnv_result.R
  atac_infer_res = read.csv(paste0(cellAtacInferDir, '/atac_infer_res.csv')),
  atac_infer = assign_infer_res_to_sr(atac_infer_res, final_hm_atac_umap),
 
 rna_infer_res = read.csv(paste0(cellRnaIcnvdir, '/rna_infer_res.csv')),
 rna_infer = assign_infer_res_to_sr(rna_infer_res, final_hm_rna_umap),
 rna_infer_noCF_res = read.csv(paste0(cellRnaIcnvdir, '/rna_remove_cf_infer_res.csv')),
 rna_noCF_infer = assign_infer_res_to_sr(rna_infer_noCF_res, final_hm_rna_umap),
 
 ## scroshi merg rna ----
  hmRna_scroshi_demo = run_scROSHI_w_demo_data(sr = final_hm_rna_umap, 
                                               cols = my_cols, pt = 1, 
                                               save_name = 'hm_rna_w_demo_marker', 
                                               save_path = CellRnaScroshiDir),
  
  hmRna_scroshi_atrt = run_scROSHI_w_cancer_marker(sr = hmRna_scroshi_demo, 
                                                   cols = my_cols, pt = 1, 
                                                   save_name = 'hm_rna_w_atrt', 
                                                   save_path = CellRnaScroshiDir),

  hmRna_infer_nocf_scroshi_demo = run_scROSHI_w_demo_data(sr = rna_noCF_infer, 
                                                          cols = my_cols, pt = 1, 
                                                          save_name = 'hm_rna_infer_nocf_w_demo_marker', 
                                                          save_path = CellRnaScroshiDir),
  
  hmRna_infer_nocf_scroshi_atrt = run_scROSHI_w_cancer_marker(sr = hmRna_infer_nocf_scroshi_demo, 
                                                              cols = my_cols, pt = 1, 
                                                              save_name = 'hm_rna_infer_nocf_w_atrt', 
                                                              save_path = CellRnaScroshiDir),
  ## scroshi merg atac ---
  hm_atacGA = get_gene_activity(final_hm_atac_umap),

  hmAtac_scroshi_demo = run_scROSHI_w_demo_data(sr = hm_atacGA, cols = my_cols, 
                                                pt = 1, save_name = 'hm_atac_w_demo_marker', 
                                                save_path = atacScroshiDir),
  
  hmAtac_scroshi_atrt = run_scROSHI_w_cancer_marker(sr = hm_atacGA, cols = my_cols, 
                                                    pt = 1, save_name = 'hm_atac_w_atrt', 
                                                    save_path = atacScroshiDir) 
  
)


# check if cluster behavior was improved after batch correction 
cluster_behavior_after_correction_plan <- drake_plan( 
  # make atac sce ----
  hmatac.sce = make_sce(final_hm_atac_umap),
  # # make rna sce -----
  hmrna.sce = make_sce(final_hm_rna_umap),
  # # check cluster behavior ---- 
  # ## Silhouette widths ------
  # ### atac ----
  sil.hmatac = calculate_silhouette(hmatac.sce, reduce_method = 'LSI'),
  sil.hmatac.p =  ggplot(sil.hmatac, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley")
  + scale_colour_manual(values = my_cols) +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_silhmAtac_p = savePlot('output/batchEffect/hmatac_cluster_behavior.png', sil.hmatac.p),
  silhmAtac_tab = table(Cluster=colLabels(hmatac.sce), sil.hmatac$closest),  # cluster 0, 1 and 4 have many cells that can easily mix with other clusters
  # ### rna -----
  sil.hmrna = calculate_silhouette(hmrna.sce, reduce_method = 'PCA'),
  sil.hmrna.p = ggplot(sil.hmrna, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_sil.hmrna.p = savePlot('output/batchEffect/hmrna_silhouette_cluster_behavior.png', sil.hmrna.p),
  sil.hmrna_tab = table(Cluster=colLabels(hmrna.sce), sil.hmrna$closest),
  # # cluster 7 & 23 have many cells that can easily mix with other clusters
  # ## cluster purity ------
  # ### atac ----
  pure.hmatac = calculate_purity(hmatac.sce, reduce_method = 'LSI'),
  pure.hmatac.p = ggplot(pure.hmatac , aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
  scale_colour_manual(values = my_cols) +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_pure.hmatac.p = savePlot('output/batchEffect/hmatac_cluster_purity.png', pure.hmatac.p),
  # # To determine which clusters contaminate each other, 
  # we can identify the cluster with the most neighbors for each cell. 
  # In the table below, each row corresponds to one cluster; 
  # large off-diagonal counts indicate that its cells are easily confused with those from another cluster.
  #  ### rna ----
  pure.hmrna = calculate_purity(hmrna.sce, 'PCA'),
  pure.hmrna.p = ggplot(pure.hmrna, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_pure.hmrna.p = savePlot('output/batchEffect/hmrna_cluster_purity.png', pure.hmrna.p),

  # # remove potential doublets ---
  # # scdblFinder did not find doublets, manually inspect clusters. 
  # If a cluster contain a majority of cells from RMS samples and some from ATRT or MRT, these MRT and ATRT 
  # maybe potential doublets. 
  # # these cells were identified manually from clean_contaminated_cells.R
  # # potential_db =  read.csv('output/cell_type/sc_rna/clean_unknown/all_cluster_remove.csv'),
  manual_db = read.csv(file  = 'output/cell_type/sc_rna/clean_unknown/clean_cluster_potential_doublets.csv'),
  dbfinder_db = read.csv(file = 'output/cell_type/sc_rna/scdbfinder_doublet_default_setting.csv'),
  potential_db = c(manual_db[,1], dbfinder_db$x),
  potential_singlet = setdiff(colnames(hmRna_scroshi_atrt), potential_db),
  rna_wo_db =  subset(hmRna_scroshi_atrt, subset = m_barcode %in% potential_singlet),
  rna_wo_db_p = DimPlot(rna_wo_db, group.by = 'Subtype', cols = my_cols),
  save_rna_wo_db_p = savePlot(paste0(cleanUnknownRnaDir,'/subtype_after_cleaning.png'), rna_wo_db_p),
  rna_wo_db_ident = change_indent(rna_wo_db, 'RNA_snn_res.0.8'),
  # # add general subtype metadata ---
  healthy_clusters = c( 8, 18, 20, 21, 13, 29, 28, 26 ),
  rna_wo_db_general_sub = generalize_subtype(rna_wo_db_ident, 
                                             healthy_clusters, cluster_col = 'RNA_snn_res.0.8'),
  # # add infercnv result ---
  #  # infer_res was generated manually from get_infercnv_result.R
  rna_nocf_infer_res = read.csv('output/cell_type/sc_rna/infercnv/rna_remove_cf_infer_res.csv'),
  rna_nodb_infer = assign_infer_res_to_sr(rna_nocf_infer_res, rna_wo_db_general_sub),
  # # check cluster tree ---
  rna_wo_db_tree = change_tree_label(rna_wo_db_general_sub, by = 'general_subtype',
                                     save_name = 'output/cell_type/sc_rna/clean_unknown/tree_after_cleaning_general_subtype.png',assay.name = 'RNA', dims= 1:30, reduction.method = 'HARMONY', cluster.col = 'RNA_snn_res.0.8'),

  # # check cluster sil and purity after removign doublet ---
  # ### rna -----
  hmrna_wodb.sce = make_sce(rna_wo_db_general_sub),
  sil.hmrnawodb = calculate_silhouette(hmrna_wodb.sce, reduce_method = 'HARMONY'),
  sil.hmrnawodb.p = ggplot(sil.hmrnawodb, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_sil.hmrnawodb.p = savePlot('output/batchEffect/hmrna_wodb_silhouette_cluster_behavior.png', sil.hmrnawodb.p),
  # ## purity ---
  pure.hmrna_wodb = calculate_purity(hmrna_wodb.sce, 'HARMONY'),
  pure.hmrna_wodb.p = ggplot(pure.hmrna_wodb, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    theme( panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size =20),
          axis.title.y = element_text(size = 20)),
  save_pure.hmrna_wodb.p = savePlot('output/batchEffect/hmrna_wodb_cluster_purity.png', pure.hmrna_wodb.p),

  # # scroshi after removing db--
  hmRna_wodb_scroshi_demo = run_scROSHI_w_demo_data(sr = rna_nodb_infer, cols = my_cols, 
                                                    pt = 1, save_name = 'hm_rna_wodb_w_demo_marker', 
                                                    save_path = CellRnaScroshiDir),
  
  hmRna_wodb_scroshi_atrt = run_scROSHI_w_cancer_marker(sr = rna_nodb_infer, cols = my_cols, 
                                                        pt = 1, save_name = 'hm_rna_wodb_w_atrt', 
                                                        save_path = CellRnaScroshiDir)
)

assign_tumor_cell_plan <- drake_plan(
  # tumor cells were identified manually by 
  # inspecting singleR, infercnv, scroshi and markers
  # check identify_tumor_cells.R for more details 
  rna_w_tumor_label = assign_tumor_cells(hmRna_wodb_scroshi_atrt),
  rna_hm_new_bc = paste0(rna_w_tumor_label$barcodes, '_', rna_w_tumor_label$library),
  rna_w_tumor_label_newbc = RenameCells(rna_w_tumor_label, new.names = rna_hm_new_bc),
  # atac --
  atachm_new_bc = paste0(hmAtac_scroshi_atrt$barcodes, '_', hmAtac_scroshi_atrt$library),
  atac_hm_newbc = RenameCells(hmAtac_scroshi_atrt, new.names = atachm_new_bc),
  atac_hm_w_tumor_label = assign_cross_labels(des_sr = atac_hm_newbc, source_sr = rna_w_tumor_label_newbc, 
  label_col = 'cell_identity'),
  atac_hm_tumor_nona = remove_na_cells_in_metadata(atac_hm_w_tumor_label, 'cell_identity'),

    # atac harmony with gene activity
  
  atachmGA_new_bc = paste0(hm_atacGA$barcodes, '_', hm_atacGA$library),
  atac_hmGA_newbc = RenameCells(hm_atacGA, new.names = atachmGA_new_bc),
  atac_hmGA_w_tumor_label = assign_cross_labels(des_sr = atac_hmGA_newbc, source_sr = rna_w_tumor_label_newbc, 
  label_col = 'cell_identity'),
  atac_hmGA_tumor_nona = remove_na_cells_in_metadata(atac_hmGA_w_tumor_label, 'cell_identity')
)

# assign tumor cells and remove doublet from
# rna without removing patient-effect w harmony ---
no_harmony_plan <- drake_plan(
  rna_nohm_nodb = subset(rna_group_sgr, subset =  m_barcode %in% potential_singlet),
  newbc_rna_nohm = paste0(rna_nohm_nodb$barcodes, '_', rna_nohm_nodb$library),
  rna_nohm_nodb_newbc = RenameCells(rna_nohm_nodb, new.names = newbc_rna_nohm), 
  rna_nohm_tumor_label = assign_cross_labels(des_sr = rna_nohm_nodb_newbc, 
                                             source_sr = rna_w_tumor_label_newbc, 
                                             label_col = 'cell_identity'),
  # atac ----
  atac_new_bc = paste0(atac_group_sgr$barcodes, '_', atac_group_sgr$library),
  atac_newbc = RenameCells(atac_group_sgr, new.names = atac_new_bc ),
  atac_nohm_tumor_label = assign_cross_labels(des_sr = atac_newbc, 
                                              source_sr = rna_w_tumor_label_newbc, 
                                              label_col = 'cell_identity'),
  atac_nohm_tumor_nona = remove_na_cells_in_metadata(atac_nohm_tumor_label, 'cell_identity'),
  atac_nohm_tumor_ga = get_gene_activity(atac_nohm_tumor_nona)


  

  )

# identify markers ----
marker_plan <- drake_plan(
  # # calculate markers ---
  # ref https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
   rna_markers = FindAllMarkers(object =  rna_w_tumor_label_newbc, only.pos = T, logfc.threshold = 0.25),
  ## find markers that distinguish clusters from the same cancer type to others
  p3w_marker = FindMarkers(rna_nodb_infer, ident.1 = 6),
   ## show it visually:
# DoHeatmap(subset(srat, downsample = 50),
#           features = top10$gene,
#           group.colors=cluster.colors,
#           assay='RNA',
          # slot='scale.data'
          # ) + theme(axis.text.y=element_text(size=6)) # try here https://github.com/scgenomics/scgenomics-public.github.io/blob/main/docs/14-enrich/14-enrich.R,
atac_hm_cell_typeIdent = change_indent(atac_hm_tumor_nona, by = 'cell_identity'),
atac_markers = FindAllMarkers(object = atac_hm_cell_typeIdent, 
                              only.pos = T, logfc.threshold = 0.25)
)


# process healthy data to compare ----
# healthy data are obtained online from descartes
healthy_plan <- drake_plan(
    hthysr = target(readRDS(paste0(hthy_dataDir, '/', ts, '_filtered.seurat.for_website.RDS')),
                transform = map(ts = !!hthytissue_list ,
                                id.vars = !!hthytissue_list ,
                                .id = id.vars)),
    # process atac hthy ---
    atac_hthysr = target(set_default_assay(hthysr, assay = 'peaks'),
        transform = map(hthysr ,
                        id.vars = !!hthytissue_list ,
                        .id = id.vars)),
    atac_hthysrFill = target(filter_outliers_healthyAtac(atac_hthysr, healthyFigDir),
                transform = map(atac_hthysr,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    atac_hthyNor = target(sc_atac_normalize(atac_hthysrFill),
                transform = map(atac_hthysrFill,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    atac_hthyDim = target(sc_atac_dim_redu(atac_hthyNor),
                transform = map(atac_hthyNor,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    atac_hthySubset = target(sampling_sr(atac_hthyDim, 
                                         percent_to_keep = 800, 
                                         type = 'number', 
                                         class_col = 'cell_type'),
                transform = map(atac_hthyDim,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    atac_hthymrg = target(merge_sr_list(c(atac_hthySubset), healthyDir),
                transform = combine(atac_hthySubset,
                id.vars = !!hthytissue_list ,
                .id = id.vars)),
    atac_hthymrgNor = sc_atac_normalize(atac_hthymrg),
    atac_hthymrgDim = sc_atac_dim_redu(atac_hthymrgNor),
    # atac_hthymrgGA = get_gene_activity_wo_fragFile(atac_hthymrgDim, hg38),

    # mrg descartes atac download ---
    desc_atac_download = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/all_celltypes.downsampled.filtered.RDS'),
    # atac downloaded online have no fragment file, so use custom-made function to get gene activity 
    # desc_atac_download_ga = get_gene_activity_wo_fragFile(desc_atac_download, hg38),

    # create the same gene activity for multiome data to compare with reference data --
    # atac_ga_nofrag = get_gene_activity_wo_fragFile(atac_group_sgr, hg38),

    # merg rna hthy ---
    rna_hthySubset = target(sampling_sr(hthysr, percent_to_keep = 800, 
                                        type = 'number', class_col = 'cell_type'),
                transform = map(hthysr,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    rna_hthymrg = target(merge_sr_list(c(rna_hthySubset), healthyDir),
                transform = combine(rna_hthySubset,
                id.vars = !!hthytissue_list ,
                .id = id.vars)),
    rna_hthymrg_nor = normalize_dim_plot_sr(rna_hthymrg, save_path = healthyDir, 
                                            lib_name = 'hthy_all_rna' ),
    rna_hthymrg_clus = clustering_rna_data(rna_hthymrg_nor)
    
)

logistic_rna_plan <- drake_plan(
  # with descardes data ---
  train_rna_40k = trainModel(GetAssayData(rna_hthymrg_clus),
                            classes = rna_hthymrg_clus$cell_type, 
                            maxCells = 40000),
  predict_rna_40k = predictSimilarity(train_rna_40k,
        GetAssayData(rna_w_tumor_label_newbc),
        classes = rna_w_tumor_label_newbc$cell_identity,
        minGeneMatch = 0.7, logits = FALSE),
   train_rna_allCells = trainModel(GetAssayData(rna_hthymrg_clus),
                            classes = rna_hthymrg_clus$cell_type, 
                            maxCells = 82300),
  predict_rna_allCells = predictSimilarity(train_rna_allCells,
        GetAssayData(rna_w_tumor_label_newbc),
        classes = rna_w_tumor_label_newbc$cell_identity,
        minGeneMatch = 0.7, logits = FALSE),
  # xi 2020 ---
  xi_2020 = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Jeff_rf/xi_2020.rds'),
  xi_2020_nor = normalize_dim_plot_sr(xi_2020, save_path = healthyDir, lib_name = 'xi_2020' ),
  xi_2020_dim = clustering_rna_data(xi_2020_nor),

  xi_train = trainModel(GetAssayData(xi_2020_dim),
    classes = xi_2020_dim$cell_type,
    maxCells = 40000),
  p_xi = predictSimilarity(xi_train,
    GetAssayData(rna_w_tumor_label_newbc),
    classes= rna_w_tumor_label_newbc$cell_identity,
    minGeneMatch=0.70, logits = FALSE)

)

logistic_atac_plan <- drake_plan(
  # # convert peak coordinate to be comparable ---
  dsc_atacchr = createSrWChromatinAssay(atac_hthymrgDim, hg38),
  atacchr = createSrWChromatinAssay(atac_hm_tumor_nona, hg38),
  atac_gr = granges(atacchr),
  dsc_gr = granges(dsc_atacchr),
  atac_hm_count = GetAssayData(atac_hm_tumor_nona),
  # reduce data size to speed up the make matrix func
  # remove features with most count 0 
  freq_0 = rowSums(atac_hm_count !=0),
  # keep only features with at least 100 cells non 0 
  atac_non0 = atac_hm_count[freq_0 > 300,], # reduce to 972263 features from 1760372 in atac_hm 
  feature_tokeep = rownames(atac_hm_count)[freq_0 >300],
  atac_gr_df = as.data.frame(atac_gr),
  index_tokeep = paste0(atac_gr_df$seqnames, '-', atac_gr_df$start, '-', atac_gr_df$end) %in% feature_tokeep,
  new_atac_gr = atac_gr[index_tokeep],
  # new_atachm_mx = makeCountMx_withSamePeaks_optimized3(dsc_gr,new_atac_gr, atac_non0),
  new_atachm_mx = readRDS('output/logistic_regression/atac_hm_features_above_300cells.RDS'),
  new_atachmMx_colname = assign_colname_newMx(new_atachm_mx, atac_hm_tumor_nona),
  # train atac dsc ---
  # train_dsc_atac40k = trainModel(GetAssayData(atac_hthymrgDim), classes = atac_hthymrgDim$cell_type, maxCells = 40000),
  # # train only overlap features ---
  # sub_dsc_atac = subset(atac_hthymrgDim, features = rownames(new_atachm_mx)),
  # tran_sub_dsc_atac = trainModel(GetAssayData(sub_dsc_atac), classes = atac_hthymrgDim$cell_type, maxCells = 40000),

  # p_sub_dsc_atac = predictSimilarity(tran_sub_dsc_atac, new_atachmMx_colname, classes = atac_hm_tumor_nona$cell_identity,
  # minGeneMatch = 0.2, logits = FALSE ), # out of memory 350Gb

  # add bc metadata ---
  dsc_atac = add_barcode_metadata(atac_hthyDim),
  dsc_atac_ident = change_indent(dsc_atac, by = 'cell_type'),
  dsc_markers = FindAllMarkers(dsc_atac_ident, only.pos = T, logfc.threshold = 0.25),
  atac_features = rownames(new_atachm_mx),

  topfeatures7k = dsc_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 7000, 
          wt = avg_log2FC),
  
  features_to_keep7k = topfeatures7k$gene,

  train_feature7k = intersect(features_to_keep7k, atac_features),
  sub_dsc7k = subset(dsc_atac_ident, features = train_feature7k),
  sub_dsc7k75 = sampling_sr(sub_dsc7k, 75, class_col = 'cell_type', type = 'percent'),
  
  sub_dsc_25bc = setdiff(colnames(sub_dsc7k), colnames(sub_dsc7k75)),
  sub_dsc7k25 =  subset(sub_dsc7k, subset = cell_bc %in% sub_dsc_25bc),
  
  # train ----
  train_dsc7k = trainModel(GetAssayData(sub_dsc7k75), class = sub_dsc7k75$cell_type, maxCell = 82300),
  # test ----
  p_dsc7k_test25 = predictSimilarity(train_dsc7k, GetAssayData(sub_dsc7k25), 
                                     classes = sub_dsc7k25$cell_identity, 
                                     logits = F, minGeneMatch = 0.0),
  # predict ----
  sub_atac7k = new_atachm_mx[rownames(new_atachm_mx) %in% train_feature7k,], # 9760 features
  p_dsc7k = predictSimilarity(train_dsc7k, sub_atac7k, classes = atac_hmIdent$cell_identity, 
                              logits = F, minGeneMatch = 0.0)
  # 10k features ---
  topfeatures10k = dsc_markers %>% 
    group_by(cluster) %>% 
    top_n(n = 10000, 
          wt = avg_log2FC),
  
  features_to_keep10k = topfeatures10k$gene,
  train_feature10k = intersect(features_to_keep10k, atac_features),
  sub_dsc10k = subset(dsc_atac_ident, features = train_feature10k),
  sub_dsc10k75 = sampling_sr(sub_dsc10k, 75, class_col = 'cell_type', type = 'percent'),
  
  sub_dsc_25bc = setdiff(colnames(sub_dsc10k), colnames(sub_dsc10k75)),
  sub_dsc10k25 =  subset(sub_dsc10k, subset = cell_bc %in% sub_dsc_25bc),
  
  # train ----
  train_dsc10k = trainModel(GetAssayData(sub_dsc10k75), class = sub_dsc10k75$cell_type, maxCell = 2000),
  # test ----
  p_dsc10k_test25 = predictSimilarity(train_dsc10k, GetAssayData(sub_dsc10k25), 
                                      classes = sub_dsc10k25$cell_identity, 
                               logits = F, minGeneMatch = 0.0),
  
  # predict ----
  sub_atac10k = new_atachm_mx[rownames(new_atachm_mx) %in% train_feature10k,], # 9760 features
  p_dsc10k = predictSimilarity(train_dsc10k, sub_atac10k, classes = atac_hmIdent$cell_identity, 
                             logits = F, minGeneMatch = 0.0)
  
)



plan <- bind_plans(combine_peak_plan, process_special_lib_plan, 
                  process_plan, cell_annotation_plan, 
                  cluster_behavior_plan, batch_detection_plan, 
                  batch_correction_plan, 
                  cluster_behavior_after_correction_plan, 
                  marker_plan, assign_tumor_cell_plan,
                  no_harmony_plan,healthy_plan, logistic_rna_plan,
                  logistic_atac_plan)

drake_config(plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
