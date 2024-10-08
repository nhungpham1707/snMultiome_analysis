# ---------------------------------------
# This script is used to defined and 
# create variables that will be called 
# in the analysis.
# These include: directory where data 
# and results are stored. Variables 
# such as filter thredshold, and
# figures configuration
# 
# Nhung 23 11 2023
# -------------------------------------------
#  Directory explain: 
# base_sc_dir: a base location of the project 
# where data and output locate
# base_data_dir: where data is (h5 count matrix, fragment files, souporcell) 
# output directories are organized base
# on the type of outputs: 
#   - report_dir: any statistic csv files, 
#   e.g. cells numbers before 
#   and after processing
#   - rnaProcessDir/
#    atacProcessDir: process files 
#    for scATAC and scRNA 
#   - atacDemultiFigDir/
#     atacProcessFigDir figures, 
#     e.g. dimplot, vlnplot of data
#   - atcMrgDir : 
#     merge seurats of all samples.
#     this file is the amin file to use 
#     for downstream analysis 
# ---------------------------------------

# define dir -------
base_sc_dir <- "/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu" # modify if change project
base_data_dir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/analyses/analyses'  # modify if change project
output_dir <- paste0(base_sc_dir, '/output')
report_dir <- paste0(output_dir,'/report')
dataset_overview_fig_dir <- paste0(output_dir,'/figures/dataset')
cell_type_dir <- paste0(output_dir,'/cell_type')
cell_type_fig_dir <- paste0(output_dir,'/figures/cell_type')

rnaProcessDir <-paste0(output_dir,'/sc_RNA/processing') 
rnaFigDir <- paste0(output_dir, '/figures/sc_RNA/processing')
rnaMrgDir <- paste0(output_dir,'/sc_RNA/merge_all')
rnaMrgFigDir <- paste0(output_dir,'/figures/sc_RNA/merge_all')
sc_rna_demultiplex_fig_dir <- paste0(output_dir,'/figures/sc_RNA/demultiplex')
cell_type_rna_dir <- paste0(cell_type_dir, '/sc_rna')
cell_type_rna_fig_dir <- paste0(cell_type_fig_dir, '/sc_rna')

cellRNAsingRdir <- paste0(cell_type_rna_dir, '/singler')
cell_type_rna_singler_fig_dir <- paste0(cell_type_rna_fig_dir, '/singler')

cellRnaIcnvdir <- paste0(cell_type_rna_dir, '/infercnv')
rnaInferInputDir <- paste0(cellRnaIcnvdir, '/input')
cell_type_rna_infercnv_fig_dir <- paste0(cell_type_rna_fig_dir, '/infercnv')
cell_type_rna_scroshi_dir <- paste0(cell_type_rna_dir, '/scROSHI')

atacProcessDir<-paste0(output_dir,'/sc_atac/processing') 
atacProcessFigDir <- paste0(output_dir,'/figures/sc_atac/processing')
atacnoRProcessFigDir <- paste0(output_dir,'/figures/sc_atac/processing_noRelapse')
healthyFigDir <- paste0(output_dir, '/figures/healthy_data')
healthyDir <- paste0(output_dir, '/healthy_data')
atacDemultiFigDir <- paste0(output_dir,'/figures/sc_atac/demultiplex')
atacMrgFigDir <- paste0(output_dir,'/figures/sc_atac/merge_all')
atacnoRMrgFigDir <- paste0(output_dir,'/figures/sc_atac/merge_all_no_relapse')
atcMrgDir <- paste0(output_dir,'/sc_atac/merge_all')
atcnoRMrgDir <- paste0(output_dir,'/sc_atac/merge_all_no_relapse')
cell_type_atac_dir <- paste0(cell_type_dir, '/sc_atac')
cell_type_atac_fig_dir <- paste0(cell_type_fig_dir, '/sc_atac')
cellAtacInferDir <- paste0(cell_type_atac_dir, '/infercnv')
cellAtacnoRInferDir <- paste0(cell_type_atac_dir, '/infercnv/no_relapse')
AtacInferInputDir <- paste0(cellAtacInferDir, '/Input')
AtacnoRInferInputDir <- paste0(cellAtacInferDir, '/InputnoRelapse')
cell_type_atac_infercnv_fig_dir <- paste0(cell_type_atac_fig_dir, '/infercnv')
cellAtacInferMergDir <- paste0(cellAtacInferDir, '/merge_sr')
atacCellSngRDir <- paste0(cell_type_atac_dir, '/singler')
atacCellSngRFigDir <- paste0(cell_type_atac_fig_dir, '/singler')
atacnoRCellSngRFigDir <- paste0(cell_type_atac_fig_dir, '/singlerNoRelap')

atacScroshiFigDir <- paste0(cell_type_atac_fig_dir, '/scROSHI')
atacScroshiDir <- paste0(cell_type_atac_dir, '/scROSHI')
hyperparameter_atac_dir <- paste0(atcMrgDir, '/hyperparameter')
clustering_atac_dir <- paste0(atcMrgDir, '/cluster')
batchDir <- paste0(output_dir, '/batchEffect')
batchAtacDir <- paste0(batchDir, '/atac')
batchRnaDir <- paste0(batchDir, '/rna')
# create dir -----
# figures
dir.create("output/figures", recursive = TRUE)
dir.create(report_dir, recursive = TRUE)
dir.create(dataset_overview_fig_dir , recursive = TRUE)
dir.create(cell_type_dir, recursive = TRUE)
dir.create(cell_type_fig_dir, recursive = TRUE)
# output for scRNA
dir.create(rnaProcessDir, recursive = TRUE)
dir.create(rnaFigDir, recursive = TRUE)
dir.create(sc_rna_demultiplex_fig_dir, recursive = TRUE)
dir.create(cell_type_rna_dir, recursive = TRUE)
dir.create(cell_type_rna_fig_dir, recursive = TRUE)
dir.create(cellRnaIcnvdir, recursive = TRUE)
dir.create(cell_type_rna_infercnv_fig_dir, recursive = TRUE)
dir.create(cellRNAsingRdir, recursive = TRUE)
dir.create(cell_type_rna_singler_fig_dir, recursive = TRUE)
dir.create(cell_type_rna_scroshi_dir, recursive = TRUE)
# output for scATAC
dir.create(atacProcessDir, recursive = TRUE)
dir.create(atacProcessFigDir, recursive = TRUE)
dir.create(atacDemultiFigDir, recursive = TRUE)
dir.create(atcMrgDir, recursive = TRUE)
dir.create(atacMrgFigDir, recursive = TRUE)

dir.create(cell_type_atac_dir, recursive = TRUE)
dir.create(cell_type_atac_fig_dir, recursive = TRUE)
dir.create(cellAtacInferDir, recursive = TRUE)
dir.create(cell_type_atac_infercnv_fig_dir, recursive = TRUE)
dir.create(cellAtacInferMergDir, recursive = TRUE)
dir.create(atacCellSngRDir, recursive = TRUE)
dir.create(atacCellSngRFigDir, recursive = TRUE)
dir.create(atacScroshiFigDir, recursive = TRUE)
dir.create(hyperparameter_atac_dir, recursive = TRUE)
dir.create(clustering_atac_dir, recursive = TRUE)
dir.create(AtacInferInputDir, recursive = TRUE)
dir.create(healthyFigDir, recursive = TRUE)
dir.create(healthyDir, recursive = TRUE)
dir.create(atacnoRProcessFigDir, recursive = TRUE)
dir.create(atacnoRCellSngRFigDir, recursive = TRUE)
dir.create(AtacnoRInferInputDir, recursive = TRUE)
dir.create(atcnoRMrgDir, recursive = TRUE)
dir.create(atacnoRMrgFigDir, recursive = TRUE)
dir.create(cellAtacnoRInferDir, recursive = TRUE)
dir.create(atacScroshiDir, recursive = TRUE)
dir.create(rnaMrgFigDir, recursive = TRUE)
dir.create(rnaMrgDir, recursive = TRUE)
dir.create(batchDir, recursive = TRUE)
dir.create(batchAtacDir, recursive = TRUE)
dir.create(batchRnaDir, recursive = TRUE)
dir.create(rnaInferInputDir, recursive = TRUE)
