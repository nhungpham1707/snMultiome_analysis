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
#   - sc_rna_processing_dir/
#    atacProcessDir: process files 
#    for scATAC and scRNA 
#   - atacDemultiFigDir/
#     atacProcessFigDir figures, 
#     e.g. dimplot, vlnplot of data
#   - sc_atac_merge_sr_dir : 
#     merge seurats of all samples.
#     this file is the amin file to use 
#     for downstream analysis 
# ---------------------------------------

# define dir -------
base_sc_dir <- "/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code" # modify if change project
base_data_dir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/analyses/analyses'  # modify if change project
output_dir <- paste0(base_sc_dir, '/output')
report_dir <- paste0(output_dir,'/report')
dataset_overview_fig_dir <- paste0(output_dir,'/figures/dataset')
cell_type_dir <- paste0(output_dir,'/cell_type')
cell_type_fig_dir <- paste0(output_dir,'/figures/cell_type')

sc_rna_processing_dir <-paste0(output_dir,'/sc_RNA/processing') 
sc_rna_processing_fig_dir <- paste0(output_dir, '/figures/sc_RNA/processing')
sc_rna_demultiplex_fig_dir <- paste0(output_dir,'/figures/sc_RNA/demultiplex')
cell_type_rna_dir <- paste0(cell_type_dir, '/sc_rna')
cell_type_rna_fig_dir <- paste0(cell_type_fig_dir, '/sc_rna')

cell_type_rna_singler_dir <- paste0(cell_type_rna_dir, '/singler')
cell_type_rna_singler_fig_dir <- paste0(cell_type_rna_fig_dir, '/singler')

cell_type_rna_infercnv_dir <- paste0(cell_type_rna_dir, '/infercnv')
cell_type_rna_infercnv_fig_dir <- paste0(cell_type_rna_fig_dir, '/infercnv')
cell_type_rna_scroshi_dir <- paste0(cell_type_rna_dir, '/scROSHI')

atacProcessDir<-paste0(output_dir,'/sc_atac/processing') 
atacProcessFigDir <- paste0(output_dir,'/figures/sc_atac/processing')
atacDemultiFigDir <- paste0(output_dir,'/figures/sc_atac/demultiplex')
sc_atac_merge_sr_fig_dir <- paste0(output_dir,'/figures/sc_atac/merge_all')
sc_atac_merge_sr_dir <- paste0(output_dir,'/sc_atac/merge_all')

cell_type_atac_dir <- paste0(cell_type_dir, '/sc_atac')
cell_type_atac_fig_dir <- paste0(cell_type_fig_dir, '/sc_atac')
cellAtacInferDir <- paste0(cell_type_atac_dir, '/infercnv')
cell_type_atac_infercnv_fig_dir <- paste0(cell_type_atac_fig_dir, '/infercnv')
cellAtacInferMergDir <- paste0(cellAtacInferDir, '/merge_sr')
MergAtacInferInputDir <- paste0(cellAtacInferMergDir, '/Input')

atacCellSngRDir <- paste0(cell_type_atac_dir, '/singler')
cell_type_atac_singler_fig_dir <- paste0(cell_type_atac_fig_dir, '/singler')

cell_type_atac_scroshi_dir <- paste0(cell_type_atac_dir, '/scROSHI')

hyperparameter_atac_dir <- paste0(sc_atac_merge_sr_dir, '/hyperparameter')
clustering_atac_dir <- paste0(sc_atac_merge_sr_dir, '/cluster')
# create dir -----
# figures
dir.create("output/figures", recursive = TRUE)
dir.create(report_dir, recursive = TRUE)
dir.create(dataset_overview_fig_dir , recursive = TRUE)
dir.create(cell_type_dir, recursive = TRUE)
dir.create(cell_type_fig_dir, recursive = TRUE)
# output for scRNA
dir.create(sc_rna_processing_dir, recursive = TRUE)
dir.create(sc_rna_processing_fig_dir, recursive = TRUE)
dir.create(sc_rna_demultiplex_fig_dir, recursive = TRUE)
dir.create(cell_type_rna_dir, recursive = TRUE)
dir.create(cell_type_rna_fig_dir, recursive = TRUE)
dir.create(cell_type_rna_infercnv_dir, recursive = TRUE)
dir.create(cell_type_rna_infercnv_fig_dir, recursive = TRUE)
dir.create(cell_type_rna_singler_dir, recursive = TRUE)
dir.create(cell_type_rna_singler_fig_dir, recursive = TRUE)
dir.create(cell_type_rna_scroshi_dir, recursive = TRUE)
# output for scATAC
dir.create(atacProcessDir, recursive = TRUE)
dir.create(atacProcessFigDir, recursive = TRUE)
dir.create(atacDemultiFigDir, recursive = TRUE)
dir.create(sc_atac_merge_sr_dir, recursive = TRUE)
dir.create(sc_atac_merge_sr_fig_dir, recursive = TRUE)

dir.create(cell_type_atac_dir, recursive = TRUE)
dir.create(cell_type_atac_fig_dir, recursive = TRUE)
dir.create(cellAtacInferDir, recursive = TRUE)
dir.create(cell_type_atac_infercnv_fig_dir, recursive = TRUE)
dir.create(cellAtacInferMergDir, recursive = TRUE)
dir.create(atacCellSngRDir, recursive = TRUE)
dir.create(cell_type_atac_singler_fig_dir, recursive = TRUE)
dir.create(cell_type_atac_scroshi_dir, recursive = TRUE)
dir.create(hyperparameter_atac_dir, recursive = TRUE)
dir.create(clustering_atac_dir, recursive = TRUE)
dir.create(MergAtacInferInputDir, recursive = TRUE)
# define variable that will be used in all analysis 
reso <- 300 # figures resolution

# thredshold for scRNA 
sc_rna_dims <- c(1:15) # default pca. 
#check elbow plot of each library for more specific one. 
#those before the plot flatten

# thredshold for outlier detection in isoutlier 
# and in singler prunescore filter
nmad_threshold <- 3
# male and stress genes list to remove. 
# from SCutils package 
data(refdata_gex_GRCh38_2020_A) 
male.genes <- genelists$male
stress_gene_list <- genelists$stress

# normal cells from singler human dex for infercnv
normal_cells <- c('Pre-B_cell_CD34-', 'Pro-B_cell_CD34+', 
                  'B_cell', 'T_cells', 'Macrophage', 
                  'Monocyte', 'Neutrophils', 'NK_cells') 
# colors for plot
my_cols <- c("pink1", '#ccb1f1', "violet", "slateblue1", "purple3",
            "turquoise2", "skyblue", "steelblue",'#25aff5', "blue2", 
            "navyblue", "orange", "tomato", "coral2", "palevioletred", 
            "violetred", "red2", "springgreen2", "yellowgreen", 
            "palegreen4", '#1fa195', "wheat2", "#d4d915", "yellow4", 
            "yellow", "wheat4", "tan", "tan3", "brown",
            "grey70", "grey30", 'maroon4', 'maroon2', 'maroon', 
            'lightsalmon', 'lightsalmon3', 'lightsalmon4', 
            'pink4', 'peru', 'mediumpurple4', 'mediumpurple', 'lightblue4',
            'aquamarine4', 'aquamarine', 
            'cornsilk4', 'black', 'bisque4', 'darkred')
            
my_cols2 <- c('#25aff5', "blue2", "navyblue", "turquoise2", "skyblue",
              "steelblue", "wheat4", "tan", "tan3", "brown",
             "pink1", '#ccb1f1', "violet", "slateblue1", "purple3",
             "orange", "tomato", "coral2", "palevioletred", "violetred", 
             "red2", "springgreen2", "yellowgreen", "palegreen4", 
             '#1fa195', "wheat2", "#d4d915", "yellow4", "yellow", 
             "grey70", "grey30", 'maroon4', 'maroon2', 'maroon', 
             'lightsalmon', 'lightsalmon3',
             'lightsalmon4', 'pink4', 'peru', 'mediumpurple4', 
             'mediumpurple', 'lightblue4', 'aquamarine4', 'aquamarine' )