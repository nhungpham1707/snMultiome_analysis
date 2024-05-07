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
base_sc_dir <- getwd()
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
cellRnaDir <- paste0(cell_type_dir, '/sc_rna')
cellRnaFigDir <- paste0(cell_type_fig_dir, '/sc_rna')
cellRnaMarkerDir <- paste0(cellRnaDir, '/markers')
cellRNAsingRdir <- paste0(cellRnaDir, '/singler')
rnaCellSngRFigDir <- paste0(cellRnaFigDir, '/singler')

cellRnaIcnvdir <- paste0(cellRnaDir, '/infercnv')
rnaInferInputDir <- paste0(cellRnaIcnvdir, '/input')
cell_type_rna_infercnv_fig_dir <- paste0(cellRnaFigDir, '/infercnv')
CellRnaScroshiDir<- paste0(cellRnaDir, '/scROSHI')

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
batchAtacHarmonyDir <- paste0(batchAtacDir, '/harmony')
atacSysviDir <- paste0(batchAtacDir, '/sysvi')
batchAtacScanormaDir <- paste0(batchAtacDir, '/scanorama')

batchRnaDir <- paste0(batchDir, '/rna')
batchRnaHarmonyDir <- paste0(batchRnaDir, '/harmony')
batchRnaSysviDir <- paste0(batchRnaDir, '/sysvi')
batchRnaScanoramaDir <- paste0(batchRnaDir, '/scanorama')
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
dir.create(cellRnaDir, recursive = TRUE)
dir.create(cellRnaFigDir, recursive = TRUE)
dir.create(cellRnaIcnvdir, recursive = TRUE)
dir.create(cell_type_rna_infercnv_fig_dir, recursive = TRUE)
dir.create(cellRNAsingRdir, recursive = TRUE)
dir.create(rnaCellSngRFigDir, recursive = TRUE)
dir.create(CellRnaScroshiDir, recursive = TRUE)
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
dir.create(CellRnaScroshiDir, recursive = TRUE)
dir.create(atacSysviDir, recursive = TRUE)
dir.create(batchAtacHarmonyDir, recursive = TRUE)
dir.create(batchRnaHarmonyDir, recursive = TRUE)
dir.create(batchRnaSysviDir, recursive = TRUE)
dir.create(batchAtacScanormaDir, recursive = TRUE)
dir.create(batchRnaScanoramaDir , recursive = TRUE)
dir.create(cellRnaMarkerDir, recursive = TRUE)
# define variable that will be used in all analysis 
reso <- 300 # figures resolution

# thredshold for scRNA 
sc_rna_dims <- c(1:15) # default pca. 
#check elbow plot of each library for more specific one. 
#those before the plot flatten
mt_percent_limit <- 5
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

# Define cell groups: (from primary human atlas)
# ref https://phys.org/news/2021-07-cell-human-tissues.html
immune_cells <- c("Monocyte", "T_cells", "Neutrophils", "Macrophage", "B_cell", "NK_cell", "DC", "Pro-B_cell_CD34+", "BM & Prog.", "CMP", "Erythroblast", "GMP", "Myelocyte" , "Pro-Myelocyte", "MEP", "HSC_CD34+", "HSC_-G-CSF", "Platelets", "Pre-B_cell_CD34-", "BM")
brain_cells <- c("Astrocyte", "Neurons", "Neuroepithelial_cell")
bone_cells <- c("Chondrocytes", "Osteoblasts")
stem_cells <- c("Embryonic_stem_cells", "iPS_cells", "Tissue_stem_cells", "MSC", "Gametocytes")
stroma_cells <- c("Fibroblasts", "Endothelial_cells", "Keratinocytes", "Smooth_muscle_cells" )
# normal cells from literature https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9886402/

nonMalignant_cells = c('B_cell', 'BM & Prog.', 'DC', 'CMP', 'Monocyte', 'Neutrophils', 'NK_cell', 'Pro-B_cell_CD34+', 'T_cells', 'Endothelial_cells', 'Erythroblast', 'Gametocytes', 'GMP', 'HSC_-G-CSF', 'HSC_CD34+', 'Macrophage', 'MEP', 'MSC', 'Myelocyte', 'Osteoblasts', 'Platelets', 'Pre-B_cell_CD34−', 'Pro-Myelocyte')
# immune_cells = c('B_cell', 'BM & Prog.', 'DC', 'CMP', 'Monocyte', 'Neutrophils', 'NK_cell', 'Pro-B_cell_CD34+', 'T_cells', 'Gametocytes', 'GMP', 'HSC_-G-CSF', 'HSC_CD34+', 'Macrophage', 'MEP', 'MSC', 'Platelets', 'Pre-B_cell_CD34−', 'Pro-Myelocyte')

## markers ----

ATRT_TYR <- c('MITF', 'OTX2', 'TYR', 'PDGFRB', 'JAK1', 'BMP4')
ATRT_SHH <- c('NOTCH1', 'GLI2', 'MYCN', 'ASCL1', 'HES1', 'DTX1', 'PTCH1', 'BOC')
ATRT_MYC <- c('HOXC10', 'CCND3', 'MYC')
RMS <- c('DESMIN', 'MYOG', 'MYOD1',"CD24", "CD44", "CD133", "ALDH1A1")
rms_csc <- c("CD24", "CD44", "CD133", "ALDH1A1" ) # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9406733/
rms_uri <- c("EPS8L2", "SPARC", "HLA-DRB1", "ACAN", "CILP") # https://clinicalproteomicsjournal.biomedcentral.com/articles/10.1186/s12014-023-09401-4
Sysa <- c('TLE1', 'BMP5', 'BMP7', 'TNFRS19')
t_nk <- c("CD3D", "CD8A", "GNLY")
myeloid <- c("APOE", "CD14", "C1QB")
b_cell <- c("PAX5", " MS4A1", "CD19")
endo <- c("VWF", "PECAM1", "CLDN5")
sysa_fusion_targets <- c('SS18', 'SSX1', 'SSX2')
rms_fusion_targets <- c('PAX3')

# fibroblast is one of the most abundant stroma cells
# it also in the microenvironment of the tumor. 
# fibroblast associate w tumor can be 
# identified with these markers https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5983719/
fibroblast <- c('PDGFRA','PDGFRB', "SMA", 'FAP', 'FSP1', 'S100A4')
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7719961/
fibroblast_2 <- c('ACTA2', 'PDGFRA', 'PDGFRB', 'ASPN', 'DCN')
hepatocyte_malig <- c('ALB', 'EPCAM', 'KRT19', 'ALDH1A1', 'CD24')
endothelial <- c('PECAM1', 'VWF', 'ENG', 'CDH5', 'RAMP3')
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
            
my_cols2 <- c('#25aff5',  "navyblue", 
               "tan", "tan3", "brown",
              "yellow", "slateblue1","violet", "purple3",
             "orange", "tomato", "palevioletred", "violetred", "skyblue",
             "red2", "blue2","springgreen2", '#ccb1f1',"yellowgreen", "palegreen4", 
             '#1fa195', "wheat2",  "turquoise2","#d4d915", "coral2", "yellow4",  
             "grey70", "grey30", 'maroon4', 'maroon2', "pink1",  "wheat4", "steelblue",'maroon', 
             'lightsalmon', 'lightsalmon3',
             'lightsalmon4', 'pink4', 'peru', 'mediumpurple4', 
             'mediumpurple', 'lightblue4', 'aquamarine4', 'aquamarine' )
my_cols3 <- c("grey",'#25aff5',  "navyblue", 
               "tan", "tan3", "brown",
              "yellow", "slateblue1","violet", "purple3",
             "orange", "tomato", "palevioletred", "violetred", "skyblue",
             "red2", "blue2","springgreen2", '#ccb1f1',"yellowgreen", "palegreen4", 
             '#1fa195', "wheat2",  "turquoise2","#d4d915", "coral2", "yellow4",  
             "grey70",  'maroon4', 'maroon2', "pink1",  "wheat4", "steelblue",'maroon', 
             'lightsalmon', 'lightsalmon3',
             'lightsalmon4', 'pink4', 'peru', 'mediumpurple4', 
             'mediumpurple', 'lightblue4', 'aquamarine4', 'aquamarine' )
# colors for scpurb
scpu_colors <- c("ATRT_SHH" = "skyblue" ,                     
            "ecMRT"  = "blue2",                       
            "ecMRT_BrainMet"  = "#805100",                 
            "ATRT_TYR"     = "turquoise2",                   
            "ATRT_MYC"  = '#a925f5',                      
            "FN-eRMS" = "tan",                        
            "SySa"    = "pink1",                         
            "ATRT-MYC (replapse from ATRT21)"='maroon4', 
            "FP-RMS (P3W)"  = "tan3",                  
            "FP-RMS (P3F)" = "brown",                    
            "FP-RMS" = "slateblue1")

scpu_colors2 = c( "ATRT_MYC" = "pink1",
            "ATRT_SHH" = '#ccb1f1', 
            "ATRT_TYR" = "violet", 
            "ATRT-MYC (replapse from ATRT21)" = "slateblue1", 
            "ecMRT" = "purple3",
            "ecMRT_BrainMet" = "turquoise2", "FN-eRMS" = "skyblue", 
            'FP-RMS' = "steelblue",
            "FP-RMS (P3F)" = '#25aff5', 
            "FP-RMS (P3W)" = "blue2", 
            "SySa" = "navyblue")