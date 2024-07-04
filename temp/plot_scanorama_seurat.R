source('R/utils.R')
library(Seurat)
library(SeuratDisk)
library(SCpubr)
data = 'rna'
save_dir = '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/batchEffect/', data, '/scanorama/'

# Convert(paste0(save_dir,"/atac_scanorama.h5ad"), dest = "h5seurat", overwrite = TRUE, assay = 'peaks')
Convert(paste0(save_dir,"/", data, "_scanorama.h5ad"), dest = "h5seurat", overwrite = TRUE)

sr <- LoadH5Seurat(paste0(save_dir,"/", data, "_scanorama.h5seurat"))

# colors for scpurb
colors <- c("ATRT_SHH" = "skyblue" ,                     
            "ecMRT"  = "blue2",                       
            "ecMRT_BrainMet"  = "navyblue",                 
            "ATRT_TYR"     = "turquoise2",                   
            "ATRT_MYC"  = '#25aff5',                      
            "FN-eRMS" = "tan",                        
            "SySa"    = "pink1",                         
            "ATRT-MYC (replapse from ATRT21)"='maroon4', 
            "FP-RMS (P3W)"  = "tan3",                  
            "FP-RMS (P3F)" = "brown",                    
            "FP-RMS" = "slateblue1")

scpur_p_lib = do_DimPlot(sr , group.by = 'library')
scpur_p_type = do_DimPlot(sr , group.by = 'Subtype', colors.use = colors)
  scpur_p = scpur_p_lib | scpur_p_type, 
savePlot(paste0(save_dir,'/umap_type_lib.png'), scpur_p)
savePlot(paste0(save_dir,'/umap_lib.png'), scpur_p_lib)
