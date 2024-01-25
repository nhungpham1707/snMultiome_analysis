setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
lapply(list.files("./R", full.names = TRUE), source)

## read metadata ----
filename <- '30_11_2023_all_multiome_libraries.csv' 
metadata <- getData(filename, delim = ',')
orimet <- metadata
metadata <- metadata[1:3,]
lbLst <- unique(metadata$name)
idLst <- lbLst %>% map(splitName)

save_path <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/output_19_12_2023/cell_type/sc_atac/infercnv/input'

make_anno_count_matrix <- function(lib){
  singleR_res = readRDS(paste0('/hpc/pmc_drost/PROJECTS/cell_origin_NP/output_19_12_2023/cell_type/sc_atac/singler/singleR_sc_atac_', lib, '.RDS'))
  labels <- data.frame(singleR_labels = singleR_res$labels,
                     singleR_pruned.labels = singleR_res$pruned.labels)

  rownames(labels) <- rownames(singleR_res)
  sr <- readRDS(paste0('/hpc/pmc_drost/PROJECTS/cell_origin_NP/output_19_12_2023/sc_atac/processing/', lib, '_scATAC_isoutlier_normalize_gene_activity.RDS'))
  sr_sngr <- AddMetaData(sr, metadata=labels)
  count_matrix <- GetAssayData(sr, slot = 'counts')
  saveRDS(count_matrix,paste0(save_path,"/", lib, "_sr_count_matrix.RDS") )
  df <- data.frame(cell = colnames(sr_sngr), type = sr_sngr@meta.data$singleR_labels)
  
  write.table(df, file = paste0(save_path,"/", lib, "_sr_cell_annotation.txt"), sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE)
  
}


prep_infercnv_plan <- drake_plan(
  count_matrix = target(make_anno_count_matrix(lib),
                        transform = map(lib = !!lbLst,
                                        id_var = !!idLst,
                                        .id = id_var))
)
make(prep_infercnv_plan,lock_envir = FALSE, lock_cache = FALSE)

vis_drake_graph(prep_infercnv_plan, targets_only = TRUE, lock_cache = FALSE, file = 'prep_infercnv_pipeline.png', font_size = 20 )

