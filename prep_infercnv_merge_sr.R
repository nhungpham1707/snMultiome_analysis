#!/usr/bin/env Rscript
library(Signac)
library(Seurat)
singleR_res = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/output/cell_type/sc_atac/singler/singleR_merge_all_w_disjoin_w_fail_demul.RDS')
  labels <- data.frame(singleR_labels = singleR_res$labels,
                     singleR_pruned.labels = singleR_res$pruned.labels)

  rownames(labels) <- rownames(singleR_res)
  sr <- readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/output/sc_atac/merge_all_w_disjoin/merge_sr_w_gene_activity_w_fail_demul.RDS')
  sr_sngr <- AddMetaData(sr, metadata=labels)
  count_matrix <- GetAssayData(sr, slot = 'counts')

  save_path <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/output/cell_type/sc_atac/infercnv/merge_sr/input'
dir.create(save_path, recursive = TRUE)  

saveRDS(count_matrix,paste0(save_path,"/sr_w_dj_w_fail_demul_count_matrix.RDS") )
  df <- data.frame(cell = colnames(sr_sngr), type = sr_sngr@meta.data$singleR_labels)
  
  write.table(df, file = paste0(save_path,"/sr_w_dj_w_fail_demul_annotation.txt"), sep = "\t", quote = FALSE, na = "", row.names = FALSE, col.names = FALSE)