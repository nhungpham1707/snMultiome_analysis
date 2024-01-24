#!/usr/bin/env Rscript

library(infercnv)
save_path <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/output/cell_type/sc_atac/infercnv/merge_sr/input/'

# save_path <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/output_19_12_2023/cell_type/sc_atac/infercnv/input'

normal_cells <- c('B_cell', 'T_cells', 'Macrophage', 'Monocyte', 'NK_cell')
lib <- 'LX049_LX050_an_127'
# cell_annotation_link <- paste0(save_path,"/", lib, "_sr_cell_annotation.txt")
cell_annotation_link <- paste0(save_path,"/sr_w_dj_w_fail_demul_annotation.txt")
link_to_gene_order_file = '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/analyses/analyses/gencode_v19_gene_pos.txt'
# count_matrix <-  readRDS(paste0(save_path,"/", lib, "_sr_count_matrix.RDS"))
count_matrix <- readRDS(paste0(save_path,"/sr_w_dj_w_fail_demul_count_matrix.RDS"))
output_link <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code/output/cell_type/sc_atac/infercnv/merge_sr'
gene_order_file_delim = '\t'
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= count_matrix,
                                    annotations_file= cell_annotation_link,
                                    gene_order_file=link_to_gene_order_file,
                                    ref_group_names= normal_cells)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.3, 
                             out_dir= output_link, 
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=F,
                             analysis_mode = "samples",
                             num_threads = 3,
                             output_format = "pdf",
                             window_length = 101,
                             save_final_rds = T,
                             plot_steps=F,
                             sd_amplifier = 2)

message('finish infercnv without error')