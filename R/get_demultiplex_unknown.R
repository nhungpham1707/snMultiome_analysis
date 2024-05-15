# check if those unknown or unassign cells are these weird cells in umap merg (cells with a different subtype from the cluster they are in)
 
 # sr are gexNoDb
get_unknown_unassign_cells <- function(sr){
  lb <- unique(sr$library)
  # identify gender for 0, 1 in souporcell
  cells <- colnames(sr)
  sex <- sr$sex
  genotype <- sr$genotype
  genoGendDf <- identify_gender_for_souporcell(sr)

  unassign_fgenotype_index <- which(genotype == 'unassigned' & 
                                      sex == 'female')
  unassign_fgenotype_cells <- cells[unassign_fgenotype_index]
  
  genotype[unassign_fgenotype_index] <- genoGendDf$female_genotype
  unassigned_male_genotype_index <- which(genotype == 'unassigned' & 
                                            sex == 'male')

  unassigned_male_genotype_cells <- cells[unassigned_male_genotype_index]
  genotype[unassigned_male_genotype_index] <- genoGendDf$male_genotype
  
  unknown_sex_fgenotype_index <- which(sex == 'unknown' & 
                                 genotype == genoGendDf$female_genotype)
  unknown_sex_fgenotype_cells <- cells[unknown_sex_fgenotype_index]

  unknown_sex_mgenotype_index <- which(sex == 'unknown' & 
                                 genotype == genoGendDf$male_genotype)
  unknown_sex_mgenotype_cells <- cells[unknown_sex_mgenotype_index]

  all_unknown <- c(unknown_sex_fgenotype_cells, unknown_sex_mgenotype_cells, unassigned_male_genotype_cells, unassign_fgenotype_cells )
  return (all_unknown)
}


get_bc_in_mrg <- function(mrg_sr, bc, lib){
    lib_in_mrg <- mrg_sr$barcodes[grep(lib, mrg_sr$library)]
    bc_in_mrg <- lib_in_mrg[lib_in_mrg %in% bc]
    bc_in_mrg <- names(bc_in_mrg)
    return (bc_in_mrg)
}


plot_unknow_cell <- function(mrg_sr, sr, lib){
    bc <- get_unknown_unassign_cells(sr)
    bc_in_mrg <- get_bc_in_mrg(mrg_sr, bc, lib)
    DimPlot(mrg_sr, cells.highlight = bc_in_mrg)
    return(bc_in_mrg)
}

remove_cells_from_cluster <- function(mrg_sr, cluster, subtype){
    index <- which( mrg_sr$RNA_snn_res.0.8 == cluster)
    cluster_cells <- colnames(mrg_sr)[index]

    to_keep <- setdiff(colnames(mrg_sr), cluster_cells)
    subSr <- subset(mrg_sr, subset = m_barcode %in% to_keep)
}