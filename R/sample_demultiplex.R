get_h5_link <- function(lb, metadata){
  link <- unique(metadata$data_link[which(metadata$name == lb)])
  h5Link <- paste0(base_data_dir, '/', link, 
                   '/outs/filtered_feature_bc_matrix.h5')
  h5LinkDf <- data.frame(link_to_h5 = h5Link,
                        lb = lb)# needed for later
}

create_GEX_seurat <- function(h5LinkDf){
  message ('-----------------creating rna seurat----------------')
  
  counts <- Read10X_h5(h5LinkDf$link_to_h5)
  sr <- CreateSeuratObject(counts = counts[["Gene Expression"]], 
                             assay = "RNA")
  sr$barcodes <- colnames(sr)
  sr$library <- h5LinkDf$lb
  return(sr)
}

remove_stress_genes <- function(sr, gene_list){
  message ('-----------------removing stress gene from rna seurat---')
  gene_to_retain <- setdiff(rownames(sr), gene_list)
  message (paste(length(setdiff(rownames(sr), gene_to_retain)), 
                 'stress genes were removed'))
  sr_rm <- subset(sr, feature = gene_to_retain)
  return(sr_rm)
}

normalize_dim_sr <- function(sr){
  message ('-----------------normalizing rna seurate----------------')
  sr <- NormalizeData(sr, normalization.method = "LogNormalize")
  sr <- ScaleData(sr, features = rownames(sr), verbose = FALSE)
  sr <- FindVariableFeatures(sr)
  sr <- SCTransform(sr,verbose = FALSE, variable.features.n = 3000)
  sr <- RunPCA(sr, verbose = TRUE, npcs = 50)
  sr <- RunUMAP(sr, dims = 1:30, n.neighbors = 30)
  return(sr)
}

calculate_sex_genes_module <- function(sr){
  message ('-----------------calculating sex genes module-------------')
  message('---set default assay-----')
  DefaultAssay(sr) <- "RNA"
  message('-----addmodulescore female --------')  
  sr <- AddModuleScore(sr, features = list(c("XIST", "TSIX")),
                       name = "female")
  message('-----addmodulescore male --------')
  sr <- AddModuleScore(sr, features = list(male.genes), 
                        name = "male")
  return(sr)
}

add_sex_metadata <- function(sr){
  message ('-----------------adding sex metadata----------------')
  sex <- list()
  for (i in 1:nrow(sr@meta.data)){
    if(sr$female1[i] > 0 & sr$male1[i] < 0){
      sex[i]<-"female"      
    } else {
      if(sr$female1[i] < 0 & sr$male1[i] > 0){
        sex[i]<-"male"
      } else {
        if(sr$female1[i] > 0 & sr$male1[i] > 0){
          sex[i]<-"doublets"
        } else {
          sex[i]<-"unknown"
        }
      }
    }}
  sex <- as.vector(unlist(sex))
  sr$sex <- sex
  return(sr)
}

get_link_to_souporcell <- function(lb, meatadata){
  link <- unique(metadata$souporcell_link[which(metadata$name == lb)]) 
  fullLink <- paste0(base_data_dir, '/',link, '/clusters.tsv')
}

assign_metadata_from_souporcell <- function(sr, metadata){
  message ('-----------------assigning metadata from souporcell-------')
  lb <- unique(sr$library) 
  link_to_souporcell <- get_link_to_souporcell(lb, metadata)
  soc <- read.delim(file=paste0(link_to_souporcell)) 
  genotype <- with(soc, ifelse(status=="singlet", 
                               as.character(assignment), status)) 
  names(genotype) <- soc$barcode 
  sr <- AddMetaData(sr, metadata=genotype, col.name="genotype") 
  return(sr)
}

visualize_soc_gender_demulti <- function(sr, figSavePath){
  lb <- unique(sr$library)
  p <- FeaturePlot(sr, features = 'male1', reduction = "umap", 
                   order = T, min.cutoff = 0) | 
    FeaturePlot(sr, features = 'female1', reduction = "umap", 
                order = T, min.cutoff = 0)| 
    DimPlot(sr, group.by = 'genotype', reduction = "umap", order = T)
  
  savePlot(paste0(figSavePath,'/', lb,"_samples_demultiplex.png"), p)
}

# assign gender for final assessment based on combined res 
# from souporcel and gender 
# - sex == 'unknown' & genotype == 'unassigned' --> remove
# - sex == 'doublets' | genotype == 'doublet' --> remove
# - identity 1,0 to a gender based on which one get the highest count
# - if female 0 and female 1 got the same count, something is wrong, 
# need to check manually
# - contradict prediction sex == 'female' & genotype == male_genotype 
# --> remove
# - sex == 'male' & genotype == female_genotype --> remove
# - assign gender for unassigned souporcell 
# - assign gender for unknown gender 

remove_doublet_n_unknown <- function(sr){
  unknown_cells <- which(sr$sex == 'unknown' & sr$genotype == 'unassigned')
  db_sex <- which(sr$sex == 'doublets')
  db_geno <- which(sr$genotype == 'doublet')
  torm_cells <- unique(colnames(sr)[c(unknown_cells,db_sex, db_geno)])  
  to_keep <- setdiff(colnames(sr), torm_cells)
  sr_rm <- subset(x = sr, subset = barcodes %in% to_keep)
  lb <- unique(sr$library)
  write.csv(torm_cells, paste0(report_dir, '/', lb, '_doublets_n_unkown_cells_to_remove.csv'))
  return(sr_rm)
}

identify_gender_for_souporcell <- function(sr){
  message ('------------combining souporcell and gender-----------')
  lb <- unique(sr$library)
  # identify gender for 0, 1 in souporcell
  cells <- colnames(sr)
  sex <- sr$sex
  genotype <- sr$genotype
  count_0 <- length(which(sex == 'female' & genotype == 0))
  count_1 <- length(which(sex == 'female' & genotype == 1))
  female_genotype <- c()
  male_genotype <- c()
  if (count_0 > count_1){
    female_genotype <- 0
    male_genotype <- 1
  } else if (count_0 < count_1){
    female_genotype <- 1
    male_genotype <- 0
  } else { 
    message (paste('female count 0 and 1 are equal, 
                   manual inspection needed for library', lb))}
  
  genoGendDf <- data.frame(female_genotype = female_genotype,
                                   male_genotype = male_genotype)
  return(genoGendDf)
}

fix_unknown_gender_or_souporcell <- function(sr){
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
  sex[unknown_sex_fgenotype_index] <- 'female'
  unknown_sex_mgenotype_index <- which(sex == 'unknown' & 
                                 genotype == genoGendDf$male_genotype)
  unknown_sex_mgenotype_cells <- cells[unknown_sex_mgenotype_index]
  sex[unknown_sex_mgenotype_index] <- 'male'
  
  sr <- AddMetaData(sr, sex, col.name = 'curate_gender')
  sr <- AddMetaData(sr, genotype, col.name = 'curate_genotype')
  return(sr)
}

remove_contradict_cells <- function(sr){
  genoGendDf <- identify_gender_for_souporcell(sr)
  fS_mG <- which(sr$sex == 'female' & sr$genotype == genoGendDf$male_genotype)
  mS_fG <- which(sr$sex == 'male' & sr$genotype == genoGendDf$female_genotype)
  to_remove <- unique(colnames(sr)[c(fS_mG, mS_fG)])
  to_keep <- setdiff(colnames(sr), to_remove)
  sr_rm <- subset(x = sr, subset = barcodes %in% to_keep )
  lb <- unique(sr$library)
  write.csv(to_remove, paste0(report_dir, '/', lb, '_contradict_cells_to_remove.csv'))
  return(sr_rm)
}  

assign_sample_name_and_tumor_type <- function(metadata, sr){
  message ('----------assigning sample and tumor type-------------')
  lb <- unique(sr$library)
  subdata <- metadata[which(metadata$name == lb),
                c('sampleID', 'Subtype', 'Sample.type', 'Gender')]
  meta_data <- sr@meta.data
  meta_data$order <- 1: nrow(meta_data)
  mergMeta <- merge(meta_data, subdata, 
                           by.x = 'curate_gender', 
                           by.y = 'Gender', all.x = T)  
  mergMetaOrder <- mergMeta[order(mergMeta$order),]
  mergMetaSub <- mergMetaOrder[,c('sampleID', 'Subtype', 'Sample.type')]
  # save statistic
  rnaDemulStat <- data.frame(total_cells = nrow(mergMetaSub),
   unmappedCells = length(mergMetaSub$sampleID[is.na(mergMetaSub$sampleID)]),
   ncellsS1 = length(which(mergMetaSub$sampleID == subdata$sampleID[1])),
   ncellsS2 = length(which(mergMetaSub$sampleID == subdata$sampleID[2])),
   S1 = subdata$sampleID[1],
   S2 = subdata$sampleID[2],
   subtype1 = subdata$Subtype[1],
   subStype2 = subdata$Subtype[2])
  write.csv(rnaDemulStat, paste0(report_dir, '/', lb,
                                 '_scrna_demultiplexing_statistic.csv'))  
  # mark cells that are not mapped to samples. 
  # these cells can be used for further investigation if needed
  mergMetaSub$sampleID[is.na(mergMetaSub$sampleID)] <- 
    paste0(subdata$sampleID[1], '_', subdata$sampleID[2])
  mergMetaSub$Subtype[is.na(mergMetaSub$Subtype)] <- 
    paste0(subdata$Subtype[1], '_', subdata$Subtype[2])
  mergMetaSub$Sample.type[is.na(mergMetaSub$Sample.type)] <- 
    paste0(subdata$Sample.type[1], '_', subdata$Sample.type[2])
  rownames(mergMetaSub) <- rownames(sr@meta.data)
  sr <- AddMetaData(sr, mergMetaSub)
  return(sr)
}
