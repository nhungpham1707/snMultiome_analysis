# in DESCARTES data, there are cells with uncertain annotation 
# (i.e. with unknown, ? or positive in the name)
# from the prelimenary result, these cells did not perform very well. 
# hence, they will be removed from the final logistic regression training
# in addition, to reduce multicolinearlity problem, cells with closely related class
# are grouped together

# cell_file is csv file, sep by ; indicate which cells to group or to remove 
# example /hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/output/healthy_data/dsc_cell_type_count_group.csv
groupNremoveCellDsc <- function(dsc_seurat, cell_file){
group_type <- read.csv(cell_file, sep = ';')
cell_to_keep <- group_type[group_type$Group != 'Remove',]
cleanDsc <- subset(dsc_seurat, subset = cell_type %in% cell_to_keep[,1] )

# group cell type
metadata <- data.frame(original_name = cleanDsc$cell_type,
                       bc = colnames(cleanDsc))

metadat2 <- merge(metadata, group_type, by = 'original_name')
meta_toadd <- metadat2$Group
names(meta_toadd) <- metadat2$bc
cleanDsc <- AddMetaData(object = cleanDsc, metadata = meta_toadd, col.name = 'group_cell_type')
return (cleanDsc)
}
