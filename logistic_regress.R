# try logistic regression to find cell origin 
# 

loadd(hthyDim_adrenal)
ref_sr <- hthyDim_adrenal
loadd(rna_w_tumor_label)
refMat <- GetAssayData(ref_sr, slot = 'counts')
classes <- ref_sr$cell_type
t_m <- trainModel(refMat,classes)
tgtData <- GetAssayData(rna_w_tumor_label, slot = 'counts')
p_classes <- rna_w_tumor_label$cell_identity
p_m <- predictSimilarity(t_m,tgtData,classes= p_classes,minGeneMatch=0.70)
png(filename = 'output/figures/logistic_regression/trained_kidney.png')
similarityHeatmap(p_m)
dev.off()
# endothelial cells from our scRNAseq are very similar to the vascular endothelial cells in the reference, which is a good positive control 

# test on 75% ref data to train and test on the rest 25% of the ref data 
# for each cell type select 75% cells randomly 
sub_sr <- sampling_sr(ref_sr, 75)
sub_data <- GetAssayData(sub_sr, slot = 'counts')
sub_classes <- sub_sr$cell_type  
t_sub_sr <- trainModel(sub_data, sub_classes)

# test on the remaining 25% cells, if the model was trained correctly, it should predict correctly the label of these cells
ref_sr$all_bc <- colnames(ref_sr)
to_keep <- setdiff(colnames(ref_sr), colnames(sub_sr))
test_sr <- subset(ref_sr, subset = all_bc %in% to_keep)
test_classes <- test_sr$cell_type
test_data <- GetAssayData(test_sr, slot = 'counts')
test_m <- predictSimilarity(t_sub_sr,test_data,classes= test_classes,minGeneMatch=0.70)
png(filename = 'output/figures/logistic_regression/test_adrenal.png')
similarityHeatmap(test_m)
dev.off()
# yes, cell types are correctly correlated to themselve again --> model seems to be fine 

# generate merg reference dataset from descartes
# for each tissue, get a maximum 800 cells per cell type, then merg all the tissue together 

sub_sr <- sampling_sr(sr, percent_to_keep = 800, type = 'number') 



## tr Jeff reference from xi 2020 
# try on openondemand 
 hthy <- readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Jeff_rf/xi_2020.rds')
 hthy = normalize_dim_plot_sr(hthy, save_path = healthyDir, lib_name = 'xi_2020' )
hthy = clustering_rna_data(hthy)

xi_train <- trainModel(GetAssayData(hthy), classes = hthy$cell_type)

p_xi <- predictSimilarity(xi_train,GetAssayData(rna),classes= rna$cell_identity,minGeneMatch=0.70, logits = F)


similarityHeatmap(p_xi)
# xu at alas from Terezinha 
 xu_atlas = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/Terezinha_reference/xu_atlas_2023.RDS')
train_xu_atlas = trainModel(GetAssayData(xu_atlas), 
            classes = xu_atlas$final_annotation, maxCells = 2000)
predict_xu_atlas = predictSimilarity(train_xu_atlas, 
        GetAssayData(rna), 
        classes = rna$cell_identity,
        minGeneMatch = 0.7)

probCols = brewer.pal(n = 9,name = 'RdYlBu')
probCols = rev(probCols)
probCols[1] <- 'white'

heatmap_only_significant_prob(p_ka_40k, prob_cutoff  = 0.6, probCols = probCols)


## make atac ga for descartes 
desc_eye <- readRDS()
annotation <- getHg38Annotation()
fragpath <- c('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/fragment_files/sample_55_eye.fragments.txt.gz', '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/fragment_files/sample_53_eye.fragments.txt.gz')

fragpath = '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/fragment_files/mrg_eye.txt.gz'
sr_chr <-  CreateChromatinAssay(
    counts = GetAssayData(sr),
    sep = c(":", "-"),
    annotation = hg38,
    min.cells = 10,
    min.features = 200)


sr_chr <- CreateSeuratObject(counts = sr_chr,
                                assay = "peaks" )
eyeSr_ga <- get_gene_activity(eyeSr)


# create index frag file 
tabix -p bed sample_55_eye.fragments.txt.gz 
 
in bash in conda activate /hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_env/

# create similar peak coordinate between test and training set 
# make granges dataframe from rownames of atac sr
sr <- atac_group_sgr
make_bedfile_from_rownamesSr <- function(sr){
coord_df <- c()
for (i in 1:length(rownames(sr))){
        message(paste('i is', i))
coord <- strsplit(rownames(sr)[i], split = '-')
coord_df <- data.frame(chr = coord[[1]][1],
                        start = coord[[1]][2],
                        end = coord[[1]][3]) %>% rbind( coord_df, . )
}
return(coord_df)
}

sr2 <- atac_hthyDim_eye

gr_atac <- makeGRangesFromDataFrame(coord_df_atac)
gr_eye <- makeGRangesFromDataFrame(coord_df)
allpeaks <- disjoin(c(gr_atac, gr_eye))
allpeaks2 <- GenomicRanges::reduce(c(gr_atac, gr_eye))
library(Repitools)

a_df <- annoGR2DF(allpeaks)

# new idea: merg hthy atac with my atac, it will combine the peak list 
atac <- hmAtac_scroshi_atrt
atac$source <- 'nhung_etal'
hthy = readRDS('/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes/all_celltypes.downsampled.filtered.RDS')hthy$source = 'descartes'
dsc_sr <- CreateChromatinAssay(counts = 
                GetAssayData(hthy), 
                sep = c(":","-"), 
                annotation = hg38, 
                min.cells = 10, 
                min.features = 200)
dsc_sr <- CreateSeuratObject(counts = dsc_sr,
                                assay = "peaks" )
metadata <- hthy@meta.data
rownames(metadata) <- colnames(hthy)
dsc_sr <- AddMetaData(dsc_sr, metadata)
mrg_data <- merge(x= atac, y = dsc_sr, add.cell.ids = c('this_paper', 'descartes'))
saveRDS(file = 'output/logistic_regression/mrg_descartes_atac.RDS', mrg_data)
 dsc_data <- subset(mrg_data, subset = source == 'descartes')
 atac_data <- subset(mrg_data, subset = source == 'nhung_etal')

# find overlap between atac and dsc ---

dsc_chr <-  CreateChromatinAssay(
    counts = GetAssayData(hthy),
    sep = c(":", "-"),
    annotation = hg38,
    min.cells = 10,
    min.features = 200)


dsc_chr <- CreateSeuratObject(counts = dsc_chr,
                                assay = "peaks" )
metadata <- dsc_chr@meta.data
rownames(metadata) <- colnames(dsc_chr)
dsc_chr <- AddMetaData(dsc_chr, metadata)

atac_gr <- granges(atac_hm_w_tumor_label)
dsc_gr <- granges(hthy)
gr1 <- atac_gr
gr <- dsc_gr
countOverlaps(gr, gr1)
olap <- findOverlaps(gr, gr1)
sub_olap <- subsetByOverlaps(gr, gr1)