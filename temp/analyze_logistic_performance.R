# analyze logistic regression model performance 
# Nhung 23 08 2024 
# load packages and custom-made functions 
setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu/')
source('_drake.R')

# load data -----
# cells with unknown or ? in names or unclear annotation 
# in DESCARTES were removed. 
# To reduce multicolinearity effects, Cells of closely related classes are grouped, i.e. lymphoid,
# myeloid, ATC were grouped to immune cells.

# load clean dsc ----
cleanDscAtac <- readRDS('output/healthy_data/groupCelltype_nCleanDscAtac.RDS')
# cleanDscAtac were splitted, 80% to train model, 20% to test model 

## train data ----
loadd(sub_cleanDscAtac80) # loadd is a drake function to load a drake target, check _drake.R for pipeline 
# saveRDS(file = 'output/logistic_regression/traindata_cleanDscAtac80.RDS', sub_cleanDscAtac80)
# load train data -----
loadd(sub_cleanDscAtac80)
# load train model ----
loadd(train_cleanDscAtac)
# saveRDS(file = 'output/logistic_regression/train_cleanDscAtac.RDS', train_cleanDscAtac)
## test data ----
loadd(sub_cleanDscAtac20)
# saveRDS(file = 'output/logistic_regression/testdata_cleanDscAtac20.RDS', sub_cleanDscAtac20)

## test result ----
loadd(p_cleanDscAtac_test20)
# saveRDS(file = 'output/logistic_regression/p_cleanDscAtac_test20.RDS', p_cleanDscAtac_test20)

# predict data ----
loadd(atac_hmGroup)
# predict res ---
loadd(p_cleanDscAtac)
# calcp_cleanDscAtac_test20# calculate metrics -----
rocs <- calROC(train_cleanDscAtac, sub_cleanDscAtac20, y = sub_cleanDscAtac20$group_cell_type )
auc <- calAUC(train_cleanDscAtac, sub_cleanDscAtac20, y = sub_cleanDscAtac20$group_cell_type)
cfm <- calconfusionMx(train_cleanDscAtac, sub_cleanDscAtac20, y = sub_cleanDscAtac20$group_cell_type)
performanceMetric <- calAccuracyMetric(cfm)
plotRocs(filename = 'output/logistic_regression/rocs_cleanDscAtactrain80Test202.png', rocs, auc)


fail_cells <- performanceMetric[is.na(performanceMetric[,1]),4]
fail_cells
good_cells <- base::setdiff(names(rocs), fail_cells)

png(filename = 'output/logistic_regression/dscAtacClean_accuracy.png')
plot(performanceMetric[performanceMetric[,'cell_type'] %in% good_cells,'accuracy'])
dev.off()

png(filename = 'output/logistic_regression/dscAtacClean_sensitivity.png')
plot(performanceMetric[performanceMetric[,'cell_type'] %in% good_cells,'sensitivity'])
dev.off()

png(filename = 'output/logistic_regression/dscAtacClean_specificity.png')
plot(performanceMetric[performanceMetric[,'cell_type'] %in% good_cells,'specificity'])
dev.off()

# find cells with excellent performance 
sensitivity_cutoff <- 0.8
head(performanceMetric)
performanceMetric[performanceMetric[,2] > sensitivity_cutoff,4 ] # 16 cells 

# find cells that are illy performed (either similar to most of the cells or not similar to anything)

lowsensitivity_cutoff <- 0.7
performanceMetric[performanceMetric[,2] < lowsensitivity_cutoff, ] # 3 cells 
# Parietal and chief cells, Granule neurons, Inhibitory neurons, ENS neurons, Granule neurons
# it seems there are still strong multicolinearlity effects, especially for neurons cells. 

# inspect each cell 
# neurons ---
index <- which(rownames(p_cleanDscAtac_test20) == 'Neuroendocrine cells')
plot(p_cleanDscAtac_test20[index,])

# count features in predictive set atac_hmGroup --> fail, require too much memory 

# count features in test set 
data <- atac_hmGroup
cell_type <- data$cell_identity

cells <- unique(cell_type)
featureCount <- matrix(data = NA, nrow = length(cells), ncol = 2)
atacMx <- GetAssayData(data)

for (i in 6:length(cells)){
  message(paste('process', i, 'type'))
  subcells <- colnames(data)[cell_type == cells[i]]
  subMx <- atacMx[,subcells]
  df_non0 <- subMx[apply(subMx[,-1], 1, function(x) !all(x==0)),]
  
  featureCount[i,1] <- cells[i]
  featureCount[i,2] <- nrow(df_non0)
}

featureCount[,2] <- as.numeric(featureCount[,2])
colnames(featureCount) <- c('cells', 'count')
p <- ggplot(featureCount, aes(x = cells, y = count)) + 
  geom_bar(stat = 'identity', fill = 'blue') +
  coord_flip() + 
  xlab('Features')+ ylab('Count')+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1,
                                   size = 8),
        text = element_text(size = 20)) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.background = element_blank()) 

p
cells_count <- as.data.frame(table(cell_type))
colnames(cells_count) <- c('cell_type', 'count')
head(cells_count)
p2 <- ggplot(cells_count, aes(x = cell_type, y = count)) + 
  geom_bar(stat = 'identity', fill = 'lightblue') +
  coord_flip() +
  xlab('Cell types') + ylab('Count')+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 10)) 

p|p2

seurat.obj <- atac_hmGroup
colnames(atac_hmGroup@meta.data)
VlnPlot(seurat.obj, features = c("nFeature_peaks", "nCount_peaks"), ncol = 2)
FeatureScatter(seurat.obj, feature1 = "nCount_peaks", feature2 = "nFeature_peaks", group.by = 'cell_identity') +
  geom_smooth(method = 'lm') 

