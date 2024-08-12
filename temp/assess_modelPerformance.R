source('_drake.R')
loadd(dsc_atac)
dsc_40k <- readRDS('output/logistic_regression/train_dsc_atac_40k.RDS')
dsc_atac20 <- sampling_sr(dsc_atac,percent_to_keep=20, type = 'percent', class_col = 'cell_type')
p_20 <- predictSimilarity(dsc_40k, GetAssayData(dsc_atac20), classes = dsc_atac20$cell_type, logits = F)


refMat <- GetAssayData(dsc_atac20)
norm_dsc20 <- normData(refMat , lengths = rep(1,nrow(refMat)))

x <- norm_dsc20
y <- dsc_atac20$cell_type
fac = factor(y==names(dsc_40k)[[2]])

cglm <- confusion.glmnet(dsc_40k[[2]], x, fac) # https://stackoverflow.com/questions/69957177/prediction-on-new-data-with-glmnet-and-caret-the-number-of-variables-in-newx-m

train_df <- as.data.frame(x )
colnames(train_df)

# not sure why the normalize data are smaller than the raw data 
# inspect normData func 
refMat = GetAssayData(dsc_atac20)
classes = dsc_atac20$cell_type
lengths=rep(1,nrow(refMat))
maxCells=2000
minCells=100

dat = normData(refMat,lengths)

dim(refMat)
dim(dat) # here they have the same dimension - dat is transpose of refMat 

# try again with dsc_atac 
x <- dat
y <- dsc_atac20$cell_type
i <- 9

cglm_df <- c()
asglm_df <- c()
rocs_s4 <- c()
for (i in 1:length(names(dsc_40k))){
# for (i in 1:4){
fac = factor(y==names(dsc_40k)[[i]])

cglm <- confusion.glmnet(dsc_40k[[i]], x, fac, newoffset = getPopulationOffset(fac))

cglm <- as.data.frame(cglm)
cglm$cell_type <- names(dsc_40k)[[i]]
cglm_df <- rbind(cglm_df, cglm)
# 

asglm <- assess.glmnet(dsc_40k[[i]], x, fac, newoffset = getPopulationOffset(fac))
asglm <- as.data.frame(asglm)
asglm$cell_type <- names(dsc_40k)[[i]]
asglm_df <- rbind(asglm_df, asglm)
# cglm percent correct 0.9984 
rocs <- roc.glmnet(dsc_40k[[i]], x, fac, newoffset = getPopulationOffset(fac))
rocs_s4[[i]] <- rocs
names(rocs_s4)[[i]] <- names(dsc_40k)[[i]]
}

# plot res 
asglm_df2 <- asglm_df
colnames(asglm_df) <- c(1,2,3,4,5,'cell_type')

 df <- asglm_df %>% select(1,2,3,4,5 ,cell_type) %>% tidyr::gather( key = 'variable', value = 'value', -cell_type)

df$variable <- as.numeric(df$variable)

ggplot(df, aes(x = variable, y = value)) + geom_line(aes(color = cell_type)) + theme_minimal(
  base_size = 11,
  base_family = "",
  base_line_size = base_size/22,
  base_rect_size = base_size/22
) + scale_x_discrete(labels=c("1" = "Deviance", "2" = "class",
                              "3" = "auc", '4' = 'mse', '5' = 'mae'))

asglm_df2[asglm_df2$auc == 0.5,] # models do not work very well for some cells 
cell_type_count <- data.frame(table(dsc_atac20$cell_type))
cell_type_count[asglm_df2$cell_type[asglm_df2$auc == 0.5]] # these cells have low count in the dataset 

library(RColorBrewer)
n <- length(unique(df$cell_type))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x = variable, y = value)) + geom_line(aes(color = cell_type)) + theme_minimal(
  base_size = 11,
  base_family = "",
  base_line_size = base_size/22,
  base_rect_size = base_size/22
)


plot(rocs_s4[[2]]$FPR, rocs_s4[[2]]$TPR, type = 'l')

plot(rocs_s4[[1]], type = 'l')
invisible(sapply(rocs_s4, lines, col="grey"))
lines(rocs_s4, lwd = 2,col = "red")

best <- dsc_40k[[i]]$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs, lwd = 2,col = "red")



# plot for all cells 
i <- 2
train_model <- dsc_40k[[i]]
plot(train_model)

reso <- 600
png(file = 'output/logistic_regression/dsc_atac_40klambda.png',width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px" )
par(mfrow=c(5,4))
for (i in 1:20){
    cv <- dsc_40k[[i]]
    # plot(cv) + title(names(dsc_40k)[i])
    plot(cv, main = names(dsc_40k)[i], cex.main = 2, cex.lab = 2, cex.axis = 2, lwd = 4, pch = 20 )
}
dev.off()





## try with rna to be faster ---
loadd(rna_hthymrg_clus)
loadd(predict_rna_allCells)
loadd(train_rna_allCells)

rna80 <-  sampling_sr(rna_hthymrg_clus, 75, type = 'percent', class_col = 'cell_type')
rna_hthymrg_clus$cell_bc <- colnames(rna_hthymrg_clus)

rna20bc <- setdiff(colnames(rna_hthymrg_clus), colnames(rna80))
rna20 <- subset(rna_hthymrg_clus, cell_bc %in% rna20bc)

t_rna80 <- trainModel_Nhung(GetAssayData(rna80), classes = rna80$cell_type, maxCells = 2000)