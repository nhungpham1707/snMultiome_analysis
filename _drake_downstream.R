setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

## read metadata ----
filename <- '25012024_all_multiome_lib.csv'
metadata <- getData(filename, delim = ',')
ori_metadata <- metadata
specialLib <- c("LX093_LX094_an_163")
specialLibInd <- grep(specialLib, metadata$name)
nospecialMet <- metadata[-specialLibInd,]
# get libraries that only demultiplex 
# with souporcell
soclId <- specialLib %>% map(splitName)
# get multiplex libraries list 
mulLib <- unique(nospecialMet$name[nchar(nospecialMet$souporcell_link) > 0])
# mulLib <- mulLib[1]
mulId <- mulLib %>% map(splitName)
# get single libraries list 
sngLib <- unique(nospecialMet$name[nchar(nospecialMet$souporcell_link) == 0])
# sngLib <- sngLib[1]
snglId <- sngLib %>% map(splitName)
# get all samples to make combine peaks

lbLst <- unique(c(specialLib, mulLib, sngLib)) 
idLst <- lbLst %>% map(splitName)

loadd("mrgAtacNor")
downstream_plan <- drake_plan(
    mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor)),
  mrgPtype = dimplot_w_nCell_label(mrgAtacDim, by = 'Subtype',atacMrgFigDir , col = my_cols2),
  mrgPsID = dimplot_w_nCell_label(mrgAtacDim, by = 'sampleID',atacMrgFigDir , col = my_cols2),
  mrgPlibnCell = dimplot_w_nCell_label(mrgAtacDim, by = 'library',atacMrgFigDir , col = my_cols2),
  mrgP = dimplotnSave(mrgAtacDim, atacMrgFigDir, save_name = 'cluster'),
  mrgPlib = dimplotBynSave(mrgAtacDim, by = 'library',
                           atacMrgFigDir, save_name = 'lib_no_nCell', col = my_cols2))

make(downstream_plan, lock_cache = FALSE)
vis_drake_graph(downstream_plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_downstream_plan.png', font_size = 20 )
