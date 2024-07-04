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
mrgprepinfer_plan <- drake_plan(
  # prep for infercnv
  mrgAtacDim = target(sc_atac_dim_redu(mrgAtacNor)),
  mrgGA = get_gene_activity(mrgAtacDim),
  preInferMrg = make_anno_count_Mrgmx(mrgGA, save_path=AtacInferInputDir)
)

make(mrgprepinfer_plan, lock_cache = FALSE)
vis_drake_graph(mrgprepinfer_plan, targets_only = TRUE, lock_cache = FALSE, file = 'vis_mrgprepinfer_plan.png', font_size = 20 )
