setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

# read healthy data metadata
hthy_dataDir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes'
healthy_metadata <- read.csv(paste0(hthy_dataDir,'/filtered.cell_metadata.for_website.txt.gz'), sep = '\t')
all_hthytissue_list  <- unique(healthy_metadata$tissue)
hthytissue_list  <- all_hthytissue_list 

healthy_plan <- drake_plan(
    hthysr = target(readRDS(paste0(hthy_dataDir, '/', ts, '_filtered.seurat.for_website.RDS')),
                transform = map(ts = !!hthytissue_list ,
                                id.vars = !!hthytissue_list ,
                                .id = id.vars)),
    hthysrFill = target(filter_outliers_healthyAtac(hthysr, healthyFigDir),
                transform = map(hthysr,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    hthyNor = target(sc_atac_normalize(hthysrFill),
                transform = map(hthysrFill,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    hthyDim = target(sc_atac_dim_redu(hthyNor),
                transform = map(hthyNor,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    hthySubset = target(sampling_sr(hthyDim, percent_to_keep = 800, type = 'number'),
                transform = map(hthyDim,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    hthymrg = target(merge_sr_list(c(hthySubset), healthyDir),
                transform = combine(hthySubset,
                            id.vars = !!hthytissue_list ,
                            .id = id.vars)),
    # hthymrg = target(merge_pairwise(c(hthySubset), healthyDir),
    #             transform = combine(hthySubset,
    #                         id.vars = !!hthytissue_list ,
    #                         .id = id.vars)),
    hthymrgNor = sc_atac_normalize(hthymrg),
    hthymrgDim = sc_atac_dim_redu(hthymrgNor),
    hthymrgDimP = dimplot_w_nCell_label(hthymrgDim, by = 'tissue', healthyFigDir)
)

make(healthy_plan, lock_cache = FALSE, memory_strategy = 'autoclean', garbage_collection = TRUE,  lock_envir = FALSE)
