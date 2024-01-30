setwd('/hpc/pmc_drost/PROJECTS/cell_origin_NP/clean_code_bu')
## Load your packages, e.g. library(drake).
source("./packages.R")
## Load your R files ----
functions_folder <- './R'
list_files_with_exts(functions_folder, 'R') %>%
  lapply(source) %>% invisible()

# read healthy data metadata
baseDir <- '/hpc/pmc_drost/PROJECTS/cell_origin_NP/data/healthy_data_descartes'
healthy_metadata <- read.csv(paste0(baseDir,'/filtered.cell_metadata.for_website.txt.gz'), sep = '\t')
all_tissue_list <- unique(healthy_metadata$tissue)
tissue_list <- all_tissue_list

healthy_plan <- drake_plan(
    hthysr = target(readRDS(paste0(baseDir, '/', ts, '_filtered.seurat.for_website.RDS')),
                transform = map(ts = !!tissue_list,
                                id.vars = !!tissue_list,
                                .id = id.vars)),
    hthysrFill = target(filter_outliers_healthyAtac(hthysr, healthyFigDir),
                transform = map(hthysr,
                            id.vars = !!tissue_list,
                            .id = id.vars)),
    hthyNor = target(sc_atac_normalize(hthysrFill),
                transform = map(hthysrFill,
                            id.vars = !!tissue_list,
                            .id = id.vars)),
    hthyDim = target(sc_atac_dim_redu(hthyNor),
                transform = map(hthyNor,
                            id.vars = !!tissue_list,
                            .id = id.vars)),
    hthymrg = target(merge_sr_list(c(hthyDim), healthyDir),
                transform = combine(hthyDim,
                            id.vars = !!tissue_list,
                            .id = id.vars)),
    hthymrgNor = sc_atac_normalize(hthymrg),
    hthymrgDim = sc_atac_dim_redu(hthymrgNor),
    hthymrgDimP = dimplot_w_nCell_label(hthymrgDim, by = 'tissue', healthyFigDir),
    hthyanchors = target(FindIntegrationAnchors(c(hthyDim), 
                   reduction = 'rlsi'),
                  transform = combine(hthyDim,
                      id.var = !!tissue_list,
                      .id = id.)),
    hthyInt = IntegrateData(anchorset = hthyanchors, dims = 1:50)
)

make(healthy_plan, lock_cache = FALSE)
