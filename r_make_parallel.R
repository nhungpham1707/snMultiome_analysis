# make.R
library(drake)
r_vis_drake_graph(targets_only = TRUE, file = 'drake_pipeline_paralel.png', font_size = 20 )

options(
  clustermq.scheduler = "slurm"
  # Created by drake_hpc_template_file("sge_clustermq.tmpl"):
)
make(
  plan,
  parallelism = "clustermq",
  jobs = 8
)