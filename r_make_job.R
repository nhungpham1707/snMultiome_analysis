library(drake)
r_vis_drake_graph(targets_only = TRUE, file = 'drake_pipeline_before.png', font_size = 20 )
r_make() # it will start a new R session and run _drake.R using the drake_config 
# ref https://books.ropensci.org/drake/projects.html