# data_file_name is the sample_gender_meta_data 
getData <- function(filename, delim = ','){
  message('--reading metadata file----')
  data <- read.csv(paste0(base_data_dir, '/', filename), sep = delim, row.names = NULL)
  return (data)
}



# generate shorter lib name to be visible in drake graph
splitName <- function(name){
  name <- strsplit(name, split = '_')[[1]][1]
}

# save plot ---
savePlot <- function(filename,
                      a_plot) {
  
  ggsave(file = filename,
         width = 1200 * reso/72, 
         height = 700 * reso/72, units ="px", dpi = reso)
  
}