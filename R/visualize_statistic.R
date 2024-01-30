combine_statistic_report <- function(tissue_list, dataPath) {
    statistic_df <- c()
    for (tissue in tissue_list){
        df <- read.csv(paste0(dataPath, '/',tissue, 'outlier_statistic.csv'))      # read the file
        statistic_df <- rbind(statistic_df, df)    # append the current file
    }
    #writing the appended file  
    write.csv(statistic_df,paste0(dataPath, "/healthy_data.csv"), row.names=FALSE, quote=FALSE)
    return(statistic_df)
}

plot_healthy_statistic <- function(df, tissue_list, save_path){
    df <- df[,c('total_cells', 'nOutlier_count', 'nOutlier_blk')]
    df$tissue <- tissue_list
    stackdf <- cbind(df[4], stack(df[1:3]))
    ggplot(stackdf, aes(x=tissue, y=values, fill=ind)) +
    geom_bar(stat='identity', position='dodge') + coord_flip() + ggtitle('Number of cells') + scale_fill_manual(values=c('midnightblue', 'sienna1', 'grey')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size = 50)) + xlab('') + ylab('nCells')

    ggsave(filename = paste0(Sys.Date(), "healthy_data.png"), path = save_path,
       width = 1200 * reso/72, 
       height = 700 * reso/72, units ="px", dpi = reso)
}

