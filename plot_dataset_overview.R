filename <- '15042024_add_treatment_metadata.csv'
metadata <- getData(filename, delim = ',')
table(metadata$Age.at.first.diagnosis)
hist(metadata$Age.at.first.diagnosis)
metadata[metadata$Age.at.first.diagnosis == 0,]
plot_dataset(metadata, 'Topography.label',save_name = 'output/figures/dataset/topo.png' )

# subtype
plot_dataset(metadata, 'Subtype', save_name = 'output/figures/dataset/subtype.png')

plot_dataset(metadata, 'Individual.ID', save_name = 'output/figures/dataset/patientID.png')
# gender

patients <- unique(metadata$Individual.ID)
length(patients)
genders <- c()
for (i in 1:length(patients)){
  index <- which(metadata$Individual.ID == patients[i])
  genders[i] <- metadata$Gender[index[1]] 
}

length(genders)
count <- table(genders)
png(filename = 'output/figures/dataset/dataset_overview_genders.png' )
pie(table(genders), col= c('salmon', 'navyblue'), cex = 3, cex.main = 3)
dev.off()

# age
age <- metadata %>% group_by(Age.at.first.diagnosis) %>% tally()
png(filename = 'output/figures/dataset/dataset_overview_age.png', width = 1000, height = 1000 )

pie(age$n, labels = paste(age$Age.at.first.diagnosis, '(', age$n, ')'), 
    main="age", col = my_cols, cex = 3, cex.main = 3)
dev.off()

