library(drake)
loadd(rna_hm)
# manually assemble sample list of Jeff's data. from SFig4 in his paper.
jeff_s <- c('RMS127', 'RMS410', 'RMS012',
            'RMS444', 'RMS000CPU', 'RMS000ETXY',
            'RMS672', 'RMS000FWE', 'RMS000FYP',
            'RMS000EEC', 'RMS000FAN', 'RMS000DJE',
            'RMS000GKV', 'RMS000FLV', 'RMS000CXI',
            'RMS000HQC', 'RMS000HEI',
            'RMS000GRN', 'RMS000HVX') # 19 samples

# 6 samples in my multiome dataset "RMS672"    "RMS000FYP" "RMS000EEC" "RMS000GKV" "RMS000HEI" "RMS000HVX" 

tz_s <- read.csv('/hpc/pmc_drost/rstudio/multiome_terezinha/metadata/multiome_experiment_metadata.csv') # 10 libraries
unique(rna_hm$sampleID[rna_hm$RNA_lib %in% tz_s$RNA_lib]) # all 10 libraries are in my multiome dataset
#  [1] "SS220a"    "ATRT24"    "SS000MVX"  "JD166T"    "RMS000BLE" "RMS672"   
#  [7] "RMS000BXG" "RMS000JPD" "RMS000FOE" "SS000ZKK"  "RMS000FYP" "SS001DUK" 
# [13] "RMS000GKV" "RMS000OFM" "RMS000HEI" "SS000DAZ"  "SS459"     "SS001AAJ"  # 7 sysa samples 