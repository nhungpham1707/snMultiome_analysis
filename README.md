# scATAC_scRNA analysis pipeline
This repo contains code to process and analyze single-nucleus (sn) multiome data (snRNAseq + snATACseq). The analysis includes:
- QC filter 
- Normalization, Dimension reduction
- Clustering 
- Inspecting clustering quality
- Cell type annotation
- Removing patient-specific effects
- Comparing with published snRNAseq and snATACseq from normal embryonal tissues to identify potential cancer cell of origin

In this pipeline, drake, a workflow manager for R was used. drake helps to store intermediate outputs as cache and can be loaded if needed. In addition, drake helps to rerun the pipeline and skips steps that are still up-to-date. So no need to comment out steps. For more information about drake please visit the [github](https://github.com/ropensci/drake)

The main scripts are:
-  _drake.R: a master script that run all analysis steps except infercnv 
- _drake_infercnv.R: master script to run infercnv (getting input data from _drake.R)
- 
All codes were run on HPC with SLURM, using 350Gb memory for 1Tb data. Jobs was called by job.sh

Before running the code, modify the data directory in global_variables.R and global_variables_infercnv.R 

## Data:
### scATAC:

1. filter matrix.h5. Data was preprocessed with cellranger according to 10X protocol (fastq qc, trimming, alignment, duplication removal, peak calling)

2. fragment file
### scRNA:

filter matrix.h5 output from cellranger

### data for sample demultiplexing:
clsuter.tsv file from souporcell (k=2 because there are 2 samples in each multiplex library) 
## Prepare a metadata file

The file should have the following columns with these exact names (case sensitive): 
  
  | **column name**| **description**|
  |------------|------------|
  |**RNA_lib**|ID of scRNA library. Can be empty|
  |**ATAC_lib**| ID of scATAC library. Can be empty|
  |**data_link**| link to the h5 file. for example, if the path to a filter matrix for a library A is libraryA/an_127/outs/filtered_feature_bc_matrix.h5, data_link will be libraryA/an_127|
  |**name**| replace '/' in data_link with '_'. for ex, libraryA_an127. This will be used to label clusters by library. Can also be any custom name|
  |**source**| metadata for the sample. for instance, tissue location|
  |**sampleID**| ID of the sample for each library (can be different from library ID). This infomation is used for confounding correction|
  |**Subtype**| cancer type. for instance, 'MRTK', 'ATRT'|
  |**Souporcell_link**| link to souporcell result, for ex libraryA/an_207. By default the function will look for libraryA/an_207/k2/clusters.tsv file| 
  |**Date.of.Library**| sequencing date of each library. This infomation is used for confounding correction|
  |**Individual.ID**| patients ID for each sample. This infomation is used for confounding correction| 
  |**Gender**| patient genders for each sample. This infomation is used for confounding correction|


    #' @return RDS of seurat object in sc_atac_preprocessing_dir, report in report_dir and plots in sc_atac_preprocessing_fig_dir defined in global_variables.R
 
  An example of the metadata file can be found below
<p>
<img src="https://github.com/nhungpham1707/scATAC_scRNA/blob/main/example_metadata_file.png" alt>
</p>

<p>
    <em>Metadata file example. columns are separated by ';'<em>.
        </p>
        
### Data folder structure:

data_link/outs/filtered_feature_bc_matrix.h5

data_link/outs/atac_fragments.tsv.gz

souporcell_link/k2/clusters.tsv  
## Tools:

The main packages are Seurat, Signac and an in-house package SCutils that can be downloaded from [bitbucket](https://bitbucket.org/princessmaximacenter/scutils/src/master/) for curated stress and male gene list. A similar conda environment can be reproduced from the yaml file using the code below:
```
conda env create -f multiome_env.yml
```
Activate conda environment by: 
```
conda activate scRNA_scATAC_env
```
 
There is a conflict to install infercnv and Seurat in the same environment. Therefore, infercnv was installed in a separate environment. For some reason, infercnv cannot be installed with conda or in R in a conda env. Fortunately, infercnv can be installed with mamba following instruction [here](https://gist.github.com/hiraksarkar/28824d9943309544a454c595ac0441f7)

In short, 
```
mamba create -n r43
mamba activate r43
mamba install r-essentials=4.3

export PKG_CONFIG_PATH=/home/user/miniconda3/envs/r43/lib/pkgconfig/:$PKG_CONFIG_PATH

# replace /home/user/miniconda3/ to the path of your conda or miniconda 

# activate env 
conda activate r43 

# start R
R

# Install R package
install.packages("rjags")
install.packages("BiocManager")
BiocManager::install("infercnv")
```

## Steps
Data was processed according to published procedure  https://stuartlab.org/signac/articles/pbmc_vignette 
The overall steps are described in Fig below
<p>
<img src="https://github.com/nhungpham1707/clean_code_bu/blob/main/github_fig/pipeline.drawio.png" width="450" alt>
</p>
<p>
    <em>A schematic of the analysis pipeline<em>.
        </p>

