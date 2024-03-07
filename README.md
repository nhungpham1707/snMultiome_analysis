# scATAC_scRNA analysis pipeline
## Data:
### scATAC:

1. filter matrix.h5. Data was preprocessed with cellranger according to 10X protocol (fastq qc, trimming, alignment, duplication removal, peak calling)

2. fragment file
### scRNA:

filter matrix.h5 output from cellranger

### data for sample demultiplexing:
clsuter.tsv file from souporcell (k=2) 
## Prepare a metadata file

The file should have following columns with these exact names (case sensitive): 
  
  | **column name**| **description**|
  |------------|------------|
  |**RNA_lib**|ID of scRNA as the same with the name of the folder containing .h5 files|
  |**ATAC_lib**| ID of scATAC as the same with the name of the folder containing .h5 and fragment files|
  |**data_link**| link to the h5 file. for example, if the path to filter matrix is LX049_LX050/an_127/outs/filtered_feature_bc_matrix.h5, data_link will be LX049_LX050/an_127|
  |**name**| replace '/' in data_link with '_'. for ex, LX049_LX050_an127|
  |**source**| metadata for the sample|
  |**sampleID**| ID of the sample (can be different from library ID)|
  |**Subtype**| tumor type|
  |**Gender**| patient gender|
  |**Souporcell_link**| link to souporcell result, for ex LX049_LX050/an_207. BY default the function will look for LX049_LX050/an_207/k2/clusters.tsv file|  

    #' @return RDS of seurat object in sc_atac_preprocessing_dir, report in report_dir and plots in sc_atac_preprocessing_fig_dir defined in global_variables.R
 
  An example of the metadata file can be found below


<p>
    <em>Metadata file example. columns are separated by ';'<em>.
        </p>
        
### Data folder structure:

data_link/outs/filtered_feature_bc_matrix.h5

data_link/outs/atac_fragments.tsv.gz

souporcell_link/k2/clusters.tsv  
## Tools:

The main packages are Seurat, Signac and an in-house package SCutils that can be downloaded from bitbucket for curated stress and male gene list. A similar conda environment can be reproduced from the yaml file using the code below:
```
conda env create -f export_environment.yaml
```
Activate conda environment by: 
```
conda activate scRNA_scATAC
```

## Steps
Data was processed according to published procedure  https://stuartlab.org/signac/articles/pbmc_vignette 
The overall steps are described in Fig below
<p>
<img src="https://github.com/nhungpham1707/clean_code_bu/blob/main/github_fig/pipeline.jpg" width="300" alt>
</p>
<p>
    <em>Module 1. Data processing and alignment<em>.
        </p>
### scATAC processing:
        - annotation
       - calculate qc metrics
       - filtering
       - normalization
       - dimension Reduction
       - clustering
       - umap
