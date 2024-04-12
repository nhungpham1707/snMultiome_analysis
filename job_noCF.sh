#!/bin/bash
#SBATCH --job-name=rna_noCF
#SBATCH --output=log_removeCF_rna.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
source /hpc/pmc_drost/nhung/anaconda3/envs/infercnv_mamba/etc/profile.d/conda.sh
conda activate scRNA_scATAC_env_copy
Rscript remove_confound_genes.R 