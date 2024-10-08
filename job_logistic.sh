#!/bin/bash
#SBATCH --job-name=4rnalogistic
#SBATCH --output=log_logistic_rna4.out
#SBATCH --time=600:0:0
#SBATCH --ntasks=1
#SBATCH --mem=350G
#SBATCH --cpus-per-task=1
#SBATCH --gres=tmpspace:10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
source /hpc/pmc_drost/nhung/anaconda3/envs/infercnv_mamba/etc/profile.d/conda.sh

# process, merge and prep for infercnv
conda activate scRNA_scATAC_env_copy

# Rscript _drake.R
Rscript _drake_logistic.R

