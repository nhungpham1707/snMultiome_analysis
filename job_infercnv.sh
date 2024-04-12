#!/bin/bash
#SBATCH --job-name=rnainfercnv
#SBATCH --output=log_rna_infercnv.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --gres=tmpspace:10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
source /hpc/pmc_drost/nhung/anaconda3/envs/infercnv_mamba/etc/profile.d/conda.sh

# run infercnv 
conda activate r43_copy
Rscript _drake_infercnv.R

# integrate infercnv res 
