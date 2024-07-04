#!/bin/bash
#SBATCH --job-name=350infernoR
#SBATCH --output=log_infercnvnoR350.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=350G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

source /hpc/pmc_drost/nhung/anaconda3/envs/infercnv_mamba/etc/profile.d/conda.sh
conda activate r43_copy
Rscript _drake_infercnv.R