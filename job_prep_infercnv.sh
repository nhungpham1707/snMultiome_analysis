#!/bin/bash
#SBATCH --job-name=prepinfer
#SBATCH --output=prepinfercnv.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
# Rscript _drake_prep_infercnv.R
Rscript prep_infercnv_merge_sr.R 
