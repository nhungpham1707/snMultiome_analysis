#!/bin/bash
#SBATCH --job-name=infercnv
#SBATCH --output=infercnv.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
Rscript run_infercnv.R
# Rscript _drake_infercnv.R