#!/bin/bash
#SBATCH --job-name=visinfer
#SBATCH --output=visinfercnv.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
Rscript _drake_vis_infercnv.R
