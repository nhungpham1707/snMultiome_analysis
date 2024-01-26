#!/bin/bash
#SBATCH --job-name=cleanHealt
#SBATCH --output=healthy_data.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=4
#SBATCH --gres=tmpspace:10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
Rscript merge_healthy_data.R
