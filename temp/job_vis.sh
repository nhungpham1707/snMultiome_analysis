#!/bin/bash
#SBATCH --job-name=visinfer
#SBATCH --output=log_visinfercnv.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
source /hpc/pmc_drost/nhung/anaconda3/envs/infercnv_mamba/etc/profile.d/conda.sh
conda activate r43_copy
Rscript _drake_vis_infercnv.R
Rscript _drake_vis_mrg_infer.R

conda activate scRNA_scATAC_env_copy
Rscript _drake_vis_no_relapse.R
Rscript _drake_vis_graph.R
Rscript _drake_vis_healthy.R
echo 'vis conda worked!'
