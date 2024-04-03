#!/bin/bash
#SBATCH --job-name=sysvi20
#SBATCH --output=log_sysvi_20epo.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --gres=tmpspace:10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl
source /hpc/pmc_drost/nhung/anaconda3/envs/infercnv_mamba/etc/profile.d/conda.sh

# process, merge and prep for infercnv
conda activate scRNA_scATAC_env_copy
python sysvi_20epoch.py