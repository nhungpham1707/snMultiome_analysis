#!/bin/bash
#SBATCH --job-name=bavaria
#SBATCH --output=log_bavaria.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --gres=tmpspace:10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


#conda create --name bav_env_py37_2 --file requirements.txt
# conda activate bav_env_py37_2
python 00_preparation.py
