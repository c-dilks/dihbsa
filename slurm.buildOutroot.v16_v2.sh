#!/bin/bash

#SBATCH --job-name=dihbsa_buildOutroot
#SBATCH --account=clas12
#SBATCH --partition=production

#SBATCH --mem-per-cpu=3000
#SBATCH --time=2:00:00

#SBATCH --array=0-63
#SBATCH --ntasks=1

#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err

dataList=(/w/hallb-scifs17exp/clas12/rg-a/trains/v16_v2/skim5_inclusiveHadron/*)

srun buildOutroot.exe ${dataList[$SLURM_ARRAY_TASK_ID]}
