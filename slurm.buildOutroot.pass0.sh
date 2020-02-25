#!/bin/bash

#SBATCH --job-name=dihbsa_buildOutroot
#SBATCH --account=clas12
#SBATCH --partition=production

#SBATCH --mem-per-cpu=3000
#SBATCH --time=1:30:00

#SBATCH --array=0-8
#SBATCH --ntasks=1

#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err

dataList=(/lustre19/expphy/volatile/clas12/rg-a/production/pass0/physTrain/dst/train/skim4/*)

srun buildOutroot.exe ${dataList[$SLURM_ARRAY_TASK_ID]}

