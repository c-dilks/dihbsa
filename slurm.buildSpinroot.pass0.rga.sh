#!/bin/bash

#SBATCH --job-name=dihbsa_buildSpinroot
#SBATCH --account=clas12
#SBATCH --partition=production

#SBATCH --mem-per-cpu=3000
#SBATCH --time=1:30:00

#SBATCH --array=0-8
#SBATCH --ntasks=1

#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err

dataList=(${PWD}/outroot/*.root)

srun buildSpinroot.exe -f ${dataList[$SLURM_ARRAY_TASK_ID]} -i2 -t2 -l1 -m1 -b

