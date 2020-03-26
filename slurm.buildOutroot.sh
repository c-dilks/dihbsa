#!/bin/bash

datadir="/w/hallb-scifs17exp/clas12/rg-a/trains/v16_v2/skim5_inclusiveHadron"

exe="buildOutroot"
slurm=slurm.${exe}.bat
> $slurm

function app { echo "$1" >> $slurm; }

nruns=$(ls ${datadir}/*.hipo | wc -l)
let nruns--

app "#!/bin/bash"
app "#SBATCH --job-name=dihbsa_${exe}"
app "#SBATCH --account=clas12"
app "#SBATCH --partition=production"
app "#SBATCH --mem-per-cpu=1000"
app "#SBATCH --time=3:00:00"
app "#SBATCH --array=0-${nruns}"
app "#SBATCH --ntasks=1"
app "#SBATCH --output=/farm_out/%u/%x-%j-%N.out"
app "#SBATCH --error=/farm_out/%u/%x-%j-%N.err"
app "dataList=(${datadir}/*.root)"

app "srun ${exe}.exe \${dataList[\$SLURM_ARRAY_TASK_ID]}"

echo "job script"
printf '%70s\n' | tr ' ' -
cat $slurm
printf '%70s\n' | tr ' ' -
echo "submitting to slurm..."
sbatch $slurm
squeue -u `whoami`
