#! /bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=flaq_amr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ENTER EMAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=1-00
#SBATCH --output=flaq_amr.%j.out
#SBATCH --error=flaq_amr.%j.err

module load apptainer

#Run script/command and use $SLURM_CPUS_ON_NODE
python flaq_amr.py fastqs/ 
