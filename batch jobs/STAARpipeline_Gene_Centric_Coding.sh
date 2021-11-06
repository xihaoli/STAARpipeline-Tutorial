#!/bin/bash
#SBATCH -J LDL_C
#SBATCH -p shared
#SBATCH --time=0-144:00
#SBATCH --array=1-379 --mem=6000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-type=NONE

module purge
module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01
module load intel-mkl/2017.2.174-fasrc01
module load R/3.6.1-fasrc01

export R_LIBS_USER=$HOME/R-3.6.1-MKL
echo $R_LIBS_USER

/n/home05/zilinli/R-3.6.1/bin/Rscript --slave --no-restore --no-save STAARpipeline_Gene_Centric_Coding.r ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout