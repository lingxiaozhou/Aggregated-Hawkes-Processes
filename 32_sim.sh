#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --mem=6G
#SBATCH --array=1-400
#SBATCH --time=24:00:00               # Time limit hrs:min:sec

module load R/3.6
R CMD BATCH --no-save "--args ${agg1} ${agg2} ${SLURM_ARRAY_TASK_ID}" Output/3_sims/1_agg${agg1}/2_agg${agg2}/Out_files/Routput_${SLURM_ARRAY_TASK_ID}.Rout


