#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --mem=6G
#SBATCH --array=1-400
#SBATCH --time=24:00:00               # Time limit hrs:min:sec

module load R/3.6
R CMD BATCH --no-save "--args ${tagg} ${sagg} ${mu} ${alpha} ${SLURM_ARRAY_TASK_ID}" 21_sim.R Output/2_sims/mu0${mu}/alpha0${alpha}/beta1/tagg${tagg}/sagg${sagg}/Out_files/Routput_${SLURM_ARRAY_TASK_ID}.Rout


