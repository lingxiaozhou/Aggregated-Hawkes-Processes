#!/bin/bash
#SBATCH --account=gpapadogeorgou
#SBATCH --qos=gpapadogeorgou-b
#SBATCH --mem=4G
#SBATCH --array=1-400
#SBATCH --time=6:00:00               # Time limit hrs:min:sec

module load R/3.6
R CMD BATCH --no-save "--args ${agg} ${mu} ${alpha} ${SLURM_ARRAY_TASK_ID}" 11_sim.R /orange/gpapadogeorgou/Lingxiao/Simulation/Output/1_sims/mu0${mu}/alpha0${alpha}/beta1/agg${agg}/Out_files/Routput_${SLURM_ARRAY_TASK_ID}.Rout


