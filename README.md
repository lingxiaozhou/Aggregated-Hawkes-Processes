Bayesian Inference for Aggregated Hawkes Processes
================

This GitHub repository contains the code and results associated with the Section 4 of the manuscript. Thw following files contain main code to run simulations and reproduce figures.


* 11_sim.R: code for running one simulation for tempral case
* 21_sim.R: code for running one simulation for spatio-temporal case
* 31_sim.R: code for running one simulation for multiple spatio-temporal case
* 12_sim.sh, 22_sim.sh, 32_sim.sh: corresponding shell script files
* setup.txt: shell script for creating folders to store the generated data and results
* submit.txt: shell script for submit jobs to the cluster
* combineResults.R: read and combine all simulation results
* simulation_plots.R: reproduce all the figures in section 4

Required packages are:

* coda (Version 0.19-4)
* truncnorm (Version 1.0-8)
* ggplot2 (Version 3.3.6)


The code folder contains all the untility functions, and the MCMC results for simulations are provided in the result folder.


## Instructions for use

The provided run simulations for all the scenarios that we considered (21200 simulations in total) and reproduce all the figures in the Section 4 of the manuscript.
Access to the cluster is required and it takes approximately 5 hours to run one simulations in the most complex scenario. The steps for running the code are as follows:

1.  Download all files, and transfer them to the cluster
2.  Execute setup.txt to create folders for generated data and MCMC results
3.  Execute submit.txt to submit jobs to the cluster and run simulations
4.  Run combineResults.R to read and combine all simulation results
5.  Run simulation_plots.R to reproduce all the figures in Section 4 of the manuscript.

