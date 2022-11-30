Bayesian Inference for Aggregated Hawkes Processes
================

This GitHub repository contains the code and results associated with the Section 4 of the manuscript.


## Instructions for use

The provided run simulations for all the scenarios that we considered (21200 simulations in total) and reproduce all the figures in the Section 4 of the manuscript.
Access to the cluster is required and it takes approximately 5 hours to run one simulations in the most complex scenario. The steps for running the code are as follows:

1.  Download the code folder and result folder, and transfer them to the cluster
2.  Execute setup.txt to create folders for generated data and MCMC results
3.  Execute submit.txt to submit jobs to the cluster and run simulations
4.  Run combineResults.R to read and combine all simulation results
5.  Run simulation_plots.R to reproduce all the figures in Section 4 of the manuscript.

