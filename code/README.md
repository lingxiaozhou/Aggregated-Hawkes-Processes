This folder contains all necessary files to run all the simulations for three scenarios we considered (temporal process, single spatio_temporal process and multiple spatio-temporal process).  The contents of this folder are as follows:

* aggregation.R: function for aggregating a continous Hawkes process
* combineMCMC.R: function that combines two MCMC chains
* generateHP.R: function for generating Hawkes processes
* initializalabel.R: function for initializing the branching structure
* MCMC_result.R: function that reads and combines all the simulated results
* MCMC_st.R: function that fit the model in the spatio-temporal case
* MCMC_t.R: function that fit the model in the temporal case
* simulation_plots.R: reproduce all the figures in section 4
* simulationSource: function that sources all the necessary functions
* update_simulation: function that updates variables in the MCMC
* 11_sim.R: code for running one simulation for tempral case
* 21_sim.R: code for running one simulation for spatio-temporal case
* 31_sim.R: code for running one simulation for multiple spatio-temporal case
* 12_sim.sh, 22_sim.sh, 32_sim.sh: corresponding shell script files
* setup.txt: shell script for creating folders to store the generated data and results
* submit.txt: shell sccrip for submit jobs to the cluster


Required packages are:

* coda (Version 0.19-4)
* truncnorm (Version 1.0-8)
* ggplot2 (Version 3.3.6)
