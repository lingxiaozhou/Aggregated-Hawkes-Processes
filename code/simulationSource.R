library(coda)
library(truncnorm)


# workdic <- "~/Github/Aggregated-Hawkes-Processes/"
workdic <- ""    # change if necessary

source(paste0(workdic,"code/MCMC_st.R"))
source(paste0(workdic,"code/MCMC_t.R"))
source(paste0(workdic,"code/update_simulation.R"))
source(paste0(workdic,"code/aggregation.R"))
source(paste0(workdic,"code/initializelabel.R"))
source(paste0(workdic,"code/generateHP.R"))
source(paste0(workdic,"code/MCMC_st.R"))
source(paste0(workdic,"code/combineMCMC.R"))
