library(coda)
library(truncnorm)

# the full name of the code folder
funcdict <- "~/Github/Aggregated-Hawkes-Processes/code/"

source(paste0(funcdict,"MCMC_st.R"))
source(paste0(funcdict,"MCMC_t.R"))
source(paste0(funcdict,"update_simulation.R"))
source(paste0(funcdict,"aggregation.R"))
source(paste0(funcdict,"initializelabel.R"))
source(paste0(funcdict,"generateHP.R"))
source(paste0(funcdict,"MCMC_st.R"))
source(paste0(funcdict,"combineMCMC.R"))
