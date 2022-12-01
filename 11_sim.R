args <- commandArgs(TRUE)
agg_size <- as.numeric(args[1])/100

sim <- 1

mu <- as.numeric(paste0('0.', args[2]))
alpha <- as.numeric(paste0('0.', args[3]))
index <- as.numeric(args[4])


beta <- 1

# change the working dictionary if necessary
workdic <- ""

outname <- paste0(workdic,"Output/",sim,"_sims/mu0",mu*10,"/alpha0",alpha*10,"/beta1/agg",agg_size*100,"/Results/")
datapath <-paste0(workdic,"Data_specs/mu0",mu*10,"/alpha0",alpha*10,"/beta1/")

# make sure the full name of the file is correct
source(paste0(workdic,"code/simulationSource.R")) 

cat(mu,"\n",alpha,"\n",beta,"\n",agg_size,"\n")

t <- 500;bin.length <- 1;Total_itr <- 20000; burn <- 10000; epsilon_beta <- ifelse(mu==0.3,0.2,0.3)


if(paste0("temporal_",index,".dat") %in% list.files(outname)){
  cat("Exist!")
}else{
  #------------------------------------simulate data-----------------------------------------
  
  if(paste0("data_",index,".dat") %in% list.files(datapath)){
    cat("Data Exist!")
    load(paste0(datapath,"data_",index,".dat"))
  }else{
    HP <- generateData(seed = index,remove = TRUE,mu = c(mu),alpha = array(alpha,c(1,1)),beta = array(beta,c(1,1)), t = t,gamma = array(1,c(1,1)), w.width = 100,w.length = 100)
    save(HP, file = paste0(datapath,"data_",index,".dat"))
    write.table(HP$HP,file = paste0(datapath,"data_",index,".txt"),row.names=FALSE,col.names=FALSE)
  }
  
  gr <- 5
  
  # intial values for the MCMC
  mu_1 <- 1;alpha_1 <- 1;beta_1 <- 1;label_1 <- c();time_1 <- c()
  mu_2 <- 1;alpha_2 <- 1;beta_2 <- 1;label_2 <- c();time_2 <- c()
  
  while (gr>1.1) {
    chain1 <- MCMC_t(HP$HP,HP$label,HP$type,
                     seed = index,Total_itr=Total_itr,burn = burn,
                     q=0.998,use_q=FALSE,updatelabel_ll=200,portion=1,t=t,bin.length=agg_size*bin.length,
                     mu=rep(mu_1,1),alpha=array(alpha_1,c(1,1)),array(beta_1,c(1,1)), time_0=time_1,label_0=label_1,
                     e_beta=array(epsilon_beta ,c(1,1)),silent=FALSE,summ_itr=5000, save.latent = FALSE,
                     update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.time=TRUE,update.label=TRUE)
    
    chain2 <- MCMC_t(HP$HP,HP$label,HP$type,
                     seed = index+200,Total_itr=Total_itr,burn = burn,
                     q=0.998,use_q=FALSE,updatelabel_ll=200,portion=1,t=t,bin.length=agg_size*bin.length,
                     mu=rep(mu_2,1),alpha=array(alpha_2,c(1,1)),beta=array(beta_2,c(1,1)),time_0=time_2,label_0=label_2,
                     e_beta=array(epsilon_beta ,c(1,1)),silent=FALSE,summ_itr=5000, save.latent = FALSE,
                     update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.time=TRUE,update.label=TRUE)
    
    #------------------------------check convergence of MCMC-----------------------------------------
    
    combined <- combineMCMC(chain1,chain2)
    tmp <- gelman.diag(combined,autoburnin = FALSE,multivariate = TRUE)$psrf[,1]
    gr <- max(tmp)
    
    if(gr>1.1){
      cat("parameter with convergence problem:", which(tmp>1.1),"\n")
      if(chain1$ar>0.5){
        epsilon_beta <- epsilon_beta*2
        mu_1 <- 1;alpha_1 <- 1;beta_1 <- 1;label_1 <- c();time_1 <- c()
        mu_2 <- 1;alpha_2 <- 1;beta_2 <- 1;label_2 <- c();time_2 <- c()
        Total_itr <-20000
        burn <- 10000
      }else{
        mu_1 <- chain1$mu.p[1,length(chain1$alpha.p)]
        alpha_1 <- chain1$alpha.p[1,length(chain1$alpha.p)]
        beta_1 <- chain1$beta.p[1,length(chain1$alpha.p)]
        time_1 <- chain1$lasttime
        label_1 <- chain1$lastlabel
        
        mu_2 <- chain2$mu.p[1,length(chain2$alpha.p)]
        alpha_2 <- chain2$alpha.p[1,length(chain2$alpha.p)]
        beta_2 <- chain2$beta.p[1,length(chain2$alpha.p)]
        time_2 <- chain2$lasttime
        label_2 <- chain2$lastlabel
        
        Total_itr <-10000
        burn <- 0
      }
      
    }
    
    if(gr<1.1){
      result <-list()
      for(i in 1:length(chain1)){
        result[[i]] <- cbind(chain1[[i]],chain2[[i]])
      }
      names(result) <- names(chain1)
      save(result, file = paste0(outname,"temporal_",index,".dat"))
    }
  }
}











