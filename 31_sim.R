args <- commandArgs(TRUE)

sim <- 3

agg_size1 <- as.numeric(args[1])/100
agg_size2 <- as.numeric(args[2])/100


index <- as.numeric(args[3])


workdic <- ""     # change if needed
source(paste0(workdic,"code/simulationSource.R"))

mu <- c(0.3,0.5);alpha <- array(c(0.7,0.3,0.15,0.5),c(2,2));beta <- array(1,c(2,2)); gamma <- array(1,c(2,2))
t <- 500;w.width <- 100;w.length <- 100; 
Total_itr <- 40000; burn <- 20000; epsilon_beta <- array(c(0.15,0.3,0.3,0.3),c(2,2)) 

cat(agg_size1,"\n",agg_size2,"\n")
outname <- paste0(workdic,"Output/",sim,"_sims/1_agg",agg_size1*100,"/2_agg",agg_size2*100,"/Results/")



# if the result file already exists, then the code will not run this simulation
if(paste0("st2_",index,".dat") %in% list.files(outname)){
  cat("Exist!")
}else{


  #------------------------------------simulate data-----------------------------------------
  
  HP <- generateData(seed = index,remove = TRUE,mu = mu,alpha = alpha,beta = beta, t = t,gamma = gamma, w.width = w.width,w.length = w.length)
  
  bin.length <- rep(1,nrow(HP))
  bin.length[which(HP$type==1)] <- rep(1*agg_size1,sum(HP$type==1))
  bin.length[which(HP$type==2)] <- rep(1*agg_size2,sum(HP$type==2))
  
  gr <- 5
  
  # intial values for the MCMC
  mu_1 <- rep(1,2);alpha_1 <- array(1,c(2,2));beta_1 <- array(1,c(2,2));gamma_1 <- array(1,c(2,2));label_1 <- c();time_1 <- c();x_1 <- c();y_1 <- c()
  mu_2 <- rep(1,2);alpha_2 <- array(1,c(2,2));beta_2 <- array(1,c(2,2));gamma_2 <- array(1,c(2,2));label_2 <- c();time_2 <- c();x_2 <- c();y_2 <- c()
  
  while (gr>1.1) {
    chain1 <- MCMC_st(HP$HP,HP$x,HP$y,HP$label,HP$type,w.width,w.length,
                      seed = index,Total_itr=Total_itr,burn = burn,
                      q=0.998,use_q=TRUE,updatelabel_ll=0,portion=1,t=t,bin.length=bin.length,bin.length.x=bin.length,bin.length.y=bin.length,
                      mu=mu_1,alpha=alpha_1,beta=beta_1,gamma=gamma_1,time_0=time_1,label_0=label_1,x_0=x_1,y_0=y_1,
                      e_beta=epsilon_beta,silent=FALSE,summ_itr=10000, save.latent = FALSE,
                      update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.gamma=TRUE,update.time=TRUE,update.pos=TRUE,update.label=TRUE)
    
    chain2 <- MCMC_st(HP$HP,HP$x,HP$y,HP$label,HP$type,w.width,w.length,
                      seed = index+1000,Total_itr=Total_itr,burn = burn,
                      q=0.998,use_q=TRUE,updatelabel_ll=0,portion=1,t=t,bin.length=bin.length,bin.length.x=bin.length,bin.length.y=bin.length,
                      mu=mu_2,alpha=alpha_2,beta=beta_2,gamma=gamma_2,time_0=time_2,label_0=label_2,x_0=x_2,y_0=y_2,
                      e_beta=epsilon_beta,silent=FALSE,summ_itr=10000, save.latent = FALSE,
                      update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.gamma=TRUE,update.time=TRUE,update.pos=TRUE,update.label=TRUE)
    
    #------------------------------check convergence of MCMC-----------------------------------------
    
    combined <- combineMCMC(chain1,chain2,model = "st2")
    tmp <- gelman.diag(combined,autoburnin = FALSE,multivariate = TRUE)$psrf[,1]
    gr <- max(tmp)
    
    if(gr>1.1){
      cat("parameter with convergence problem:", which(tmp>1.1),"\n")
      if(any(chain1$ar>0.5) | any(chain1$ar < 0.15)){
        epsilon_tmp <- c(epsilon_beta)
        for (i in 1:4) {
          epsilon_tmp[i] <- ifelse(chain1$ar[i]>0.5,epsilon_tmp[i]*2,epsilon_tmp[i])
          epsilon_tmp[i] <- ifelse(chain1$ar[i]<0.15,epsilon_tmp[i]/2,epsilon_tmp[i])
        }
        epsilon_beta <- array(epsilon_tmp,c(2,2))
        
        mu_1 <- rep(1,2);alpha_1 <- array(1,c(2,2));beta_1 <- array(1,c(2,2));gamma_1 <- array(1,c(2,2));label_1 <- c();time_1 <- c();x_1 <- c();y_1 <- c()
        mu_2 <- rep(1,2);alpha_2 <- array(1,c(2,2));beta_2 <- array(1,c(2,2));gamma_2 <- array(1,c(2,2));label_2 <- c();time_2 <- c();x_2 <- c();y_2 <- c()
        
        Total_itr <-40000
        burn <- 20000
      }else{
        mu_1 <- chain1$mu.p[,length(chain1$alpha.p[1,])]
        alpha_1 <- array(chain1$alpha.p[,length(chain1$alpha.p[1,])],c(2,2))
        beta_1 <- array(chain1$beta.p[,length(chain1$alpha.p[1,])],c(2,2))
        gamma_1 <- array(chain1$gamma.p[,length(chain1$alpha.p[1,])],c(2,2))
        time_1 <- chain1$lasttime
        label_1 <- chain1$lastlabel
        x_1 <- chain1$lastx
        y_1 <- chain1$lasty
        
        
        mu_2 <- chain2$mu.p[,length(chain2$alpha.p[1,])]
        alpha_2 <- array(chain2$alpha.p[,length(chain2$alpha.p[1,])],c(2,2))
        beta_2 <- array(chain2$beta.p[,length(chain2$alpha.p[1,])],c(2,2))
        gamma_2 <- array(chain2$gamma.p[,length(chain1$alpha.p[1,])],c(2,2))
        time_2 <- chain2$lasttime
        label_2 <- chain2$lastlabel
        x_2 <- chain2$lastx
        y_2 <- chain2$lasty
        
        Total_itr <-20000
        burn <- 0
      }
    }
    
    if(gr<1.1){
      result <-list()
      for(i in 1:length(chain1)){
        result[[i]] <- cbind(chain1[[i]],chain2[[i]])
      }
      names(result) <- names(chain1)
      save(result, file = paste0(outname,"st2_",index,".dat"))
    }
  }
}