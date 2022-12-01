args <- commandArgs(TRUE)

sim <- 2

tagg_size <- as.numeric(args[1])/100
sagg_size <- as.numeric(args[2])/100
mu <- as.numeric(args[3])/10
alpha <- as.numeric(args[4])/10

index <- as.numeric(args[5])


# 

workdic <- ""   # change this if needed
source(paste0(workdic,"code/simulationSource.R"))

beta <- 1; gamma <- 1
t <- 500;w.width <- 100;w.length <- 100; 
Total_itr <- 40000; burn <- 20000; epsilon_beta <- 0.3
bin.length <- 1; bin.length.x <- 1; bin.length.y <- 1
cat(mu,"\n",alpha,"\n",beta,"\n",tagg_size,"\n",sagg_size,"\n")
outname <- paste0(workdic,"Output/",sim,"_sims/mu0",mu*10,"/alpha0",alpha*10,"/beta1/tagg",tagg_size*100,"/sagg",sagg_size*100,"/Results/")


# if the result file already exists, then the code will not run this simulation
if(paste0("st_",index,".dat") %in% list.files(outname)){
  cat("Exist!")
}else{


  #------------------------------------simulate data-----------------------------------------
  
  HP <- generateData(seed = index,remove = TRUE,mu = c(mu),alpha = array(alpha,c(1,1)),beta = array(beta,c(1,1)), t = t,gamma = array(gamma,c(1,1)), w.width = w.width,w.length = w.length)
  
  
  gr <- 5
  
  # intial values for the MCMC
  mu_1 <- 1;alpha_1 <- 1;beta_1 <- 1;label_1 <- c();time_1 <- c();x_1 <- c();y_1 <- c()
  mu_2 <- 1;alpha_2 <- 1;beta_2 <- 1;label_2 <- c();time_2 <- c();x_2 <- c();y_2 <- c()
  
  while (gr>1.1) {
    chain1 <- MCMC_st(HP$HP,HP$x,HP$y,HP$label,HP$type,w.width,w.length,
                      seed = index,Total_itr=Total_itr,burn = burn,
                      q=0.998,use_q=TRUE,updatelabel_ll=0,portion=1,t=t,bin.length=tagg_size*bin.length,bin.length.x=sagg_size*bin.length.x,bin.length.y=sagg_size*bin.length.y,
                      mu=rep(1,1),alpha=array(1,c(1,1)),beta=array(1,c(1,1)),gamma=array(1,c(1,1)),time_0=time_1,label_0=label_1,x_0=x_1,y_0=y_1,
                      e_beta=array(epsilon_beta,c(1,1)),silent=FALSE,summ_itr=10000, save.latent = FALSE,
                      update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.gamma=TRUE,update.time=TRUE,update.pos=TRUE,update.label=TRUE)
    
    chain2 <- MCMC_st(HP$HP,HP$x,HP$y,HP$label,HP$type,w.width,w.length,
                      seed = index+200,Total_itr=Total_itr,burn = burn,
                      q=0.998,use_q=TRUE,updatelabel_ll=0,portion=1,t=t,bin.length=tagg_size*bin.length,bin.length.x=sagg_size*bin.length.x,bin.length.y=sagg_size*bin.length.y,
                      mu=rep(1,1),alpha=array(1,c(1,1)),beta=array(1,c(1,1)),gamma=array(1,c(1,1)),time_0=time_2,label_0=label_2,x_0=x_2,y_0=y_2,
                      e_beta=array(epsilon_beta,c(1,1)),silent=FALSE,summ_itr=10000, save.latent = FALSE,
                      update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.gamma=TRUE,update.time=TRUE,update.pos=TRUE,update.label=TRUE)
    
    #------------------------------check convergence of MCMC-----------------------------------------
    
    combined <- combineMCMC(chain1,chain2,model = "st1")
    tmp <- gelman.diag(combined,autoburnin = FALSE,multivariate = TRUE)$psrf[,1]
    gr <- max(tmp)
    
    if(gr>1.1){
      cat("parameter with convergence problem:", which(tmp>1.1),"\n")
      if(chain1$ar>0.5 | chain1$ar < 0.15){
        epsilon_beta <- ifelse(chain1$ar>0.5,epsilon_beta*2,epsilon_beta/2)
        mu_1 <- 1;alpha_1 <- 1;beta_1 <- 1;label_1 <- c();time_1 <- c()
        mu_2 <- 1;alpha_2 <- 1;beta_2 <- 1;label_2 <- c();time_2 <- c()
        Total_itr <-40000
        burn <- 20000
      }else{
        mu_1 <- chain1$mu.p[1,length(chain1$alpha.p)]
        alpha_1 <- chain1$alpha.p[1,length(chain1$alpha.p)]
        beta_1 <- chain1$beta.p[1,length(chain1$alpha.p)]
        time_1 <- chain1$lasttime
        label_1 <- chain1$lastlabel
        x_1 <- chain1$lastx
        y_1 <- chain1$lasty
        
        
        mu_2 <- chain2$mu.p[1,length(chain2$alpha.p)]
        alpha_2 <- chain2$alpha.p[1,length(chain2$alpha.p)]
        beta_2 <- chain2$beta.p[1,length(chain2$alpha.p)]
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
      save(result, file = paste0(outname,"st_",index,".dat"))
    }
  }
}