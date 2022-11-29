# Fit aggregated multivariate Hawkes model on temporal data  
# It will aggregated the data by specified bin.length first if inputs are continuous timestamps and marks

# Input:
#   time.true: time stamps (can be exact or aggregated)
#   type: a vector indicating the corresponding process of each event (1 for the first process, 2 for the second process,...)
#   seed: set the seed if not 0
#   Total_itr: total number of iterations for MCMC
#   burn: burn-in samples
#   q: required only if use_q is TRUE
#   use_q: if TRUE, then truncated the posterior of branching labels based on the "q"-quantile of current value of bata
#   updatelabel_ll: A nonnegative integer required only  if use_q is FALSE. 
#                   If 0, the exact posterior of banching labels is used in the MCMC.
#                   If not 0, then the posterior of branching labels is truncated with possible values greater than t-updatelabel_ll,
#                   where t is the current time of the event.
#   portion: a number between 0 and 1. Percentage of time and spatial marks that are updated in an iteration
#   t: length of the observed time window
#   bin.length: aggregation size in time (0 means no aggregation)
#   mu: a vector containing initial values for mu
#   alpha: an array containing initial values for alpha
#   beta: an array containing values for beta
#   e_beta: an array containing values for the standard deviations used in the proposal distribution of beta
#   silent: if FALSE, then print the acceptance rate every "summ_itr" iterations
#   summ_itr: required only when silent is FALSE
#   save.latent: if TRUE, then save values of all latent variables (exact time, exact x, exact y, branching labels)
#   update.mu: if TRUE, then update the value of mu in the MCMC
#   update.alpha: if TRUE, then update the value of alpha in the MCMC
#   update.beta: if TRUE, then update the value of beta in the MCMC
#   update.time: if TRUE, then update the value of exact time stamps in the MCMC


#Output:
#   A dataframe containing the posterior samples of all parameters
MCMC_t <- function(time.true,label.true,type,
                    seed = 1234,Total_itr=10000,burn = 3000,
                    q=0.99998,use_q=TRUE,updatelabel_ll=0,portion=1,t=200,bin.length=0.1,
                    mu=rep(1,2),alpha=array(1,c(2,2)),beta=array(1,c(2,2)),time_0=c(),label_0=c(),
                    e_beta=array(8,c(2,2)),silent=FALSE,summ_itr=1000, save.latent = FALSE,
                    update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.time=TRUE,update.label=TRUE){
  if(seed!=0){
    set.seed(seed)
  }
  
  if(bin.length==0){
    update.time=FALSE
  }
  
  dim <- length(mu)

  
  #-------------------------initialization--------------------------------------
  if(update.time==TRUE){
    ll <- agg(time.true,bin.length)
    time <- time_0
    
    if(length(time)==0){
      time <- generateST(ll,bin.length,c(0,t))
    }
  }
  
  if(update.time==FALSE){
    time <- time.true
  }
  
  
  if(update.label==TRUE){
    g <- label_0
    if(length(label_0)==0){
      g <- initializeg3(time,10)
    }
    
  }
  if(update.label==FALSE){
    g <- label.true
    time <- time.true
  }
  
  total_count <- length(time)
  
  
  #alpha <- array(c(1,0.1,0.5,0.7),c(2,2))
  #mu <- c(1.2,0.5)
  
  
  
  #-------------------Define lists to store postburn samples----------------
  muls <- list()
  alphals <- list()
  betals <- list()
  
  if(save.latent==TRUE){
    labells <- list()
    timels <- list()
  }
  
  
  
  itr <- 1
  flag <- rep(0,dim^2)
  nameflag <- rep("",dim^2)
  for (j in 1:dim) {
    for (l in 1:dim) {
      nameflag[j+l-l%%2] <- paste0("beta",j,l)
    }
    
  }
  
  
  names(flag) <- nameflag
  
  
  
  
  
  while(itr < Total_itr){
    itr <- itr + 1
    
    for (j in 1:dim) {
      Ij <- which(g == 0 & type==j)
      if(length(Ij) < 10){
        cat("in itr",itr, paste0("number of I",j," is"),length(Ij),"\n")
      }
      
      Oj<- which(g > 0 & type==j)
      
      #------------------------ update mu---------------------------
      if(update.mu==TRUE){
        mu[j] <- rgamma(1,1+length(Ij),t+0.01)
      }
      
      
      
      for (l in 1:dim) {
        Ol <- which(g>0 & type==l)
        Ojl <- Ol[type[g[Ol]]==j]
        
        if(length(Ojl) < 10){
          cat("in itr",itr, paste0("number of O",j,l," is"),length(Ojl),"\n")
        }
        
        
        #----------------------------------update alpha---------------------------------
        if(update.alpha==TRUE){
          alpha[j,l] <- rgamma(1,1+length(Ojl),sum(1-exp(-beta[j,l]*(t-time[type==j])))+0.01)
        }
        
        
        #------------------------------------update beta-------------------------------
        if(update.beta==TRUE){
          betac <- rtruncnorm(1,mean = beta[j,l], sd = e_beta[j,l], a = 0, b = Inf)
          

          
          # use exp(0.1) as prior
          D3 <- 0.1*(beta[j,l]-betac)+sum(dexp(time[Ojl]-time[g[Ojl]],rate = betac,log = TRUE)-dexp(time[Ojl]-time[g[Ojl]],rate = beta[j,l],log = TRUE))+alpha[j,l]*sum(pexp(t-time[type==j],rate = beta[j,l],log.p = FALSE)-pexp(t-time[type==j],rate = betac,log.p = FALSE))+dtruncnorm(beta[j,l],mean = betac, sd = e_beta[j,l], a = 0, b = Inf)-dtruncnorm(betac,mean = beta[j,l], sd = e_beta[j,l], a = 0, b = Inf)
      
          u <- runif(1)
          if(log(u) < D3){
            beta[j,l] <- betac
            flag[j+l-l%%2] <- flag[j+l-l%%2]+1
          }
        }
      }
      
      
      
    }
    
    
    #-------------------------------update branching structure----------------------------------------
    if(update.label==TRUE){ 
      
      g <- updatelabel_t(time,portion,mu,alpha,beta,g,type,q,use_q,updatelabel_ll)
      
    }
    
    #-------------------------------update time--------------------------------------------------------
    if(update.time == TRUE){
      time <- updatetime(time,ll,alpha,beta,g,type,bin.length,t,dim)
    }
    
    
    
    
    
    
    
    
    
    
    
    #----------------------------save postburn samples-----------------------------------------
    if(itr > burn){
      muls[[itr - burn]] <- mu
      betals[[itr - burn]] <- beta
      alphals[[itr - burn]] <- alpha
      
      if(save.latent == TRUE){
        labells[[itr - burn]] <- g
        timels[[itr - burn]] <- time
      }
    }
    
    
    
    #------------------------------compute acceptance rate----------------------------------------
    ar <- flag / itr #Acceptance rate
    if(itr %% summ_itr==0){
      if (silent==FALSE){
        print(ar)
      }
      
    }
    
  }
  
  
  
  mu.p <- array(unlist(muls),c(dim,itr - burn))
  alpha.p <- array(unlist(alphals),c(dim^2,itr - burn))
  beta.p <- array(unlist(betals),c(dim^2,itr - burn))

  
  if(save.latent==TRUE){
    timearr <- array(unlist(timels),c(length(time),itr - burn))
    labelarr <- array(unlist(labells),c(length(time),itr - burn))
    
    return(list(mu.p=mu.p,alpha.p=alpha.p,beta.p=beta.p,labelarr=labelarr,timearr=timearr,ar=ar))
  }
  
  
  return(list(mu.p=mu.p,alpha.p=alpha.p,beta.p=beta.p,lasttime=time, lastlabel=g, ar=ar))
  
}


