# Fit aggregated multivariate Hawkes model with spatial marks.  
# It will aggregated the data by specified bin.length first if inputs are continuous timestamps and marks.

# Input:
#   time.true: time stamps (can be exact or aggregated)
#   x.true: x coodinates (can be exact or aggregated)
#   y.true: y coodinates (can be exact or aggregated)
#   type: a vector indicating the corresponding process of each event (1 for the first process, 2 for the second process,...)
#   w.width: width of the observed window
#   w.length: length of the observed window
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
#   bin.length.x: aggregation size in x coordinates (0 means no aggregation)
#   bin.length.y: aggregation size in y coordinates (0 means no aggregation)
#   mu: a vector containing initial values for mu
#   alpha: an array containing initial values for alpha
#   beta: an array containing values for beta
#   gamma: an array containing values for gamma
#   e_beta: an array containing values for the standard deviations used in the proposal distribution of beta
#   silent: if FALSE, then print the acceptance rate every "summ_itr" iterations
#   summ_itr: required only when silent is FALSE
#   save.latent: if TRUE, then save values of all latent variables (exact time, exact x, exact y, branching labels)
#   update.mu: if TRUE, then update the value of mu in the MCMC
#   update.alpha: if TRUE, then update the value of alpha in the MCMC
#   update.beta: if TRUE, then update the value of beta in the MCMC
#   update.gamma: if TRUE, then update the value of gamma in the MCMC
#   update.time: if TRUE, then update the value of exact time stamps in the MCMC
#   update.pos: if TRUE, then update the values of exact x and exact y in the MCMC
#   update.label: if TRUE, then update the values of branching labels in the MCMC

#Output:
#   A dataframe containing the posterior samples of all parameters.
                    
MCMC_st <- function(time.true,x.true,y.true,label.true,type,w.width,w.length,
                    seed = 1234,Total_itr=10000,burn = 3000,
                    q=0.9998,use_q=TRUE,updatelabel_ll=0,portion=1,t=200,bin.length=0.1,bin.length.x=0.1,bin.length.y=0.1,
                    mu=rep(1,2),alpha=array(1,c(2,2)),beta=array(1,c(2,2)),gamma=array(1,c(2,2)),
                    e_beta=array(8,c(2,2)),silent=FALSE,summ_itr=1000, save.latent = FALSE,
                    update.beta=TRUE,update.alpha=TRUE,update.mu=TRUE,update.gamma=TRUE,update.time=TRUE,update.pos=TRUE,update.label=TRUE){
  if(seed!=0){
    set.seed(seed)
  }
  
  if(bin.length == 0){
    update.time=FALSE
  }
  if(bin.length.x == 0 | bin.length.y==0){
    update.pos=FALSE
  }
  
  dim <- length(mu)

  
  #-------------------------initialization--------------------------------------
  if(update.time==TRUE){
    ll <- agg(time.true,bin.length)
    time <- generateST(ll,bin.length)
  }
  if(update.time==FALSE){
    time <- time.true
  }
  
  

  if(update.pos==TRUE){
    ll.x <- agg(x.true,bin.length.x)
    ll.y <- agg(y.true,bin.length.y)
    x <- generateST(ll.x,bin.length.x)
    y <- generateST(ll.y,bin.length.y)

    # x <- x[order(time)]
    # y <- y[order(time)]
  }
  if(update.pos==FALSE){
    x <- x.true
    y <- y.true
  }
  
  
  # time <- time[order(time)]
  
  if(update.label==TRUE){
    g <- initializeg3(time,10)
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
  gammals <- list()
  
  if(save.latent==TRUE){
    labells <- list()
    xls <- list()
    yls <- list()
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
        
        
        
        #----------------------------------update gamma----------------------------------
        if(update.gamma == TRUE){
          
          Ol <- which(g>0 & type==l)
          Ojl <- Ol[type[g[Ol]]==j]
          
          gamma[j,l] <- sqrt(1/rgamma(1,0.001+length(Ojl),0.001+1/2*(sum((x[Ojl]-x[g[Ojl]])^2+(y[Ojl]-y[g[Ojl]])^2))))
          
        }
        
        
      }
      
      
      
    }
    
    
    #-------------------------------update branching structure----------------------------------------
    if(update.label==TRUE){ 
      
      g <- updatelabel_st(time,portion,mu,alpha,beta,gamma,g,type,x,y,q,use_q,updatelabel_ll,w.width,w.length)
      
    }
    
    #-------------------------------update time--------------------------------------------------------
    if(update.time == TRUE){
      time <- updatetime(time,ll,alpha,beta,g,type,bin.length,t,dim)
    }
    
    
    
    #-------------------------------update position---------------------------------------------------
    if(update.pos == TRUE){
      tmp <- updatepos(x,y,ll.x,ll.y,gamma,g,type,bin.length.x,bin.length.y,t)
      x <- tmp[[1]]
      y <- tmp[[2]]
    }
    
    
    
    
    

    
    

    #----------------------------save postburn samples-----------------------------------------
    if(itr > burn){
      muls[[itr - burn]] <- mu
      betals[[itr - burn]] <- beta
      alphals[[itr - burn]] <- alpha
      gammals[[itr - burn]] <- gamma
      
      if(save.latent == TRUE){
        labells[[itr - burn]] <- g
        xls[[itr - burn]] <- x
        yls[[itr - burn]] <- y
        timels[[itr - burn]] <- time
      }
    }
    
    
    
    #------------------------------compute acceptance rate----------------------------------------
    if(itr %% summ_itr==0){
      ar <- flag / itr #Acceptance rate
      if (silent==FALSE){
        print(ar)
      }
      
    }
    
  }
  
  
  #--------------------------------save postburn samples in all iterations-------------------------------------
  mu.p <- array(unlist(muls),c(dim,itr - burn))
  alpha.p <- array(unlist(alphals),c(dim^2,itr - burn))
  beta.p <- array(unlist(betals),c(dim^2,itr - burn))
  gamma.p <- array(unlist(gammals),c(dim^2,itr - burn))
  
  if(save.latent==TRUE){
    xarr <- array(unlist(xls),c(length(time),itr - burn))
    yarr <- array(unlist(yls),c(length(time),itr - burn))
    timearr <- array(unlist(timels),c(length(time),itr - burn))
    labelarr <- array(unlist(labells),c(length(time),itr - burn))
    
    return(list(mu.p=mu.p,alpha.p=alpha.p,beta.p=beta.p,gamma.p=gamma.p,labelarr=labelarr,x=xarr,y=yarr,timearr=timearr))
  }

  
  return(list(mu.p=mu.p,alpha.p=alpha.p,beta.p=beta.p,gamma.p=gamma.p))
  
}


