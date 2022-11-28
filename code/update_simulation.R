
#-----------------------------------update label-----------------------------------------

# use gibbs sampling to update labels for temporal data
updatelabel_t <- function(time,portion,mu,alpha,beta,g,type,q,use_q,updatelabel_ll){
  ind <- sample(2:length(time),ceiling(portion*length(time))-1)
  for(i in ind){
    # lower bound for time of events that we will consider
    if(use_q==TRUE){
      t.lb <- time[i]-max(qexp(q,beta[,type[i]]))
    }
    
    if(use_q==FALSE){
      t.lb <- ifelse(updatelabel_ll==0,0,time[i]-updatelabel_ll)
    }
    
    
    candidate.index <-which(time>t.lb & time<time[i])
    candidate.type <- type[candidate.index]
    
    
    Lik.tmp <-log(alpha[candidate.type,type[i]])+dexp(time[i]-time[candidate.index],rate = beta[candidate.type,type[i]],log = TRUE)
    
    
    candidate.index <- c(0,candidate.index)
    Lik <- c(log(mu[type[i]]),Lik.tmp)
    
    weight <- exp(Lik)/sum(exp(Lik))
    
    
    g[i] <- sample(candidate.index,1,prob = weight)
  }
  
  return(g)
}





# use gibbs sampling to update labels for spatio-temporal data
updatelabel_st <- function(time,portion,mu,alpha,beta,gamma,g,type,x,y,q,use_q,updatelabel_ll,w.width,w.length){
  ind <- sample(1:length(time),ceiling(portion*length(time)))
  for(i in ind){
    # lower bound for time of events that we will consider
    if(use_q==TRUE){
      t.lb <- time[i]-max(qexp(q,beta[,type[i]]))
    }
    
    if(use_q==FALSE){
      t.lb <- ifelse(updatelabel_ll==0,0,time[i]-updatelabel_ll)
    }
    
    
    
    
    candidate.index <-which(time>t.lb & time<time[i])
    candidate.type <- type[candidate.index]
    
    
    Lik.tmp <-log(alpha[candidate.type,type[i]])+dexp(time[i]-time[candidate.index],rate = beta[candidate.type,type[i]],log = TRUE)+dnorm(x[i],x[candidate.index],gamma[candidate.type,type[i]],log = TRUE)+dnorm(y[i],y[candidate.index],gamma[candidate.type,type[i]],log = TRUE)
    
    
    candidate.index <- c(0,candidate.index)
    Lik <- c(log(mu[type[i]])-log(w.width*w.length),Lik.tmp)
    
    weight <- exp(Lik)/sum(exp(Lik))
    
    
    g[i] <- sample(candidate.index,1,prob = weight)
  }
  
  return(g)
}




#---------------------------------update time-----------------------------------------------------------
# Update time for aggregated HP
updatetime <- function(time,ll,alpha,beta,g,type,bin.length,t,dim){
  # Update ith time
  for(i in 1:length(time)){
    Oi <- which(g==i)
    betai <- beta[type[i],type[Oi]]
    upper <- ifelse(is.na(Oi[1]),ll[i]+bin.length,min(ll[i]+bin.length,time[which(g==i)]))
    lower <- ifelse(g[i]==0,ll[i],max(ll[i],time[g[i]]))
    tc <- runif(1, lower,upper)
    
    
    D8 <- 0
    for (j in 1:dim) {
      D8 <-D8+alpha[type[i],j]*(pexp(t-time[i],rate = beta[type[i],j],log.p = FALSE)-pexp(t-tc,rate = beta[type[i],j],log.p = FALSE))
    }
    
    if(g[i]>0){
      D8 <- D8-dexp(time[i]-time[g[i]],rate=beta[type[g[i]],type[i]],log=TRUE)+dexp(tc-time[g[i]],rate=beta[type[g[i]],type[i]],log=TRUE)
    }
    
    if(is.na(Oi[1])==FALSE){
      D8 <- D8-sum(dexp(time[Oi]-time[i],rate=betai,log=TRUE))+sum(dexp(time[Oi]-tc,rate = betai,log = TRUE))
    }
    
    
    u <- runif(1)
    
    if(is.na(D8)){
      print(lower)
      print(upper)
      print(itr)
    }
    
    if(log(u) < D8){
      time[i] <- tc
    }
  }
  
  return(time)
}




#---------------------------------update position-------------------------------------
updatepos <- function(x,y,ll.x,ll.y,gamma,g,type,bin.length.x,bin.length.y,t){
  
  for(i in 1:length(time)){
    Oi <- which(g==i)
    pai <- g[i]
    gammai <- gamma[type[i],type[Oi]]
    
    x.c <- runif(1,ll.x[i],ll.x[i]+bin.length.x)
    y.c <- runif(1,ll.y[i],ll.y[i]+bin.length.y)
    
    
    D9 <- 0
    
    # if this event is not an immigrant 
    if(pai!=0){
      D9 <- dnorm(x.c,x[pai],gamma[type[pai],type[i]],log = TRUE)+dnorm(y.c,y[pai],gamma[type[pai],type[i]],log = TRUE)-dnorm(x[i],x[pai],gamma[type[pai],type[i]],log = TRUE)-dnorm(y[i],y[pai],gamma[type[pai],type[i]],log = TRUE)
    }
    
    # if this event has offspring
    if(length(Oi)!=0){
      D9 <- D9+sum(dnorm(x[Oi],x.c,gammai,log = TRUE))+sum(dnorm(y[Oi],y.c,gammai,log = TRUE))-sum(dnorm(x[Oi],x[i],gammai,log = TRUE))-sum(dnorm(y[Oi],y[i],gammai,log = TRUE))
    }
    
    
    u <- runif(1)
    
    if(log(u) < D9){
      x[i] <- x.c
      y[i] <- y.c
    }
  }
  
  return(list(x,y))
}
