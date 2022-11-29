# generate a spatial temporal Hawkes process on [0,t]x[0,w.length]x[0,w.width]
# input:
#   mu,alpha,beta,gamma,t,w.length,w.width: model parameters
#   remove: if True, remove events outside the window
generateData <- function(mu=c(1,1),alpha=array(c(0.5,0.1,0.2,0.5),dim=c(2,2)),beta=array(50,dim=c(2,2)),t = 100,gamma=array(0.3,dim=c(2,2)),seed = 0,w.width = 10, w.length = 10, remove = TRUE){
  
  if(seed!=0){
    set.seed(seed)
  }
  
  HP <- c()
  dim <- length(mu)
  k <- rpois(dim,mu*t)  # number of immigrants
  C <- list()
  
  tmp <- array(numeric(),c(0,5)) 
  for(j in 1:dim){
    if(k[j]>0){
      C[[j]] <- runif(k[j],min = 0,max = t)  #immigrants
      label <- rep(0,length(C[[j]]))
      HP <- C[[j]]
      
      x <- runif(length(C[[j]]),0,w.length)
      y <- runif(length(C[[j]]),0,w.width)
      type <- rep(j,k[j])
      tmp <- rbind(tmp,cbind(HP,label,x,y,type))
    }
    df <- as.data.frame(tmp)
  }  
  while (min(k)>0){  
    temp <- cbind.data.frame(HP = c(),label = c(),x = c(),y = c(),type=c())
    for (j in 1:dim) {
      
      
      if(k[j]>0){
        for (l in 1:dim){
          D <- rpois(k[j],alpha[j,l]) # number of type l children for each parent
          for (i in 1:k[j]){
            if(D[i]>0){
              E <- rexp(D[i],beta[j,l])
              HP <- C[[j]][i]+E
              label <- rep(C[[j]][i],length(HP))
              x <- rnorm(D[i], mean = df[df$HP == C[[j]][i],]$x, sd = gamma[j,l])
              y <- rnorm(D[i], mean = df[df$HP == C[[j]][i],]$y, sd = gamma[j,l])
              type <- rep(l,length(HP))
              temp <- rbind(temp,cbind(HP,label,x,y,type))
            }
          }
        }
      }
      
      
      
      
      temp <- temp[temp$HP<t,]  # remove descendants outside [0, T]
      

      
    }
    df <- rbind(df,temp)
    for (j in 1:dim) {
      k[j] <- sum(temp$type==j)
      C[[j]] <- temp[temp$type==j,]$HP
    }
  } 
  
  
  df <- df[order(df$HP),]
  HP <- df$HP

  x <- df$x
  y <- df$y
  
  

  
  if(remove == TRUE){
    df <-df[df$x<w.length & df$x>0 & df$y<w.width & df$y>0,]
  }
  
  label <- df$label
  for (i in 1: length(label)){
    if(label[i]>0){
      tmp <- which(round(df$HP,12)==round(label[i],12))
      label[i] <- ifelse(length(tmp)==1,tmp,NA)
    }
  }
  df$label <- label
  
  return(df)
}
