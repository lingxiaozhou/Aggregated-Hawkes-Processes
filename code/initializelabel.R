# functions that initialize the branching structure
# time is a sorted array
# If i is the index of current point, 0, i-3, i-2,i-1 will be the possible labels
initializeg1 <- function(time){
  g <- rep(0,length(time))
  g[2] <- sample(c(0,1),1)
  g[3] <- sample(c(0,1,2),1)
  for (i in 4:length(g)){
    g[i] <- sample(c(0,i-1,i-2,i-3),1)
  }
  
  return(g)
}



# functions that initialize the branching structure
# time is a sorted array
# If i is the index of current point, 0, i-r,.., i-2,i-1 will be the possible labels
initializeg2 <- function(time,r){
  g <- rep(0,length(time))
  
  for (i in 2:length(g)){
    lb <- max(0,i-r)
    g[i] <- sample(c(0,lb:(i-1)),1)
  }
  
  return(g)
}


# functions that initialize the branching structure
# time need not be sorted
# If i is the index of current point, 0, i-r,.., i-2,i-1 will be the possible labels
initializeg3 <- function(time,r){
  g <- rep(0,length(time))
  ord <- rank(time)
  
  for (i in 1:length(time)){
    if(ord[i]==1){
      g[i] <- 0
    }
    if(ord[i]>1){
      lb <- max(0,ord[i]-r)
      g[i] <- sample(c(0,which(ord%in%lb:(ord[i]-1))),1)
    }
    
  }
  
  return(g)
}
