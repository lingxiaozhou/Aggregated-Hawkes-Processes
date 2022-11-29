

# aggregate the data by the specified bin.length
# a: a vector containing data
# bin.length: can be a value or a vector
agg <- function(a,bin.length){
  if(length(bin.length)==1){
    bin.length <- rep(bin.length,length(a))
  }
  result <- c()
  for (i in 1:length(a)) {
    tmp <- ifelse(bin.length[i]==0,a[i],a[i]%/%bin.length[i]*bin.length[i])
    result <- c(result,tmp)
  }
  
  return(result)
}




# generate time/space based on the aggregated value
#ll: vector containing the lower limit of the associated bin
#bin.length: can be a value or a vector
generateST <- function(ll,bin.length,constrain = c()){
  x <- c()
  
  if(length(bin.length)==1){
    bin.length <- rep(bin.length,length(ll))
  }
  
  for (i in 1:length(ll)){
  
    if(length(constrain)!=2){
      x <- c(x,runif(1,ll[i],ll[i]+bin.length[i]))
    }
    if(length(constrain)==2){
      x <- c(x,runif(1,max(ll[i],constrain[1]),min(ll[i]+bin.length[i],constrain[2])))
      if(is.na(runif(1,max(ll[i],constrain[1]),min(ll[i]+bin.length[i],constrain[2])))){
        print(ll[i])
        print(ll[i+bin.length])
      }
    }
    
  }
  

  return(x)
}