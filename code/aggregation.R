# aggregate the data by the specified bin.length
# a: a vector containing data
# bin.length: can be a value or a vector
agg <- function(a,bin.length){
  return(a%/%bin.length*bin.length)
}




# generate time/space based on the aggregated value
#ll: vector containing the lower limit of the associated bin
#bin.length: can be a value or a vector
generateST <- function(ll,bin.length){
  x <- c()
  if(length(bin.length)==1){
    for (i in 1:length(ll)){
      x <- c(x,runif(1,ll[i],ll[i]+bin.length))
    }
  }
  
  if(length(bin.length)>1){
    for (i in 1:length(ll)){
      x <- c(x,runif(1,ll[i],ll[i]+bin.length[i]))
    }
  }
  
  
  return(x)
}