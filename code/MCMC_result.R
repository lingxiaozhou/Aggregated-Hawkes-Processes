# compute the mean, median, upper and lower limits of 95% credible interval 
# and coverage based on the posterior samples of a parameter

# input:
# x: posterior samples
# x1: true value of the paramter

computeRes <- function(x, x1){
  ub <- quantile(x,0.975)
  lb <- quantile(x,0.025)
  
  return(list(mean(x),median(x), ub-lb,ub>= x1&& lb<=x1))
  
}





pos.sum <- function(ls){
  tmp <- array(0,c(0,4))
  name.tmp1 <- names(ls)
  name.tmp2 <- c(11,21,12,22)
  name.r <- c()
  for (i in 1:length(ls) ) {
    var <- ls[[i]]
    for (j in 1:nrow(var)) {
      m <- mean(var[j,])
      ub <- quantile(var[j,],0.975)
      lb <- quantile(var[j,],0.025)
      l <- ub-lb
      tmp <- rbind(tmp,c(m,lb,ub,l))
      name.r <- c(name.r,paste0(name.tmp1[i],name.tmp2[j]))
    }
    
  }
  colnames(tmp) <- c("pos.mean","CI.L","CI.U","CI.length")
  rownames(tmp) <- name.r
  return(tmp)
}




# Read and combine the results of all simulations
# input:
# truevalue: a vector containing the true value of all parameters
# start: the first simulation index
# end: the last simulation index
# prefix: the full name of the result file without the substring indicating the index of the simulation
# chain.name: the name of the saved dataframe
# model: possible values are:
#   "t": temporal process
#   "st": single spatio-temporal process
#   "st2": two mutually exciting spatio-temporal processes

read.sum <- function(truevalue, start = 1, end = 400,prefix="~/Simulation/Output/1_sims/mu03/alpha07/beta1/agg0/Results/temporal_",chain.name="result",model = "t"){
  nsim <- end-start+1
  if(model=="t"){
    npar <- 3
  }
  if(model == "st1"){
    npar <- 4
  }
  if(model=="st2"){
    npar <- 14
  }
  
  s <- array(NA,c(npar,4,nsim))
  valid <- c()
  for(i in start:end){
    tryCatch({
      load(paste0(prefix, i, '.dat'))
      valid <- c(valid,i-(start-1))
      df <- get(chain.name)
      if(model=="t"){
        datadf <- t(rbind(df$mu.p[1,],df$alpha.p[1,], df$beta.p[1,]))
      }
      if(model=="st1"){
        datadf <- t(rbind(df$mu.p[1,],df$alpha.p[1,], df$beta.p[1,],df$gamma.p[1,]))
      }
      if(model=="st2"){
        datadf <- t(rbind(df$mu.p[1:2,],df$alpha.p[1:4,], df$beta.p[1:4,],df$gamma.p[1:4,]))
      }
      
      
      tmp <- array(0,dim = c(0,4))
      for (j in 1:length(truevalue)) {
        tmp <- rbind(tmp, unlist(computeRes(datadf[,j],truevalue[j])))
      }
      
      s[1:npar,,i-(start-1)] <- tmp
      
      
    }, error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
    })
  }
  
  #s <- s[,,valid]
  
  return(s)
}