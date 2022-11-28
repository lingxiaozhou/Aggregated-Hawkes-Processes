combineMCMC <- function(df1,df2,model = "t1"){
  if(model=="t1"){
    c1 <- mcmc(cbind(df1$mu.p[1,],df1$alpha.p[1,],df1$beta.p[1,]))
    c2 <- mcmc(cbind(df2$mu.p[1,],df2$alpha.p[1,],df2$beta.p[1,]))
  }
  if(model=="t2"){
    tmp <- cbind(df1$mu.p[1,],df1$mu.p[2,])
    for (i in 1:4) {
      tmp <- cbind(tmp,df1$alpha.p[i,],df1$beta.p[i,])
    }
    c1 <- mcmc(tmp)
    
    tmp <- cbind(df2$mu.p[1,],df2$mu.p[2,])
    for (i in 1:4) {
      tmp <- cbind(tmp,df2$alpha.p[i,],df2$beta.p[i,])
    }
    c2 <- mcmc(tmp)
  }
  
  if(model=="st1"){
    c1 <- mcmc(cbind(df1$mu.p[1,],df1$alpha.p[1,],df1$beta.p[1,],df1$gamma.p[1,]))
    c2 <- mcmc(cbind(df2$mu.p[1,],df2$alpha.p[1,],df2$beta.p[1,],df2$gamma.p[1,]))
  }
  
  
  combined <- mcmc.list(c1,c2)
  
  return(combined)
}