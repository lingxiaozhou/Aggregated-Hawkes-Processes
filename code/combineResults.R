# --------------------read result files for temporal simulations--------
load("~/Github/Aggregated-Hawkes-Processes/result/Simulation/MCEM.dat")

for (agg in c(0,50,75,100,125,150)) {
  tmp <- read.sum(truevalue = c(0.3,0.7,1),start=1,end=400,prefix = paste0("~/Simulation/Output/1_sims/mu03/alpha07/beta1/agg",agg,"/Results/temporal_"))
  assign(paste0("t_",agg,"37"),tmp)
}

for (agg in c(0,50,75,100,125,150)) {
  tmp <- read.sum(truevalue = c(0.5,0.5,1),start=1,end=400,prefix = paste0("~/Simulation/Output/1_sims/mu05/alpha05/beta1/agg",agg,"/Results/temporal_"))
  assign(paste0("t_",agg,"55"),tmp)
}

# make a dataframe for mu=0.3
t_df1 <- data.frame(matrix(ncol = 6, nrow = 0))
for (i in c(0,50,75,100,125,150)) {
  
  cat(i,"\n")
  tmp <- get(paste0("t_",i,"37"))
  tmp.mcem <- get(paste0("arr_",i))
  tmpdf <-data.frame(apply(tmp,2,"c"))  
  colnames(tmpdf) <- c("mean","median","CI.length","Coverage")
  tmpdf$parameter <- c("mu","alpha","beta")
  tmpdf$agg <- rep(i,nrow(tmp))
  tmpdf$mu <- "0.3"
  tmpdf$MLE <- c(t(tmp.mcem))
  tmpdf$truevalue <- c(0.3,0.7,1)
  tmpdf$bias <- tmpdf$mean-tmpdf$truevalue
  
  t_df1 <- rbind.data.frame(t_df1,tmpdf)
}
t_df1$agg <- t_df1$agg/100
t_df1$agg <- factor(t_df1$agg,c("0","0.5","0.75","1","1.25","1.5"))
t_df1$parameter <- factor(t_df1$parameter,levels = c("mu","alpha","beta"))


# make a dataframe for mu=0.5
t_df2 <- data.frame(matrix(ncol = 6, nrow = 0))
for (i in c(0,50,75,100,125,150)) {
  
  cat(i,"\n")
  tmp <- get(paste0("t_",i,"55"))
  tmp.mcem <- get(paste0("arr_",i,"_55"))
  tmpdf <-data.frame(apply(tmp,2,"c"))  
  colnames(tmpdf) <- c("mean","median","CI.length","Coverage")
  tmpdf$parameter <- c("mu","alpha","beta")
  tmpdf$agg <- rep(i,nrow(tmp))
  tmpdf$mu <- "0.5"
  tmpdf$MLE <- c(t(tmp.mcem))
  tmpdf$truevalue <- c(0.3,0.7,1)
  tmpdf$bias <- tmpdf$mean-tmpdf$truevalue
  
  t_df2 <- rbind.data.frame(t_df2,tmpdf)
}
t_df2$agg <- t_df2$agg/100
t_df2$agg <- factor(t_df2$agg,c("0","0.5","0.75","1","1.25","1.5"))
t_df2$parameter <- factor(t_df2$parameter,levels = c("mu","alpha","beta"))


# make a dataframe for all temporal simulation
t_df <- rbind.data.frame(t_df1,t_df2)
# save the dataframe
save(t_df, file = "~/Github/Aggregated-Hawkes-Processes/result/Simulation/temporal.dat")



# ------------------read result files for single spatio-temporal simulations--------

# for parameter set 1(mu = 0.3) 
for (tagg in c(0,50,100,150)) {
  for (sagg in c(0,50,100,150)) {
    print(tagg)
    print(sagg)
    tmp <- read.sum(model = "st1", truevalue = c(0.3,0.7,1,1),start=1,end=400,prefix = paste0("~/Simulation/Output/2_sims/mu03/alpha07/beta1/tagg",tagg,"/sagg",sagg,"/Results/st_"))
    assign(paste0("st_",tagg,sagg),tmp)
  }
}


st_df1 <- data.frame(matrix(ncol = 6, nrow = 0))
for (i in c(0,50,100,150)) {
  for (j in c(0,50,100,150)) {
    cat(i,j,"\n")
    tmp <- get(paste0("st_",i,j))
    
    tmpdf <-data.frame(apply(tmp,2,"c"))  
    colnames(tmpdf) <- c("mean","median","CI.length","Coverage")
    tmpdf$parameter <- c("mu","alpha", "beta","gamma")
    tmpdf$tagg <- rep(i,nrow(tmp))
    tmpdf$sagg <- rep(j,nrow(tmp))
    tmpdf$truevalue <- c(0.3,0.7,1,1)
    tmpdf$bias <- tmpdf$mean-tmpdf$truevalue
    
    st_df1 <- rbind.data.frame(st_df1,tmpdf)
  }
  
  
}

st_df1$parameter <- factor(st_df1$parameter,levels = c("mu","alpha", "beta","gamma"))
st_df1$tagg <- factor(st_df1$tagg/100, levels = c("0","0.5","1","1.5"))
st_df1$sagg <- factor(st_df1$sagg/100, levels = c("0","0.5","1","1.5"))
st_df1$par_set <- "1"

# for parameter set 2(mu = 0.5) 
for (tagg in c(0,50,100,150)) {
  for (sagg in c(0,50,100,150)) {
    print(tagg)
    print(sagg)
    tmp <- read.sum(model = "st1", truevalue = c(0.5,0.5,1,1),start=1,end=400,prefix = paste0("~/Simulation/Output/2_sims/mu05/alpha05/beta1/tagg",tagg,"/sagg",sagg,"/Results/st_"))
    assign(paste0("st_",tagg,sagg,"_55"),tmp)
  }
}


st_df2 <- data.frame(matrix(ncol = 6, nrow = 0))
for (i in c(0,50,100,150)) {
  for (j in c(0,50,100,150)) {
    cat(i,j,"\n")
    tmp <- get(paste0("st_",i,j,"_55"))
    
    tmpdf <-data.frame(apply(tmp,2,"c"))  
    colnames(tmpdf) <- c("mean","median","CI.length","Coverage")
    tmpdf$parameter <- c("mu","alpha", "beta","gamma")
    tmpdf$tagg <- rep(i,nrow(tmp))
    tmpdf$sagg <- rep(j,nrow(tmp))
    tmpdf$truevalue <- c(0.5,0.5,1,1)
    tmpdf$bias <- tmpdf$mean-tmpdf$truevalue
    
    st_df2 <- rbind.data.frame(st_df2,tmpdf)
  }
  
  
}

st_df2$parameter <- factor(st_df2$parameter,levels = c("mu","alpha", "beta","gamma"))
st_df2$tagg <- factor(st_df2$tagg/100, levels = c("0","0.5","1","1.5"))
st_df2$sagg <- factor(st_df2$sagg/100, levels = c("0","0.5","1","1.5"))
st_df2$par_set <- "2"

st_df <- rbind.data.frame(st_df1,st_df2)
save(st_df, file = "~/Github/Aggregated-Hawkes-Processes/result/Simulation/spatio_temporal.dat")





#--------------------read result files for multiple spatio-temporal simulations----------------------------------
for (agg1 in c(0,75,100)) {
  for (agg2 in c(0,75,100)) {
    print(agg1)
    print(agg2)
    tmp <- read.sum(model = "st2", truevalue = c(0.3,0.5,0.7,0.3,0.15,0.5,1,1,1,1,1,1,1,1),start=1,end=400,prefix = paste0("~/Simulation/Output/3_sims/1_agg",agg1,"/2_agg",agg2,"/Results/st2_"))
    assign(paste0("st2_",agg1,agg2),tmp)
  }
}



# make a dataframe st2
st2_df <- data.frame(matrix(ncol = 6, nrow = 0))
for (i in c(0,75,100)) {
  for (j in c(0,75,100)) {
    cat(i,j,"\n")
    tmp <- get(paste0("st2_",i,j))
    
    tmpdf <-data.frame(apply(tmp,2,"c"))  
    colnames(tmpdf) <- c("mean","median","CI.length","Coverage")
    tmpdf$parameter <- c("mu[1]","mu[2]","alpha[11]","alpha[21]","alpha[12]","alpha[22]","beta[11]","beta[21]","beta[12]","beta[22]","gamma[11]","gamma[21]","gamma[12]","gamma[22]")
    tmpdf$agg1 <- rep(i,nrow(tmp))
    tmpdf$agg2 <- rep(j,nrow(tmp))
    tmpdf$truevalue <- c(0.3,0.5,0.7,0.3,0.15,0.5,1,1,1,1,1,1,1,1)
    tmpdf$bias <- tmpdf$mean - tmpdf$truevalue
    
    st2_df <- rbind.data.frame(st2_df,tmpdf)
  }
  
  
}
st2_df$parameter <- factor(st2_df$parameter,levels = c("mu[1]","mu[2]","alpha[11]","alpha[21]","alpha[12]","alpha[22]","beta[11]","beta[21]","beta[12]","beta[22]","gamma[11]","gamma[21]","gamma[12]","gamma[22]"))
st2_df$agg1 <- factor(st2_df$agg1/100, levels = c("0","0.75","1"))
st2_df$agg2 <- factor(st2_df$agg2/100, levels = c("0","0.75","1"))

save(st2_df, file = "~/Github/Aggregated-Hawkes-Processes/result/Simulation/multi_spatio_temporal.dat")




