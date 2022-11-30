library(ggplot2)

# workdic <- "~/Github/Aggregated-Hawkes-Processes/"
workdic <- ""   #change if needed

load(paste0(workdic,"result/Simulation/temporal.dat"))
load(paste0(workdic,"result/Simulation/spatio_temporal.dat"))
load(paste0(workdic,"result/Simulation/multi_spatio_temporal.dat"))

#--------------Figure A1(a)---------------------------------

tmp1 <- t_df[,c(1,5,6,7)]
colnames(tmp1) <- c("estimate","parameter","agg","mu")
tmp1$method <- "Bayesian"
tmp2 <- t_df[,c(8,5,6,7)]
colnames(tmp2) <- c("estimate","parameter","agg","mu")
problemindex <- which(tmp2$parameter=="beta" & tmp2$estimate>50)
tmp2[problemindex,1] <- NA
tmp2[problemindex-1,1] <- NA
tmp2[problemindex-2,1] <- NA
tmp2$method <- "MCEM"
tmp3 <- rbind.data.frame(tmp1,tmp2)
compare_df <- tmp3


truevalue1 <- data.frame(c(0.3,0.7,1),c("mu","alpha","beta"))
colnames(truevalue1) <- c("truevalue","parameter")
truevalue1$parameter <- factor(truevalue1$parameter,levels = c("mu","alpha","beta"))

truevalue2 <- data.frame(c(0.5,0.5,1),c("mu","alpha","beta"))
colnames(truevalue2) <- c("truevalue","parameter")
truevalue2$parameter <- factor(truevalue2$parameter,levels = c("mu","alpha","beta"))


ggplot(compare_df[compare_df$mu=="0.5",], aes(x=agg,y=estimate)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=0.5, notch=FALSE,lwd = 0.3)+
  facet_grid(parameter~method,scales = "free_y",labeller = label_parsed)+
  theme_classic() + geom_hline(data = truevalue2, aes(yintercept=truevalue), linetype="dashed", color = "red") +
  theme(axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),strip.text.y = element_text(angle = 0,size=15,margin = margin(0,0.3,0,0.3, "cm")))+
  xlab("Aggregation size in time")+
  ylab("Estimate")



#--------------------Figure 2--------------------------------------

problemindex <- which(tmp2$parameter=="beta" & tmp2$estimate>50)
tmp2[problemindex,1] <- NA
tmp2[problemindex-1,1] <- NA
tmp2[problemindex-2,1] <- NA
tmp2$method <- "MCEM"
tmp3 <- rbind.data.frame(tmp1,tmp2)
compare_df <- tmp3

ggplot(compare_df[compare_df$mu=="0.3",], aes(x=agg,y=estimate)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=0.5, notch=FALSE,lwd = 0.3)+
  facet_grid(parameter~method,scales = "free_y",labeller = label_parsed)+
  theme_classic() + geom_hline(data = truevalue1, aes(yintercept=truevalue), linetype="dashed", color = "red") +
  theme(axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),strip.text.y = element_text(angle = 0,size=15,margin = margin(0,0.3,0,0.3, "cm")))+
  xlab("Aggregation size in time")+
  ylab("Estimate")




#--------------------------Figure A1(b)---------------------------------------------
ggplot(compare_df[compare_df$mu=="0.5",], aes(x=agg,y=estimate)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=0.5, notch=FALSE,lwd = 0.3)+
  facet_grid(parameter~method,scales = "free_y",labeller = label_parsed)+
  theme_classic() + geom_hline(data = truevalue2, aes(yintercept=truevalue), linetype="dashed", color = "red") +
  theme(axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),strip.text.y = element_text(angle = 0,size=15,margin = margin(0,0.3,0,0.3, "cm")))+
  xlab("Aggregation size in time")+
  ylab("Estimate")





#---------------------Figure 4----------------
fakepoints <- data.frame(c(0.05,0.3,1),c("mu","alpha","beta"))
colnames(fakepoints) <- c("addon","parameter")
fakepoints$parameter <- factor(fakepoints$parameter,c("mu","alpha","beta"))

tmp1 <- t_df[,c(1,3,5,6,7)]
colnames(tmp1) <- c("estimate","CI.length","parameter","agg","mu")
tmp1 <- tmp1[tmp1$mu=="0.3"& tmp1$agg%in%c(0,0.5,1,1.5),]
tmp1$model <- "Without spatial marks"
tmp1 <- tmp1[-c(5)]

tmp2 <- st_df[st_df$par_set=="1",c(1,3,5,6,7)]
colnames(tmp2) <- c("estimate","CI.length","parameter","agg","sagg")
levels(tmp2$agg) <- c("0","0.5","1","1.5")
tmp2 <- tmp2[tmp2$sagg=="1.5"&tmp2$parameter!="gamma",]
tmp2$model <- "With spatial marks"
tmp2 <- tmp2[-c(5)]


tmp3 <- rbind.data.frame(tmp1,tmp2)
compare_df <- tmp3

ggplot(compare_df, aes(x=agg,y=CI.length)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=0.5, notch=FALSE,lwd = 0.3)+
  geom_blank(data=fakepoints, aes(y=addon,x="0"))+
  facet_grid(parameter~model,scales = "free",labeller = labeller(parameter = label_parsed))+
  theme_classic() + 
  theme(axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.y = element_text(angle = 0,size=13,margin = margin(0,0.3,0,0.3, "cm")))+
  xlab("Aggregation size in time")+
  ylab("CI length")




#---------------------------Figure 3----------------------------------

heatmapdata.beta <- st_df[st_df$parameter=="beta" & st_df$par_set=="1",] %>%
  group_by(tagg, sagg,parameter) %>%
  summarise(CI.length.mean = mean(CI.length, na.rm = TRUE)) 

heatmapdata.gamma <- st_df[st_df$parameter== "gamma" & st_df$par_set=="1",] %>%
  group_by(tagg, sagg,parameter) %>%
  summarise(CI.length.mean = mean(CI.length, na.rm = TRUE)) 

p1 <- ggplot(heatmapdata.beta, aes(tagg, sagg, fill= CI.length.mean)) + 
  geom_tile()+
  geom_text(aes(label = round(CI.length.mean, 3)),color = "white",size  =2.8) + 
  scale_fill_continuous( name = "CI length" )+
  labs(x = "Aggregation size for process 1",y = "Aggregation size for process 2")+
  theme_bw()+
  theme(legend.key.size = unit(0.35, "cm"), strip.background = element_blank(),plot.title = element_text(hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 8,margin = margin(t = 6)),
        axis.title.y = element_text(size = 8,margin = margin(r = 6)),axis.ticks = element_blank(),axis.text.x = element_text(margin = margin(t = -3)),axis.text.y = element_text(margin = margin(r = -3)))
p2 <- ggplot(heatmapdata.gamma, aes(tagg, sagg, fill= CI.length.mean)) + 
  geom_tile()+
  geom_text(aes(label = round(CI.length.mean, 3)),color = "white",size = 2.8) + 
  scale_fill_continuous( name = "CI length" )+
  labs(x = "Aggregation size for process 1",y = "Aggregation size for process 2")+
  facet_wrap(~parameter,labeller = label_parsed, ncol = 4)+
  theme_bw()+
  theme(legend.key.size = unit(0.35, "cm"),strip.background = element_blank(),plot.title = element_text(hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 8,margin = margin(t = 6)),
        axis.title.y = element_text(size = 8,margin = margin(r = 6)),axis.ticks = element_blank(),axis.text.x = element_text(margin = margin(t = -3)),axis.text.y = element_text(margin = margin(r = -3)))
p1 <- p1 + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))

p2 <- p2 + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))
p1
p2






#---------------------------Figure 5----------------------------------

heatmapdata.beta <- st2_df[st2_df$parameter%in%c("beta[11]","beta[21]","beta[12]","beta[22]"),] %>%
  group_by(agg1, agg2,parameter) %>%
  summarise(CI.length.mean = mean(CI.length, na.rm = TRUE)) 

heatmapdata.gamma <- st2_df[st2_df$parameter%in%c("gamma[11]","gamma[21]","gamma[12]","gamma[22]"),] %>%
  group_by(agg1, agg2,parameter) %>%
  summarise(CI.length.mean = mean(CI.length, na.rm = TRUE)) 

p1 <- ggplot(heatmapdata.beta, aes(agg1, agg2, fill= CI.length.mean)) + 
  geom_tile()+
  geom_text(aes(label = round(CI.length.mean, 3)),color = "white",size  =2.8) + 
  scale_fill_continuous( name = "CI length" )+
  labs(x = "Aggregation size for process 1",y = "Aggregation size for process 2")+
  facet_wrap(~parameter,labeller = label_parsed, ncol = 4)+
  theme_bw()+
  theme(legend.key.size = unit(0.35, "cm"), strip.background = element_blank(),plot.title = element_text(hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 8,margin = margin(t = 6)),
        axis.title.y = element_text(size = 8,margin = margin(r = 6)),axis.ticks = element_blank(),axis.text.x = element_text(margin = margin(t = -3)),axis.text.y = element_text(margin = margin(r = -3)))
p2 <- ggplot(heatmapdata.gamma, aes(agg1, agg2, fill= CI.length.mean)) + 
  geom_tile()+
  geom_text(aes(label = round(CI.length.mean, 3)),color = "white",size = 2.8) + 
  scale_fill_continuous( name = "CI length" )+
  labs(x = "Aggregation size for process 1",y = "Aggregation size for process 2")+
  facet_wrap(~parameter,labeller = label_parsed, ncol = 4)+
  theme_bw()+
  theme(legend.key.size = unit(0.35, "cm"),strip.background = element_blank(),plot.title = element_text(hjust = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 8,margin = margin(t = 6)),
        axis.title.y = element_text(size = 8,margin = margin(r = 6)),axis.ticks = element_blank(),axis.text.x = element_text(margin = margin(t = -3)),axis.text.y = element_text(margin = margin(r = -3)))
p1 <- p1 + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))

p2 <- p2 + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))
p1
p2




