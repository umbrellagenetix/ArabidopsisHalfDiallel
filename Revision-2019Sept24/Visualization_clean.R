####################
## FIGURES- VISUALIZATION
####################
## visualize cq.df
# #setwd("dist-s3") # for S3 figures- before quartile cut
# setwd("../visualization_output/")
# #for (i  in (1:1)){ # incase you would like to ptocess multiple files
# eval(cq.df[,-26])->xc #transformed data
# eval(x)->d.df #original data <--- change here to change data
# d.mat<-as.matrix(d.df)
# 
# for (j in (1:25)){
#   a<-matrix(nrow=length(xc[,1]),ncol=2)
#   d.mat[,j]->a[,1]
#   xc[,j]->a[,2]
#   colnames(a)<-c("original","predict")
#   a.f <-as.data.frame(a)
#   header_list<-colnames(xc)[1:length(d.mat)]
#   
#   #write.table(summary(a),paste("F",i,as.character(lookup$V1[j]),"summary.txt",sep="."))
#   
#   p1<- ggplot(a.f, aes(sample= original)) + stat_qq() + #geom_abline(min(predict), "blue") + 
#     ggtitle("Original Data") + labs(x="Theoretical",y="Observed")
#   
#   p2<- ggplot(a.f, aes(sample= predict)) + stat_qq() + ggtitle("Treated Data") + labs(x="Theoretical",
#                                                                                       y="Observed")
#   
#   #Histograms-with ggplot2
#   p1h<- ggplot(a.f, aes(original)) + geom_histogram()  + labs(x="Standardized Trait Value", y="Frequency") 
#   p2h<- ggplot(a.f, aes(predict)) + geom_histogram() + labs(x="Standardized Trait Value", y="Frequency")
#   
#   g <- arrangeGrob(p1, p2, p1h, p2h, nrow=2) #generates g
#   
#   ggsave(file=paste(as.name(colnames(xc)[j]),"exp","pdf",sep="."), g) #saves g
# }
# #}
# setwd("../")

############################################  
## MODEL PERFORMANCE COMPARISON ACROSS TRAITS
## BOX 1. Arabidopsis Case Study
############################################
## paths are written assuming that the starting directory is the "output" folder from the repo
read.table("files/modelStats_out.csv",sep=",", header=T)->results.out
hs<-read.table("files/heritabilities.csv", sep=",",header=T)
results.out<-rbind(results.out,hs)

## 7 ML models + 2 h2 estimates visualized :90 observations
#enumerate the models and traits
modelFrame<-data.frame(
  mID= numeric(9),
  model = character(9L)
)

modelFrame$model<-unique(results.out$model)
modelFrame$mID<-c(1:9)

trtFrame<-data.frame(
  trtID= numeric(10),
  trait = character(10L)
)

trtFrame$trait<-unique(results.out$trait)
trtFrame$trtID<-c(1:10)

# add traitID and modelID to the results
results.out$trtID<-0
results.out$mID<-0

data_avgs <- results.out %>%
  group_by(trait) %>% 
  summarize(Avg.rq = mean(Rsquared),
            SD.rq = mean(RsquaredSD),
            N = n()) %>% 
  ungroup() %>% 
  mutate(L_CI = Avg.rq - SD.rq ) %>% 
  mutate(U_CI = Avg.rq  + SD.rq )

data_avgs<-data_avgs[order(data_avgs$Avg.rq),]

#initialize new df
data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.1

for( i in (1:length(trtFrame$trtID))){
  c.trait<-trtFrame$trait[i]
  c.trtID<-trtFrame$trtID[i]
  c.set<-results.out[which(results.out$trait == c.trait),]
  c.set$trtID<-c.trtID
  res.1<-rbind(res.1,c.set)
}

data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.2

for( i in (1:length(modelFrame$mID))){
  c.model<-modelFrame$model[i]
  c.mID<-modelFrame$mID[i]
  c.set<-res.1[which(res.1$model == c.model),]
  c.set$mID<-c.mID
  res.2<-rbind(res.2,c.set)
}

# order levels of the trait vars explicitly based on data_avgs
# trait.order <- factor(data_avgs$trait,
#                levels = data_avgs$trait)
## order traits based on Avg.rq from data_avgs

data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.3

for( i in (1:length(trtFrame$trtID))){
  c.trait<-trtFrame$trait[i]
  c.trtID<-trtFrame$trtID[i]
  c.set<-results.out[which(results.out$trait == c.trait),]
  c.set$Avg.rq<-data_avgs$Avg.rq[which(data_avgs$trait==c.trait)]
  res.3<-rbind(res.3,c.set)
}

res.3<-res.3[order(res.3$Avg.rq,res.3$Rsquared),]
## the trait levels needs to be ordered- not the trait itself

trait.order <- factor(res.3$trait,
                      levels = data_avgs$trait)
res.3$tro<-trait.order

####
## Add the NarrowSense and Braodsense heritabilities per trait
customPal1<-c("#FF3300", "#33CC33", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", 
             "#225EA8", "#FF9933", "#1B9E77", "#D95F02", "#7570B3","#E7298A")
customPal2<-customPal1[c(9,2,1,4:8,3)]
customPal3<-c("#1B9E77", "#33CC33", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8",
              "#C7E9B4", "#FFC733", "#FF3300")

viz.data<- res.3

bp <- ggplot(viz.data, aes(x=tro, y=Rsquared, fill=model)) 
p1<- bp + geom_bar(position=position_dodge(), stat="identity", color="black") 
p1<- p1 + scale_fill_manual(values=customPal3, 
                            name="model")
p1<- p1+ theme(plot.title = element_text(size = 12, face = "bold"), 
               axis.text.x = element_text(size = 12, face = "bold"),
               axis.text.y = element_text(size = 12, face = "bold"),
               legend.text = element_text(size = 12, face = "bold")
               ) 
p1<- p1 + geom_pointrange(aes(ymin=(Rsquared-RsquaredSD), ymax=(Rsquared+RsquaredSD)),
                    position=position_dodge(.9))
p1<- p1 + geom_errorbar(aes(ymin=(Rsquared-RsquaredSD), ymax=(Rsquared+RsquaredSD)),
                    position=position_dodge(.9))
p1 <- p1 + labs(y = "Fraction of Phenotypic Variance Explained",x=NULL,
                aes(element_text(size =11, face = bold)))#, x = "Model Rsq by Traits") 
p1 <- p1 + ggtitle("R-squared Estimates & Their Standard Deviation per Trait") +
  theme(legend.position = "top")
#p1 <- p1 + coord_flip() 
print(p1)

ggsave("MF1_Rsq.png",device="png", width = 48, height = 17, units = "cm")
quartz.save("MethodComparison.png", type = "png", device = dev.cur(), dpi = 300) 
## need to work on aesthetics a bit
## rename the traits and models and such
## also add the heritabilities from Seynour et al before plotting.
## its good for now.


### MF2- visualization

results.out<-read.table("MFPlots_MF2_fast/mf2_modelStats_out.csv",sep=",", header=T)
#results.out<-read.table("mf2_modelStats_out.csv",sep=",", header=T)

## 3 ML models + 2 h2 estimates visualized :50 observations
#enumerate the models and traits
modelFrame<-data.frame(
  mID= numeric(5),
  model = character(5L)
)

modelFrame$model<-unique(results.out$model)
modelFrame$mID<-c(1:5)

trtFrame<-data.frame(
  trtID= numeric(10),
  trait = character(10L)
)

trtFrame$trait<-unique(results.out$trait)
trtFrame$trtID<-c(1:10)

# add traitID and modelID to the results
results.out$trtID<-0
results.out$mID<-0

data_avgs <- results.out %>%
  group_by(trait) %>% 
  summarize(Avg.rq = mean(Rsquared),
            SD.rq = mean(RsquaredSD),
            N = n()) %>% 
  ungroup() %>% 
  mutate(L_CI = Avg.rq - SD.rq ) %>% 
  mutate(U_CI = Avg.rq  + SD.rq )

data_avgs<-data_avgs[order(data_avgs$Avg.rq),]

#initialize new df
data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.1

for( i in (1:length(trtFrame$trtID))){
  c.trait<-trtFrame$trait[i]
  c.trtID<-trtFrame$trtID[i]
  c.set<-results.out[which(results.out$trait == c.trait),]
  c.set$trtID<-c.trtID
  res.1<-rbind(res.1,c.set)
}

data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.2

for( i in (1:length(modelFrame$mID))){
  c.model<-modelFrame$model[i]
  c.mID<-modelFrame$mID[i]
  c.set<-res.1[which(res.1$model == c.model),]
  c.set$mID<-c.mID
  res.2<-rbind(res.2,c.set)
}

# order levels of the trait vars explicitly based on data_avgs
# trait.order <- factor(data_avgs$trait,
#                levels = data_avgs$trait)
## order traits based on Avg.rq from data_avgs

data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.3

for( i in (1:length(trtFrame$trtID))){
  c.trait<-trtFrame$trait[i]
  c.trtID<-trtFrame$trtID[i]
  c.set<-results.out[which(results.out$trait == c.trait),]
  c.set$Avg.rq<-data_avgs$Avg.rq[which(data_avgs$trait==c.trait)]
  res.3<-rbind(res.3,c.set)
}

res.3<-res.3[order(res.3$Avg.rq,res.3$Rsquared),]
## the trait levels needs to be ordered- not the trait itself

trait.order <- factor(res.3$trait,
                      levels = data_avgs$trait)
res.3$tro<-trait.order

####
## Add the NarrowSense and Braodsense heritabilities per trait
customPal1<-c("#FF3300", "#33CC33", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", 
             "#225EA8", "#FF9933", "#1B9E77", "#D95F02", "#7570B3","#E7298A")
customPal2<-customPal1[c(9,2,1,4:8,3)]
customPal3<-c("#1B9E77", "#33CC33", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8",
              "#C7E9B4", "#FFC733", "#FF3300")
customPal4<-c("#FFC733", "#1D91C0", "#225EA8","#FF3300","#C7E9B4")

viz.data<- res.3

bp <- ggplot(viz.data, aes(x=tro, y=Rsquared, fill=model)) 
p1<- bp + geom_bar(position=position_dodge(), stat="identity", color="black") 
p1<- p1 + scale_fill_manual(values=customPal4, 
                            name="model")
p1<- p1+ theme(plot.title = element_text(size = 12, face = "bold"), 
               axis.text.x = element_text(size = 12, face = "bold"),
               axis.text.y = element_text(size = 12, face = "bold"),
               legend.text = element_text(size = 12, face = "bold")
               ) 
p1<- p1 + geom_pointrange(aes(ymin=(Rsquared-RsquaredSD), ymax=(Rsquared+RsquaredSD)),
                    position=position_dodge(.9))
p1<- p1 + geom_errorbar(aes(ymin=(Rsquared-RsquaredSD), ymax=(Rsquared+RsquaredSD)),
                    position=position_dodge(.9))
p1 <- p1 + labs(y = "Fraction of Phenotypic Variance Explained",x=NULL,
                aes(element_text(size =12, face = bold)))#, x = "Model Rsq by Traits") 
p1 <- p1 + ggtitle("R-squared Estimates & Their Standard Deviation per Trait")
#p1 <- p1 + coord_flip() 
print(p1)



ggsave("MF2_Rsq.png",device="png")
quartz.save("MethodComparison.png", type = "png", device = dev.cur(), dpi = 300) 
################
#Broadsense heritability always do best- check the rankings against narorow sense

results.out<-read.table("MFPlots_MF2_fast/mf2_noH2_modelStats_out.csv",sep=",", header=T)
#results.out<-read.table("mf2_noH2_modelStats_out.csv",sep=",", header=T)

## 3 ML models + 1 h2 estimates visualized :40 observations
#enumerate the models and traits
modelFrame<-data.frame(
  mID= numeric(4),
  model = character(4L)
)

modelFrame$model<-unique(results.out$model)
modelFrame$mID<-c(1:4)

trtFrame<-data.frame(
  trtID= numeric(10),
  trait = character(10L)
)

trtFrame$trait<-unique(results.out$trait)
trtFrame$trtID<-c(1:10)

# add traitID and modelID to the results
results.out$trtID<-0
results.out$mID<-0

data_avgs <- results.out %>%
  group_by(trait) %>% 
  summarize(Avg.rq = mean(Rsquared),
            SD.rq = mean(RsquaredSD),
            N = n()) %>% 
  ungroup() %>% 
  mutate(L_CI = Avg.rq - SD.rq ) %>% 
  mutate(U_CI = Avg.rq  + SD.rq )

data_avgs<-data_avgs[order(data_avgs$Avg.rq),]

#initialize new df
data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.1

for( i in (1:length(trtFrame$trtID))){
  c.trait<-trtFrame$trait[i]
  c.trtID<-trtFrame$trtID[i]
  c.set<-results.out[which(results.out$trait == c.trait),]
  c.set$trtID<-c.trtID
  res.1<-rbind(res.1,c.set)
}

data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.2

for( i in (1:length(modelFrame$mID))){
  c.model<-modelFrame$model[i]
  c.mID<-modelFrame$mID[i]
  c.set<-res.1[which(res.1$model == c.model),]
  c.set$mID<-c.mID
  res.2<-rbind(res.2,c.set)
}

# order levels of the trait vars explicitly based on data_avgs
# trait.order <- factor(data_avgs$trait,
#                levels = data_avgs$trait)
## order traits based on Avg.rq from data_avgs

data.frame(
  Rsquared=numeric(0),
  RsquaredSD =numeric(0),
  MAE=numeric(0),
  MAESD =numeric(0),
  RMSE=numeric(0),
  RMSESD =numeric(0),
  model =character(0L),
  trait = character(0L),
  trtID = numeric(0),
  mID = numeric(0)
)->res.3

for( i in (1:length(trtFrame$trtID))){
  c.trait<-trtFrame$trait[i]
  c.trtID<-trtFrame$trtID[i]
  c.set<-results.out[which(results.out$trait == c.trait),]
  c.set$Avg.rq<-data_avgs$Avg.rq[which(data_avgs$trait==c.trait)]
  res.3<-rbind(res.3,c.set)
}

res.3<-res.3[order(res.3$Avg.rq,res.3$Rsquared),]
## the trait levels needs to be ordered- not the trait itself

trait.order <- factor(res.3$trait,
                      levels = data_avgs$trait)
res.3$tro<-trait.order

####
## Add the NarrowSense and Braodsense heritabilities per trait
customPal1<-c("#FF3300", "#33CC33", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", 
             "#225EA8", "#FF9933", "#1B9E77", "#D95F02", "#7570B3","#E7298A")
customPal2<-customPal1[c(9,2,1,4:8,3)]
customPal3<-c("#1B9E77", "#33CC33", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8",
              "#C7E9B4", "#FFC733", "#FF3300")
customPal4<-c("#FFC733", "#1D91C0", "#225EA8","#FF3300","#C7E9B4")

viz.data<- res.3

bp <- ggplot(viz.data, aes(x=tro, y=Rsquared, fill=model)) 
p1<- bp + geom_bar(position=position_dodge(), stat="identity", color="black") 
p1<- p1 + scale_fill_manual(values=customPal4, 
                            name="model")
p1<- p1+ theme(plot.title = element_text(size = 12, face = "bold"), 
               axis.text.x = element_text(size = 12, face = "bold"),
               axis.text.y = element_text(size = 12, face = "bold"),
               legend.text = element_text(size = 12, face = "bold")
               ) 
p1<- p1 + geom_pointrange(aes(ymin=(Rsquared-RsquaredSD), ymax=(Rsquared+RsquaredSD)),
                    position=position_dodge(.9))
p1<- p1 + geom_errorbar(aes(ymin=(Rsquared-RsquaredSD), ymax=(Rsquared+RsquaredSD)),
                    position=position_dodge(.9))
p1 <- p1 + labs(y = "Fraction of Phenotypic Variance Explained",x=NULL,
                aes(element_text(size =12, face = bold)))#, x = "Model Rsq by Traits") 
p1 <- p1 + ggtitle("R-squared Estimates & Their Standard Deviation per Trait")
#p1 <- p1 + coord_flip() 
print(p1)



ggsave("noH2_MF2_Rsq.png",device="png")
quartz.save("noH2_MethodComparison.png", type = "png", device = dev.cur(), dpi = 300) 



