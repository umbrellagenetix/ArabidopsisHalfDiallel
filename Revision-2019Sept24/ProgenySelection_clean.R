library(dplyr)
library(reshape2)
library(philentropy)
library(ggplot2)
library(linkcomm)
library(caret)

#################### PART 1 :  Genetic Similarity Adjacency Variable #######################
## genotype data processing
## read.table("215k_29parentGenos.hmp.txt",skip=1, header=F, sep="\t")->geno.p

read.table("files/Kmat_parents.txt", skip=3, header=F, sep="\t")->Kp
colnames(Kp)<-c("lines",as.character(Kp$V1))
melt(Kp,id.vars = "lines")->Kpm

## Centered K has negative values - add the min to all values to start at 0
Kpm$value<-Kpm$value+(min(Kpm$value)*-1)

#write.table(Kpm,"p29_K_edgeList.csv",sep=",", row.names = FALSE)

## Kpm = Complete Genetic Distance Graph : 841 edges , 29 parents
## Visualization of the Graph

set.seed(42)
customPal<-c("#FF3300","#33CC33","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#FF9933")
plotType = c("summary","dend","graph")
saveFileName = "output/IBS-20p-spencerCircle.png" # type the file Name to save here
customPal1<-c(customPal,brewer.pal(4,"Dark2"))
library(yarrr) # for piratepal() to use pirate palettes
library(RColorBrewer) # for brewer palettes

###### hierachical clustering and community network analysis
getLinkCommunities(Kpm,hcmethod="ward.D2", directed=FALSE, bipartite = FALSE)->lck.1
plot(lck.1, type = plotType[1])
plot(lck.1, type = plotType[2])
plot(lck.1, type = plotType[3], layout = "spencer.circle")
plot(lck.1, type = plotType[3], layout = layout.fruchterman.reingold)
plot(lck.1, type = plotType[3], layout = layout.kamada.kawai)
# lck.1 stores the information from full matrix- and diesn't identify any communities

## For visualization
lck.1.1 <- newLinkCommsAt(lck.1, cutat = 3.59)
plot(lck.1.1, type = "dend")
plot(lck.1.1, type = "summary")
## lck.1.1 defines 8 cluster communities
## too many edges for visualization
## drop edges that has a value smaller then 0.2

plot(lck.1.1, type = plotType[3], layout = "spencer.circle", pal = customPal, bg = 'grey')
quartz.save(saveFileName, type = "png", device = dev.cur(), dpi = 300)

## 0.20 kinship cut off
Kpm[which(abs(Kpm$value)>0.26),]->kpmr

# Need to remove Bak-2
# missing phenotypic values for Tuescha-9 , ICE61 and Bak-2

kpmr.1<-kpmr[which(kpmr$variable != "Bak-2"),]
kpmr.1<-kpmr.1[which(kpmr.1$lines != "Bak-2"),]

kpmr.1<-kpmr.1[which(kpmr.1$variable != "Tuescha-9"),]
kpmr.1<-kpmr.1[which(kpmr.1$lines != "Tuescha-9"),]

kpmr.1<-kpmr.1[which(kpmr.1$variable != "ICE61"),]
kpmr.1<-kpmr.1[which(kpmr.1$lines != "ICE61"),]

# ###### 10% or more in K cut off
 getLinkCommunities(kpmr.1,hcmethod="ward.D2", directed=FALSE, bipartite = FALSE)->lck
 lck <- newLinkCommsAt(lck, cutat = 2.28)# for 0.25 cut off
plot(lck, type = "dend")
plot(lck, type = "summary")

plot(lck, type = plotType[3], layout = layout.fruchterman.reingold, pal = customPal, bg = 'grey')
plot(lck, type = plotType[3], layout = layout.kamada.kawai, pal = customPal, bg = 'grey')
plot(lck, type = plotType[3], layout = "spencer.circle", pal = customPal, bg = 'grey')

quartz.save(saveFileName, type = "png", device = dev.cur(), dpi = 300)




################# PART 2  : Phenotypic Similarity Adjacency Variable ######################################
## parental phenotypes & hybrid phenotypes files received from Danelle

pms<-read.table("files/2019-05-03_estimated_parental_means.txt", header=T, sep="\t")
hms<-read.table("files/2019-05-03_estimated_hybrid_means.txt", header=T, sep="\t")
hms$cross<-paste(hms$female,hms$male,sep="+")

## phenotypic distance calculation between the parents
## use pms first
## calculate distance only 28 parents out of 29 genotyped

library(philentropy)
pms.1<-pms[,c(5:14)]
rownames(pms.1)<-pms$cross
as.matrix(distance(pms.1, method="euclidean"))->pms.de # euclidean distance
as.matrix(distance(pms.1, method ="jaccard"))->pms.dj # jacquard distance
rownames(pms.de)<-pms$cross
rownames(pms.dj)<-pms$cross
colnames(pms.de)<-pms$cross
colnames(pms.dj)<-pms$cross


pms.dem<-melt(pms.de)
pms.djm<-melt(pms.dj)

# filter the distance matrix for 20% or more
pms.djm.1<-pms.djm[which(pms.djm$value > 0.025),]
pms.dem.1<-pms.dem[which(pms.dem$value > 20),]

# need edgelist for all pairwise parent combinations
## Jacquard distance network
saveFileName = "output/PhenoDist-jd.025-spencerCircle.png"

lcp<-getLinkCommunities(pms.djm.1, hcmethod = "ward.D2", 
                        #                 use.all.edges = FALSE,
                        #                 edglim = 10^4, directed = FALSE, dirweight = 0.5,
                        bipartite = FALSE, dist = NULL, plot = TRUE,
                        check.duplicates = TRUE, removetrivial = TRUE,
                        verbose = TRUE)

lcp.1 <- newLinkCommsAt(lcp, cutat = 1.84) # 8 clusters
plot(lcp.1, type = "summary")
plot(lcp.1, type = "dend")
plot(lcp.1, type = "graph", layout = "spencer.circle", pal = customPal, bg = 'grey')
quartz.save(saveFileName, type = "png", device = dev.cur(), dpi = 300) 

## Euclidean distance network
saveFileName = "output/PhenoDist-ed20-spencerCircle.png"
lce<-getLinkCommunities(pms.dem.1, hcmethod = "ward.D2", 
                        #                 use.all.edges = FALSE,
                        #                 edglim = 10^4, directed = FALSE, dirweight = 0.5,
                        bipartite = FALSE, dist = NULL, plot = TRUE,
                        check.duplicates = TRUE, removetrivial = TRUE,
                        verbose = TRUE)

lce.1 <- newLinkCommsAt(lce, cutat = 2.84) # 6 clusters
plot(lce.1, type = "summary")
plot(lce.1, type = "dend")
plot(lce.1, type = "graph", layout = "spencer.circle")
plot(lce.1, type = "graph", layout = "spencer.circle", pal = customPal, bg = 'grey')
quartz.save(saveFileName, type = "png", device = dev.cur(), dpi = 300) 

##########################################################
###### PART 3 : Evaluation of Transgressiveness ########## 
##########################################################

## Phenotype File processing
phenoMat<-read.table("files/phenotype-matrix.csv",header=T, sep=",")

## phenotypicSimilarity Matrix is DM
## MPH is for mid parent heterosis
## Calculate Frequency of MPH+ and MPH- generated by each cross across across phenotypes.
## First extract out the phenotypes that have MPH tag

mphPhenos<-phenoMat[,which(grepl("MPH",colnames(phenoMat))==TRUE)] # 9 MPH phenotypes
mphPhenos$FID<-phenoMat$FID

regPhenos<-phenoMat[,which(grepl("MPH",colnames(phenoMat))==FALSE)]
regPhenos$FID<-phenoMat$FID

read.table("files/TableS1-Pedigree.csv", header=T, sep=",")->pedFile

## add the parent id'd to both
mph.1<-merge(pedFile,mphPhenos, by.x = "CID",by.y="FID")
reg.1<-merge(pedFile,regPhenos, by.x = "CID",by.y="FID")

## the scales for mph phenos are different- needs to be normalized
## use caret preProcess functions for normalization of each mph phenotype
#############

mph.2<-preProcess(mph.1, method = c("center","scale"),na.remove = TRUE)
predict(mph.2,mph.1)->mph.3

# looks alot better when centered and scaled.
reg.2<-preProcess(reg.1, method = c("center","scale"),na.remove = TRUE)
predict(reg.2,reg.1)->reg.3

### run the Quantile Conversion**** ###
library(StatMeasures)

file.c<-c("mph.3","reg.3")
i=1

# function is modified so the that the NAs are removed prior to quantile calculation
cut_quantile <- function (x) cut(x, quantile(x))
##

for (i in (1:2)){
  xq<- eval(as.name(file.c[i])) #for quantiles
  zq<- eval(as.name(file.c[i])) #for deciles
  j.start=5
  if(i==1){j.end=13}else{j.end=12}
  for( j in (j.start:j.end)){
    if(!is.null(xq[[j]])){
      xq[[j]]<-(as.numeric(cut_quantile(xq[[j]])))
    }else{ 
      next
    }
  }
  assign(paste(file.c[i],"q",sep=""),xq)
  j.start=5
  for(j  in (j.start:j.end)){
    zq[,j]<-decile(zq[,j])
  }
  assign(paste(file.c[i],"d",sep=""),zq)
}

#outputs mph.3d mph.3q reg.3d reg.3q
mph.3d->ds # decile ranking of the mph phenotypes
mph.3q->qs # quartile conversion of the mph phenotypes

## now count the number of samples in deciles across phenotypes given a cross
qs$shp<-0
for (i in(1:length(qs[,1]))){
  qs[i,]->a
  b<-sum(length(which(a[5:13]==1)),length(which(a[5:13]==4)))
  qs$shp[i]<-b/9
}

# the shp column now can be utilized as the frequency of
# obtaining a heterotic phenotype given the cross
## do the high and the low seperately- 
qs$shp.h<-0
qs$shp.l<-0

for (i in(1:length(qs[,1]))){
  qs[i,]->a
  b<-length(which(a[5:13]==1))
  c<-length(which(a[5:13]==4))
  qs$shp.l[i]<-b/9
  qs$shp.h[i]<-c/9
}
## now use the quartile distribution frequency across traits as edgeLists

edgeList.hh<-qs[,c(2:3,15)]
edgeList.lh<-qs[,c(2:3,16)]
edgeList.xh<-qs[,c(2:3,14)]
## run community link analysis

getLinkCommunities(edgeList.hh,hcmethod="ward.D2", directed=FALSE, bipartite = FALSE)->lch
lch <- newLinkCommsAt(lch, cutat = 3.55)
plot(lch, type = "summary")
plot(lch, type = "dend")
plot(lch, type = "graph", layout = "spencer.circle", pal = customPal, bg = 'grey')
quartz.save("output/mph-hhF-spencerCircle.png", type = "png", device = dev.cur(), dpi = 300) 
## heterotic - high transgressive

getLinkCommunities(edgeList.lh,hcmethod="ward.D2", directed=FALSE, bipartite = FALSE)->lcl
lcl <- newLinkCommsAt(lcl, 15cutat = 3.55)
plot(lcl, type = "summary")
plot(lcl, type = "dend")
plot(lcl, type = "graph", layout = "spencer.circle", pal = customPal, bg = 'grey')
quartz.save("output/mph-hiF-spencerCircle.png", type = "png", device = dev.cur(), dpi = 300) 
# hybrid inferior - low transgressive

edgeList.xh<-edgeList.xh[which(edgeList.xh[,3]!=0),] ## remove the edges with 0 weight

getLinkCommunities(edgeList.xh,hcmethod="ward.D2", directed=FALSE, bipartite = FALSE)->lcx
lcx <- newLinkCommsAt(lcx, cutat = 3.9) # 8 clusters
plot(lcx, type = "summary")
plot(lcx, type = "dend")
plot(lcx, type = "graph", layout = "spencer.circle", pal=(customPal))
quartz.save("output/mph-hiF+hhF-spencerCircle.png", type = "png", device = dev.cur(), dpi = 300) 
## high or low transgressive behavior

## Generate a variable to store the cluster mebership profiles based on the edgeLists
lcl$nodeclusters->lcl.nc 
dcast(lcl.nc,node~cluster)->lcl.nc.1
## use node's Cluster Membership profile as a class

for( i in (2:6)){
  mt<-lcl.nc.1[,i]
  mt[which(!is.na(mt))]<-1
  mt[which(is.na(mt))]<-0
  lcl.nc.1[,i]<-mt
}

lcl.nc.1$prf<-c("")
for( j in (1:28)){
  cline<-lcl.nc.1[j,2:6]
  cpr<-cline[1]
  for(c in (1:5)){
    cpr<-paste(cpr,cline[c], sep="")
  }
  lcl.nc.1$prf[j]<-cpr
}
## output lcl.nc.1 - has the cluster mebership of the parents for generating
## F1 with phenotype at the low-heterotic tail for 28 parents
## >unique(lcl.nc.1$prf)
## there are 6 different profiles
## ICE119,Qui-0, ICE216, Rue3-1-31,TueWa1-2 have rare phenotypes
## All the rest are the common profile.
## SAME TREATMENT FOR THE HIGH TAIL membership

getCommunityMatrix(lch)->lch.cm
lch$nodeclusters->lch.nc 
dcast(lch.nc,node~cluster)->lch.nc.1
## use node's Cluster Membership profile as a class

for( i in (2:6)){
  mt<-lch.nc.1[,i]
  mt[which(!is.na(mt))]<-1
  mt[which(is.na(mt))]<-0
  lch.nc.1[,i]<-mt
}

lch.nc.1$prf<-c("")
for( j in (1:28)){
  cline<-lch.nc.1[j,2:6]
  cpr<-cline[1]
  for(c in (1:5)){
    cpr<-paste(cpr,cline[c], sep="")
  }
  lch.nc.1$prf[j]<-cpr
}
## output lch.nc.1
## generates 6 profile classes
## Koch-1,Yeg-1 are the most distant - Mer-6,ICE79, ICE29 ,  Fei-0 ,Cdm-0
## rest are the same

################################
# PART 4a : Progeny Selection - Model family 1 - Continous variables for predictors
################################

Kpm->gd
edgeList.hh ->hh
edgeList.lh->lh
edgeList.xh->xh

gd$cross<-paste(gd$lines,gd$variable,sep="+")
hh$cross<-paste(hh$P1,hh$P2,sep="+")
lh$cross<-paste(lh$P1,lh$P2,sep="+")
xh$cross<-paste(xh$P1,xh$P2,sep="+")

gdh.1<-merge(gd,hh,by="cross")
gdh.2<-merge(gdh.1,lh, by="cross")
gdh.3<-merge(gdh.2,xh, by ="cross")

gdh.4<-gdh.3[,c(1,4,7,10,13)]
## merge that with phenoMat
## Can evaluate the predictive properties of these nodeProfiles
## by merging them with the phenotypes

merge(pedFile,phenoMat, by.x="CID",by.y ="FID")->ppm
pedFile$cross<-paste(pedFile$P1,pedFile$P2,sep="+")

### *** actual merge below gd- the precursor to the gdh.4 is Kpm
### *** this analysis uses the G matrix directly- as the value

merge(pedFile,gdh.4, by="cross")->inputReg
merge(inputReg,ppm, by = "CID")->irp
# write.table(irp,"regressionInputFile.csv",sep=",",row.names = FALSE)
irp.1<-irp[,c(1:4,6:8)]
irp.1 %>% 
  rename(
    P1 = P1.x,
    P2 = P2.x,
    Gs = value,
    HHF = shp.h,
    LHF = shp.l
  ) ->irp.2


## irp.2 contains all variables for the model family 1 except for phenotypic distances
## add pms.djm and pms.dem to irp

pms.djm$cross<-paste(pms.djm$Var1,pms.djm$Var2,sep="+")
pms.dem$cross<-paste(pms.dem$Var1,pms.dem$Var2,sep="+")

pms.djm %>% 
  rename(
    pDist.j = value,
    P1 = Var1,
    P2 = Var2
  ) ->pms.djm

pms.dem %>% 
  rename(
    pDist.e = value,
    P1 = Var1,
    P2 = Var2
  ) ->pms.dem

merge(pms.dem[,3:4],pms.djm[,3:4], by ="cross")->pms.ds
# then merge with the irp.2 by cross

merge(irp.2, pms.ds, by="cross")->inputPred
## this is the input file for Model family 1

dss.1<-inputPred
rownames(dss.1)<-dss.1$CID
dss.1<-dss.1[,-2]

# create the trait list with the column headers from the data files
# hms is the hybrid phenotype file from Danelle Seymour
# pms is the parent phenotype file

traitList<-colnames(hms)[6:15]
parentList<-as.character(pms$female)

# This creates a series of objects named <trait_id>.out for each phenotype
for ( i in (1:length(traitList))){
  c.trait<-traitList[i]
  hp1<-hms[,c(1:3,which(colnames(hms)==c.trait))]
  hp1$p1<-0
  hp1$p2<-0
  # for each trait - merge the parental phenotypes to the  progeny phenotype 
  # according to the parent column value in the dss.1
  for( j in (1:length(parentList))){
    c.pid<-parentList[j]
    pp1<-pms[which(pms$female == c.pid),]
    pp1a<-pp1[,which(colnames(pp1)==c.trait)]
    pp1a->hp1$p1[which(hp1$female == c.pid)]
    
    pp2<-pms[which(pms$male == c.pid),]
    pp2a<-pp1[,which(colnames(pp1)==c.trait)]
    pp2a->hp1$p2[which(hp1$male == c.pid)]
    p.out<-merge(dss.1,hp1,by="cross")
    assign(paste(c.trait,"out",sep="."),p.out)
    # need to store the p.out in its own object
    # the file to be used for prediction the merged dss.1 and hp1 files
    # merge by cross column
  }
}

# fast Methods
methodList1<-c("lm", # linear model
#               "kknn", # k- nearest neighbors
#               "avNNet", # neural Net Multilayer perceptron
#               "elm", #Single Hidden Layer Feedforward Networks - https://rdrr.io/cran/elmNN/man/elmNN-package.html
               "brnn", #Bayesian Regularization for Feed-Forward Neural Networks - https://rdrr.io/cran/brnn/man/brnn.html
               "parRF", # parallel Random Forest
               "cubist", #rulebased model - https://cran.r-project.org/web/packages/Cubist/vignettes/cubist.html - https://rdrr.io/rforge/Cubist/man/cubist.html
               "treebag", #Bootstrap Aggregating Classification-RegressionTrees #https://www.datatechnotes.com/2018/04/classification-with-bagging-treebag.html
               "svmRadialSigma", # Support Vector Regression - Radial Kernel-https://www.kdnuggets.com/2017/03/building-regression-models-support-vector-regression.html
               "svmLinear"
)

# slow methods need to do 3-4 iterations max
methodList2<-c("evtree", # Tree Models from Genetic Algorithms
               "xgbDART" #eXtreme Gradient Boosting 
)


## To test inside the loop do
## ds<-Area_day29_trans.out
setwd("output/MFPlots_MF1_fast")

for( trt in (1:length(traitList))){
  ds<-eval(as.name(paste(traitList[trt],"out",sep=".")))
  list1 <- list(methodList1[1:7])
  list2 <- as.list(methodList1[1:7])
  
  outList<-list(list1,list2)
  
  set.seed(42)
  TestRows     <- c(sample(50,15), sample(50,15)+50, sample(50,15)+100)
  for ( m in (1:length(methodList1))){
    #including phenotypic distance
    TrainData    <- ds[-TestRows,c(4:6,8,12:13)]
    TrainClasses <- ds[-TestRows,11]
    TestData     <- ds[TestRows,c(4:6,8,12:13)]
    TestClasses  <- ds[TestRows,11]
    ## no phenotypic distance - doesn't improve and doesn't reduce the accuracy either
    ## its weird.
    # TrainData    <- ds[-TestRows,c(4:6,12:13)]
    # TrainClasses <- ds[-TestRows,11]
    # TestData     <- ds[TestRows,c(4:6,12:13)]
    # TestClasses  <- ds[TestRows,11]
    
    Fit <- train(x=TrainData, y=TrainClasses,
                 method = methodList1[m],
                 tuneLength = 11,
                 trace = TRUE,
                 maxit = 100)
    outList[[2]][[m]]<-as.list(Fit$results)
    pred.1<-as.data.frame(predict(Fit,TestData))
    rownames(pred.1)->pred.1$CID
    rownames(TestData)->testSampleList
    out.1<-ds[which(rownames(ds) %in% testSampleList),]
    rownames(out.1)->out.1$CID
    merge(out.1,pred.1,by="CID")->xx
    # plot(xx[,13]~xx[,15])
    ## ggplot sctter plot with Error Ribbon
    p1<-ggplot(xx, aes(x=xx[,13], y=xx[,15])) + 
      geom_point()+
      geom_smooth(method=lm, se=TRUE)
    xx$set<-"test"
    colnames(xx)[15]<-"prediction"
    
    pred.2<-as.data.frame(predict(Fit,TrainData))
    rownames(pred.2)->pred.2$CID
    rownames(TrainData)->trainSampleList
    out.2<-ds[which(rownames(ds) %in% trainSampleList),]
    rownames(out.2)->out.2$CID
    merge(out.2,pred.2,by="CID")->yy
    
    p2<-ggplot(yy, aes(x=yy[,13], y=yy[,15])) + 
      geom_point()+
      geom_smooth(method=lm, se=TRUE)
    yy$set <- "train"
    colnames(yy)[15]<-"prediction"
    
    rbind(xx,yy)->zz
    p3<-  ggplot(data=zz, 
             mapping=aes(x=zz[,13], 
                         y=zz[,15])) + 
        geom_point(mapping=aes(colour = set, shape =set)) + 
        geom_smooth(method=lm, se=TRUE) +
        scale_colour_manual(values=customPal[c(1,5)]) +
        labs(title= paste(traitList[trt],methodList1[m],"Predicted vs. Observed Cross Values",sep="-"),
           x=traitList[trt], y = "model_predicted_Values")
      
    ggsave(paste(traitList[trt],methodList1[m],".png",sep=""))
  }
  # store the outList in a named list
  assign(paste(traitList[trt],"outlist",sep="."),outList)
}
#setwd("../")

#tt<-eval(as.name(paste(traitList[j],"outlist",sep=".")))

## Create table for Model Comparison
## Metrics are RMSE, Rquared, MEA and their SDs

results.out<-data.frame(
  trait = character(0L),
  model =character(0L),
  RMSE = numeric(0),
  RMSESD =numeric(0),
  Rsquared =numeric(0),
  RsquaredSD =numeric(0),
  MAE =numeric(0),
  MAESD =numeric(0)
)

for( j in (1:length(traitList))){
  tt<-eval(as.name(paste(traitList[j],"outlist",sep=".")))
  print(traitList[j])
  c.trait<-traitList[j]
  temp.2<-data.frame()
  for( m in 1:length(methodList1)){
    c.model<-methodList1[m]
    c.stats<- as.data.frame(t(as.data.frame(do.call(rbind,tt[[2]][[m]]))))
    ## each models stats are saved in tt[[2]][[i]] per trait
    ## needs to be retrived one-by-one
    vars<-c("Rsquared","RsquaredSD","MAE","MAESD","RMSE","RMSESD")
    temp.1<-(dplyr::select(c.stats,vars))[1,]
    temp.1$model<-c.model
    print(temp.1)
    rbind(temp.2,temp.1)->temp.2
  }
  temp.2$trait<-c.trait
  rbind(results.out,temp.2)->results.out  
}

write.table(results.out,"modelStats_out.csv",sep=",", row.names = FALSE)

## visualize these model comparisons -- Visualization_clean.R script
## To create table/matrices to check for Individual performance indicators
RsqMat.f1m<-matrix(0,nrow=10,ncol=10)
colnames(RsqMat.f1m)<-traitList
rownames(RsqMat.f1m)<-methodList1

MAEmat.f1m<-matrix(0,nrow=8,ncol=10)
colnames(MAEmat.f1m)<-traitList
rownames(MAEmat.f1m)<-methodList1

for( j in (1:length(traitList))){
  tt<-eval(as.name(paste(traitList[j],"outlist",sep=".")))
  print(traitList[j])
  for( i in (1:length(methodList1))){
    print(max(tt[[2]][[i]]$Rsquared, na.rm = TRUE))
    MAEmat.f1m[i,j]<-min(tt[[2]][[i]]$MAE, na.rm = TRUE)
  }
}

################################
# PART 4b : Progeny Selection - Model family 2 - Parent Community Classes for predictors
################################

setwd("output/MFPlots_MF2_fast")
# MF2 will use classes for all variables - to predict actual phenos with regression.
graphList1<-c("lch", # heterosis
             "lcl" # hybrid inferiority
)# 5 clust each
graphList2<-c(
             "lce.1", # pheno- euclidean 
             "lcp.1", # pheno- jacq
             "lck.1.1" # geno
             ) # 8 classes each
## lch and lcl was processed previously into lch.nc.1 and lcl.nc.1
lch.nc.1$var<-"lch"
lcl.nc.1$var<-"lcl"

# select only node, prf and var
lcl.nc.1 %>% 
  select(
    node,
    prf,
    var
  ) ->lcl.nc.2

lch.nc.1 %>% 
  select(
    node,
    prf,
    var
  ) ->lch.nc.2

out<-rbind(lcl.nc.2,lch.nc.2)
  
for(cL in (1:length(graphList2))){
  c.graph=graphList2[cL]
  c.val<-eval(as.name(c.graph))
  c.val$nodeclusters->c.val.nc 
  cluster.count<-length(c.val$clusters)
  dcast(c.val.nc,node~cluster)->c.val.nc.1
  ## use node's Cluster Membership profile as a class
  for(i in (2:(cluster.count+1))){
    mt<-c.val.nc.1[,i]
    mt[which(!is.na(mt))]<-1
    mt[which(is.na(mt))]<-0
    c.val.nc.1[,i]<-mt
  }
  c.val.nc.1$prf<-c("")
  for( j in (1:length(c.val.nc.1$node))){ # number of parents in the list
    cline<-c.val.nc.1[j,2:(cluster.count+1)]
    cpr<-cline[1]
    for(c in (1:cluster.count)){
      cpr<-paste(cpr,cline[c], sep="")
    }
    c.val.nc.1$prf[j]<-cpr
  }
  c.val.nc.1$var<-graphList2[cL]
  c.val.nc.1 %>% select(node,prf,var)->nc.2
  assign(paste(c.graph,".nc.1",sep=""),c.val.nc.1)
  out = rbind(out,nc.2)
}

## then cast out
dcast(out,node~var , value.var = "prf")->out.dc
#missing values for Tuescha-9 , ICE61 and Bak-2
## drop those would create problems in training
out.dc.1<-na.exclude(out.dc) # 27 parents here
colnames(out.dc.1)<-c("node","pde","hh","gs","hl","pdj")  
# start with this add the pheno columns from pms file
merge(out.dc.1,pms, by.x = "node",by.y="female")->pcr
# use this to add variables to the cross table
## write.table(pcr,"classVariables.csv",sep=",")
## now add parent phenos into the cross phenos file
traitList<-colnames(hms)[6:15]
parentList<-as.character(pcr$node)

# This creates a series of objects named <trait_id>.mf2.out for each cross phenotype
for ( i in (1:length(traitList))){
  c.trait<-traitList[i]
  hp1<-hms[,c(1:3,which(colnames(hms)==c.trait))]
  hp1$p1.pheno<-NA # parent 1 pheno
  hp1$p2.pheno<-NA # parent 2 pheno
  hp1$p1.euc<-NA # Parent 1 euclidean phenotypic distance class
  hp1$p1.jdc<-NA # Parent 2 euclidean phenotypic distance class
  hp1$p2.euc<-NA # Parent 1 jacquard phenotypic distance class
  hp1$p2.jdc<-NA # Parent 2 jacquard phenotypic distance class
  hp1$p1.gsc<-NA # Parent 1 genetic similarity class
  hp1$p2.gsc<-NA # Parent 2 genetic similarity class
  hp1$p1.lclc<-NA # Parent 1 hybrid inferiority class 
  hp1$p2.lclc<-NA # Parent 2 hybrid inferiority class
  hp1$p1.lchc<-NA # Parent 1 heterosis class
  hp1$p2.lchc<-NA # Parent 2 heterosis class
  
  # for each trait - merge the parental phenotypes to the  progeny phenotype 
  # according to the parent column value in the pcr
  ## need to go over the parent list twice once for p1 once for p2
  for( j1 in (1:length(parentList))){
    c.pid<-parentList[j1]
    pp1<-pcr[which(pcr$node == c.pid),]
    pp1a<-pp1[,which(colnames(pp1)==c.trait)]
    pp1a->hp1$p1.pheno[which(hp1$female == c.pid)]
    pp1a->hp1$p2.pheno[which(hp1$male == c.pid)]
    
    pp1$pde->hp1$p1.euc[which(hp1$female == c.pid)]
    pp1$pdj->hp1$p1.jdc[which(hp1$female == c.pid)]
    pp1$hl->hp1$p1.lclc[which(hp1$female == c.pid)]
    pp1$hh->hp1$p1.lchc[which(hp1$female == c.pid)]
    pp1$gs->hp1$p1.gsc[which(hp1$female == c.pid)]
    pp1$pde->hp1$p2.euc[which(hp1$male == c.pid)]
    pp1$pdj->hp1$p2.jdc[which(hp1$male == c.pid)]
    pp1$hl->hp1$p2.lclc[which(hp1$male == c.pid)]
    pp1$hh->hp1$p2.lchc[which(hp1$male == c.pid)]
    pp1$gs->hp1$p2.gsc[which(hp1$male == c.pid)]
  }
  hp1.1<-na.exclude(hp1)
  assign(paste(c.trait,"mf2.out",sep="."),hp1.1)
}
## now do the model fit with the output files

## To test inside the loop do
## ds<-Area_day29_trans.mf2.out
methodList1<-c("lm", # linear model
               #               "kknn", # k- nearest neighbors
               #               "avNNet", # neural Net Multilayer perceptron
               #               "elm", #Single Hidden Layer Feedforward Networks - https://rdrr.io/cran/elmNN/man/elmNN-package.html
               #"brnn", #Bayesian Regularization for Feed-Forward Neural Networks - https://rdrr.io/cran/brnn/man/brnn.html
               "parRF", # parallel Random Forest
               "cubist" #rulebased model - https://cran.r-project.org/web/packages/Cubist/vignettes/cubist.html - https://rdrr.io/rforge/Cubist/man/cubist.html
#               "treebag" #Bootstrap Aggregating Classification-RegressionTrees #https://www.datatechnotes.com/2018/04/classification-with-bagging-treebag.html
#               "svmRadialSigma" # Support Vector Regression - Radial Kernel-https://www.kdnuggets.com/2017/03/building-regression-models-support-vector-regression.html
#               "svmLinear"
)

for( trt in (1:length(traitList))){
  ds<-eval(as.name(paste(traitList[trt],"mf2.out",sep=".")))
  list1 <- list(methodList1[1:4])
  list2 <- as.list(methodList1[1:4])
  
  outList<-list(list1,list2)
  
  set.seed(42)
  TestRows     <- c(sample(50,15), sample(50,15)+50, sample(50,15)+100)
  for ( m in (1:length(methodList1))){
    #including phenotypic distance
    TrainData    <- ds[-TestRows,c(5:16)]
    TrainClasses <- ds[-TestRows,4]
    TestData     <- ds[TestRows,c(5:16)]
    TestClasses  <- ds[TestRows,4]
    ## no phenotypic distance - doesn't improve and doesn't reduce the accuracy either
    ## its weird.
    # TrainData    <- ds[-TestRows,c(4:6,12:13)]
    # TrainClasses <- ds[-TestRows,11]
    # TestData     <- ds[TestRows,c(4:6,12:13)]
    # TestClasses  <- ds[TestRows,11]
    
    Fit <- train(x=TrainData, y=TrainClasses,
                 method = methodList1[m],
                 tuneLength = 11,
                 trace = TRUE,
                 maxit = 100)
    outList[[2]][[m]]<-as.list(Fit$results)
    pred.1<-as.data.frame(predict(Fit,TestData))
    rownames(pred.1)->pred.1$CID
    rownames(TestData)->testSampleList
    out.1<-ds[which(rownames(ds) %in% testSampleList),]
    rownames(out.1)->out.1$CID
    merge(out.1,pred.1,by="CID")->xx
    # plot(xx[,13]~xx[,15])
    ## ggplot sctter plot with Error Ribbon
    p1<-ggplot(xx, aes(x=xx[,5], y=xx[,18])) + 
      geom_point()+
      geom_smooth(method=lm, se=TRUE)
    xx$set<-"test"
    colnames(xx)[18]<-"prediction"
    
    pred.2<-as.data.frame(predict(Fit,TrainData))
    rownames(pred.2)->pred.2$CID
    rownames(TrainData)->trainSampleList
    out.2<-ds[which(rownames(ds) %in% trainSampleList),]
    rownames(out.2)->out.2$CID
    merge(out.2,pred.2,by="CID")->yy
    
    p2<-ggplot(yy, aes(x=yy[,5], y=yy[,18])) + 
      geom_point()+
      geom_smooth(method=lm, se=TRUE)
    yy$set <- "train"
    colnames(yy)[18]<-"prediction"
    
    rbind(xx,yy)->zz
    p3<-  ggplot(data=zz, 
                 mapping=aes(x=zz[,5], 
                             y=zz[,18])) + 
      geom_point(mapping=aes(colour = set, shape =set)) + 
      geom_smooth(method=lm, se=TRUE) +
      scale_colour_manual(values=customPal[c(1,5)]) +
      labs(title= paste(traitList[trt],methodList1[m],"Predicted vs. Observed Cross Values",sep="-"),
           x=traitList[trt], y = "model_predicted_Values")
    
    ggsave(paste(traitList[trt],methodList1[m],".png",sep=""))
  }
  # store the outList in a named list
  assign(paste(traitList[trt],"outlist",sep="."),outList)
}

results.out<-data.frame(
  trait = character(0L),
  model =character(0L),
  RMSE = numeric(0),
  RMSESD =numeric(0),
  Rsquared =numeric(0),
  RsquaredSD =numeric(0),
  MAE =numeric(0),
  MAESD =numeric(0)
)

for( j in (1:length(traitList))){
  tt<-eval(as.name(paste(traitList[j],"outlist",sep=".")))
  print(traitList[j])
  c.trait<-traitList[j]
  temp.2<-data.frame()
  for( m in 1:length(methodList1)){
    c.model<-methodList1[m]
    c.stats<- as.data.frame(t(as.data.frame(do.call(rbind,tt[[2]][[m]]))))
    ## each models stats are saved in tt[[2]][[i]] per trait
    ## needs to be retrived one-by-one
    vars<-c("Rsquared","RsquaredSD","MAE","MAESD","RMSE","RMSESD")
    temp.1<-(dplyr::select(c.stats,vars))[1,]
    temp.1$model<-c.model
    print(temp.1)
    rbind(temp.2,temp.1)->temp.2
  }
  temp.2$trait<-c.trait
  rbind(results.out,temp.2)->results.out  
}

write.table(results.out,"mf2_modelStats_out.csv",sep=",", row.names = FALSE)
