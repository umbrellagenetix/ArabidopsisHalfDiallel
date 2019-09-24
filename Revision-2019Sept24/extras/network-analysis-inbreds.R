library(linkcomm)

setwd("~/Desktop/Archive_April2019/DL282g2x-validation/inbred_mlm_k/stats")
#setwd("/home/elhan/trifecta/DL282g2x-validation/round1_mlm_k/stats")
#⁨Desktop⁩ ▸ ⁨Archive_April2019⁩ ▸ ⁨DL282g2x-valida⁨round1_mlm_k⁩ ▸ ⁨stats⁩
#setwd("D:\\Elhan\\GECCO")
list.files(pattern="outAnn-HsG.csv")->fileList
fileList<-fileList[c(2:10,1)]

# inbred results are in "/home/elhan/trifecta/DL282g2x-validation/round1_mlm_k/stats
## need to be parsed for the chXoutAnn-HsG.csv format before can be processed.
## process the hybrids first.

# read them all and extract columns, 1,2,3,23,and 24
# process each chromosome seperately- then merge results as needed

#initialize for joining chromosomes
all_chrom<-data.frame(
  Trait=character(0),
  Marker=character(0),
  Chr=integer(0),
  HSID=character(0),
  GeneID=character(0)
)



for ( i in c(1:10)){
  read.table(fileList[i], header=T, sep=",")->a 
  a1<-a[,c(1:3,23:24)]
  rbind(all_chrom,a1)->all_chrom
}

all_chrom$cHSID<-paste(all_chrom$Chr,all_chrom$HSID,sep=".")
all_chrom<-all_chrom[,c(2,6,5,4,3,1)]

##only c10
read.table(fileList[10], header=T, sep=",")->a 
a1<-a[,c(1:3,23:24)]
c10<-a1
c10$cHSID<-paste(c10$Chr,c10$HSID,sep=".")
c10<-c10[,c(2,6,5,4,3,1)]

for ( j in 1:3){
  mta<-unique(all_chrom[,c(j,6)])
  #mta<-unique(c10[,c(j,6)])
  
  edgeList<-mta
  
  getLinkCommunities(edgeList,hcmethod="ward.D2", directed=TRUE, bipartite = TRUE)->lc1
  lc1.1 <- newLinkCommsAt(lc1, cutat = 1.3)
  plot(lc1.1, type = "summary")
  plot(lc1.1, type = "dend")
  plot(lc1.1, type = "graph", layout = "spencer.circle")
  lc1.1$edges[which(lc1.1$edges["node2"]=="PltYield" |
                      lc1.1$edges["node2"]=="LFWDT"
  ),]
  #query by hs ex-10.H23
  lc1.1$edges[which(lc1.1$edges["node1"]=="10.H17"),]->s1
  cList.1<-as.numeric(as.character(s1$cluster))
  
  plot(lc1.1, type = "graph", layout = layout.fruchterman.reingold, shownodesin =3)
  plot(lc1, type = "graph", shownodesin = 2, node.pies = TRUE)
  plot(lc1.1, type = "graph", layout="spencer.circle", node.pies = TRUE, 
       clusterids=cList.1)
  
  #query by cluster that contains "PltYield"
  lc1.1$edges[which(lc1.1$edges["node2"]=="PltYield"),]->s2
  cList.2<-as.numeric(as.character(s2$cluster))
  plot(lc1.1, type = "graph", layout="spencer.circle", node.pies = TRUE, 
       clusterids=cList.2)
  
  # generate list of HS contributing to this trait
  lc1.1$edges->eLT
  for( c in 1:3){
    eLT[,c]<-as.character(eLT[,c])
  }
  unique(as.character(s2$node1))->hList
  
  eLT[which(eLT[,1] %in% hList),]->s2a
  
  # Find Given trait- what the regions are-
  # found it there is 42 regions= #nodes in node1 in s2a
  # then ask- given a region, for those 42 regions- what are the other traits 
  # they are associated with?
  for( h in 1:length(hList)){
    hList[h]->roi
    c.hs<-s2a[which(s2a$node1==roi),]
    print(c.hs)
  }
  
  getLinkCommunities(s2a,hcmethod="ward.D2", directed=TRUE, bipartite = TRUE)->lc1a
  lc1a.1 <- newLinkCommsAt(lc1a, cutat = 2.65)
  plot(lc1a.1, type = "summary")
  plot(lc1a.1, type = "dend")
  plot(lc1a.1, type = "graph", layout = "spencer.circle")
  plot(lc1a.1, type = "graph", layout="spencer.circle", node.pies = TRUE, 
       clusterids=c(12))
  
  
  #query byquery by cluster that contains "PltYield" then by chromosome eg.10.
  s2[which(grepl("10.",s2$node1)==TRUE),]->s3
  cList.3<-as.numeric(as.character(s3$cluster))
  plot(lc1.1, type = "graph", layout="spencer.circle", node.pies = TRUE, 
       clusterids=cList.3)
  
  edgeList<-unique(res.2)
  getLinkCommunities(edgeList,hcmethod="ward", directed=TRUE, bipartite = TRUE)->lc2
  plot(lc2, type = "graph", layout = "spencer.circle", shownodesin = 2)
  plot(lc2, type = "graph", layout = layout.fruchterman.reingold, shownodesin =2)
  plot(lc2, type = "graph", shownodesin = 2, node.pies = TRUE)
  plot(lc2, type = "graph", layout="spencer.circle",shownodesin = 2, node.pies = TRUE, clusterids=c(8,10,18,31,35,36))
  
  edgeList<-res.3
  getLinkCommunities(edgeList,hcmethod="ward.D2", directed=TRUE, bipartite = TRUE)->lc3
  lc3.1 <- newLinkCommsAt(lc3, cutat = 2.51)
  
  plot(lc3, type = "graph", layout = "spencer.circle", shownodesin = 2)
  plot(lc2, type = "graph", layout = layout.fruchterman.reingold, shownodesin =2)
  plot(lc2, type = "graph", shownodesin = 2, node.pies = TRUE)
  plot(lc2, type = "graph", layout="spencer.circle",shownodesin = 2, node.pies = TRUE, clusterids=c(8,10,18,31,35,36))
  
}
