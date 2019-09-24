########################
#### QUANTILE CONVERTION
########################
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