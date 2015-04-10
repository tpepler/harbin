# Function to check convergence of boundary values

convergence.check<-function(datavec){#,graph=c("density","histogram")){

  # datavec: vector of normalised GOI values
  # graph: type of graph to be used for representation of distribution of normalised GOI values
  
  nGOIvals<-length(datavec)
  nbounds<-7
  boundvals<-quantile(datavec,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  if(nGOIvals>9){
    par(mfrow=c(3,1))
  }
  else {
    par(mfrow=c(1,1))
  }
  #if(graph[1]=="histogram"){
    #hist(datavec,breaks=20,col="darkgrey",main="Gene of interest: Distribution",xlab="Normalised GOI value",xlim=c(0,max(datavec)+0.2))
  #}
  #if(graph[1]=="density"){
    d<-density(datavec,na.rm=TRUE)
    hist(datavec,
         col="darkgray",
         main="Gene of interest: Distribution",
         xlab="Normalised GOI value",
         xlim=c(0,max(d$x)),
         freq=FALSE)
    #plot(d,main="Gene of interest: Distribution",xlab="Normalised GOI value",xlim=c(0,max(d$x)))
    lines(d,lwd=2)
    #polygon(d,col="darkgrey")
  #}
  for(i in 1:nbounds){
    abline(v=boundvals[i],lty=3)
  }
    
  if(nGOIvals>9){
    class1UL<-rep(NA,times=nGOIvals-9)
    class2UL<-rep(NA,times=nGOIvals-9)
    class3UL<-rep(NA,times=nGOIvals-9)
    class4UL<-rep(NA,times=nGOIvals-9)
    class5UL<-rep(NA,times=nGOIvals-9)
    class1mean<-rep(NA,times=nGOIvals-9)
    class2mean<-rep(NA,times=nGOIvals-9)
    class3mean<-rep(NA,times=nGOIvals-9)
    class4mean<-rep(NA,times=nGOIvals-9)
    class5mean<-rep(NA,times=nGOIvals-9)
    for(j in 10:nGOIvals){
      tempdata<-datavec[1:j]
      tempdata<-tempdata[tempdata>=0]
      class1UL[j-9]<-quantile(tempdata,probs=0.2,na.rm=TRUE)
      class2UL[j-9]<-quantile(tempdata,probs=0.4,na.rm=TRUE)
      class3UL[j-9]<-quantile(tempdata,probs=0.6,na.rm=TRUE)
      class4UL[j-9]<-quantile(tempdata,probs=0.8,na.rm=TRUE)
      class5UL[j-9]<-quantile(tempdata,probs=1,na.rm=TRUE)
      tempboundvals<-quantile(tempdata,probs=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE)
      class1data<-tempdata[(tempdata>=tempboundvals[1]) & (tempdata<=tempboundvals[2])]
      class2data<-tempdata[(tempdata>tempboundvals[2]) & (tempdata<=tempboundvals[3])]
      class3data<-tempdata[(tempdata>tempboundvals[3]) & (tempdata<=tempboundvals[4])]
      class4data<-tempdata[(tempdata>tempboundvals[4]) & (tempdata<=tempboundvals[5])]
      class5data<-tempdata[(tempdata>tempboundvals[5]) & (tempdata<=tempboundvals[6])]
      class1mean[j-9]<-mean(class1data,na.rm=TRUE)
      class2mean[j-9]<-mean(class2data,na.rm=TRUE)
      class3mean[j-9]<-mean(class3data,na.rm=TRUE)
      class4mean[j-9]<-mean(class4data,na.rm=TRUE)
      class5mean[j-9]<-mean(class5data,na.rm=TRUE)
    }
    
    plot(x=(10:nGOIvals),y=class1UL,type="l",
         ylim=c(min(rbind(class1UL,class2UL,class3UL,class4UL,class5UL)),
                max(rbind(class1UL,class2UL,class3UL,class4UL,class5UL))),
         ylab="Class upper limits",
         xlab="Number of observations",
         main="Convergence plot"
    )
    lines(x=(10:nGOIvals),y=class2UL,type="l")
    lines(x=(10:nGOIvals),y=class3UL,type="l")
    lines(x=(10:nGOIvals),y=class4UL,type="l")
    lines(x=(10:nGOIvals),y=class5UL,type="l")
    
    plot(x=(10:nGOIvals),y=class1mean,type="l",
         ylim=c(min(rbind(class1mean,class2mean,class3mean,class4mean,class5mean)),
                max(rbind(class1mean,class2mean,class3mean,class4mean,class5mean))),
         ylab="Class means",
         xlab="Number of observations",
         main="Convergence plot"
    )
    lines(x=(10:nGOIvals),y=class2mean,type="l")
    lines(x=(10:nGOIvals),y=class3mean,type="l")
    lines(x=(10:nGOIvals),y=class4mean,type="l")
    lines(x=(10:nGOIvals),y=class5mean,type="l")
  }
  
  par(mfrow=c(1,1))
  return(boundvals)
}
