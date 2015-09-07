harbin.test<-function(x, y, reps=1000){
  
  x.n<-length(x)
  y.n<-length(y)
  all.n<-x.n+y.n
  x.bound<-quantile(x,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  y.bound<-quantile(y,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  all.bound<-quantile(c(x,y),probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  
  changeprop.h0<-rep(NA,times=reps)
  for(i in 1:reps){
    #all.h0<-c(x,y)[sample(1:all.n,size=all.n,replace=FALSE)]
    all.h0<-c(x,y)[sample(1:all.n,size=all.n,replace=TRUE)]
    x.h0<-all.h0[1:x.n]
    y.h0<-all.h0[(x.n+1):all.n]
    
    x.bound.h0<-quantile(x.h0,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
    y.bound.h0<-quantile(y.h0,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
    
    x.ind.h0<-rep(1,times=x.n)
    x.ind.h0[((x.h0>x.bound.h0[1]) & (x.h0<=x.bound.h0[2]))]<-2
    x.ind.h0[((x.h0>x.bound.h0[2]) & (x.h0<=x.bound.h0[3]))]<-3
    x.ind.h0[((x.h0>x.bound.h0[3]) & (x.h0<=x.bound.h0[4]))]<-4
    x.ind.h0[x.h0>x.bound.h0[4]]<-5

    y.ind.h0<-rep(1,times=y.n)
    y.ind.h0[((y.h0>y.bound.h0[1]) & (y.h0<=y.bound.h0[2]))]<-2
    y.ind.h0[((y.h0>y.bound.h0[2]) & (y.h0<=y.bound.h0[3]))]<-3
    y.ind.h0[((y.h0>y.bound.h0[3]) & (y.h0<=y.bound.h0[4]))]<-4
    y.ind.h0[y.h0>y.bound.h0[4]]<-5

    all.ind.h0<-rep(1,times=all.n)
    all.ind.h0[((all.h0>all.bound[1]) & (all.h0<=all.bound[2]))]<-2
    all.ind.h0[((all.h0>all.bound[2]) & (all.h0<=all.bound[3]))]<-3
    all.ind.h0[((all.h0>all.bound[3]) & (all.h0<=all.bound[4]))]<-4
    all.ind.h0[all.h0>all.bound[4]]<-5

    changevec<-c(x.ind.h0,y.ind.h0)-all.ind.h0
    changeprop.h0[i]<-length(changevec[changevec!=0])/all.n
  }
  
  x.ind<-rep(1,times=x.n)
  x.ind[((x>x.bound[1]) & (x<=x.bound[2]))]<-2
  x.ind[((x>x.bound[2]) & (x<=x.bound[3]))]<-3
  x.ind[((x>x.bound[3]) & (x<=x.bound[4]))]<-4
  x.ind[x>x.bound[4]]<-5
  
  y.ind<-rep(1,times=y.n)
  y.ind[((y>y.bound[1]) & (y<=y.bound[2]))]<-2
  y.ind[((y>y.bound[2]) & (y<=y.bound[3]))]<-3
  y.ind[((y>y.bound[3]) & (y<=y.bound[4]))]<-4
  y.ind[y>y.bound[4]]<-5
  
  alldata<-c(x,y)
  all.ind<-rep(1,times=all.n)
  all.ind[((alldata>all.bound[1]) & (alldata<=all.bound[2]))]<-2
  all.ind[((alldata>all.bound[2]) & (alldata<=all.bound[3]))]<-3
  all.ind[((alldata>all.bound[3]) & (alldata<=all.bound[4]))]<-4
  all.ind[alldata>all.bound[4]]<-5
  
  changevec<-c(x.ind,y.ind)-all.ind
  changeprop<-length(changevec[changevec!=0])/all.n
  p.value<-length(changeprop.h0[changeprop.h0>=changeprop])/reps # p-value for null hypothesis of no change
  
  return(list(statistic=changeprop,
              crit.val=quantile(changeprop.h0,prob=0.95,na.rm=TRUE),
              p.value=p.value))
}

#********************************************************************
