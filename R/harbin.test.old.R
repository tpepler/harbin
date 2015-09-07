harbin.test.old<-function(x, y, reps=1000){
  
  x.n<-length(x)
  y.n<-length(y)
  bound<-quantile(x,probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  indvec<-rep(1,times=length(x))
  indvec[((x>bound[1]) & (x<=bound[2]))]<-2
  indvec[((x>bound[2]) & (x<=bound[3]))]<-3
  indvec[((x>bound[3]) & (x<=bound[4]))]<-4
  indvec[x>bound[4]]<-5
  
  changeprop<-rep(NA,times=reps)
  for(i in 1:reps){
    bootsamp<-x[sample(1:x.n,size=y.n,replace=TRUE)]
    bootbound<-quantile(c(x,bootsamp),probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
    bootindvec<-rep(1,times=x.n)
    bootindvec[((x>bootbound[1]) & (x<=bootbound[2]))]<-2
    bootindvec[((x>bootbound[2]) & (x<=bootbound[3]))]<-3
    bootindvec[((x>bootbound[3]) & (x<=bootbound[4]))]<-4
    bootindvec[x>bootbound[4]]<-5
    changevec<-bootindvec-indvec
    changeprop[i]<-(x.n-length(changevec[changevec==0]))/x.n
  }
  
  newbound<-quantile(c(x,y),probs=c(0.2,0.4,0.6,0.8),na.rm=TRUE)
  newindvec<-rep(1,times=x.n)
  newindvec[((x>newbound[1]) & (x<=newbound[2]))]<-2
  newindvec[((x>newbound[2]) & (x<=newbound[3]))]<-3
  newindvec[((x>newbound[3]) & (x<=newbound[4]))]<-4
  newindvec[x>newbound[4]]<-5
  newchangevec<-newindvec-indvec
  newchangeprop<-(x.n-length(newchangevec[newchangevec==0]))/x.n
  p.value<-length(changeprop[changeprop>=newchangeprop])/reps # p-value for null hypothesis of no change
  
  return(list(statistic=newchangeprop,crit.val=quantile(changeprop,prob=0.95,na.rm=TRUE),p.value=p.value))
}

#********************************************************************
