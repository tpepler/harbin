#**********************************************************************************
# Function to normalise the gene of interest with the geometric means of the reference gene(s)

GOI.normalise<-function(GOIfiles, refgenefiles, refdatabase=NULL, write.output=FALSE, 
                        pathname='', filename='GOI_output.csv', graph=c("density","histogram")){
  # GOIfiles: list of gene of interest file names
  # refgenefiles: list of file names for one or more reference genes
  # refdatabase: file name of reference database
  # write.output: logical, indicating whether GOI output should be written to a CSV file
  # pathname: location of reference database file; CSV output file (if applicable) will also be written to this location
  # filename: name of GOI CSV output file to be written (if applicable)
  # graph: type of graph to be used for represenation of distribution of normalised GOI values
  
  # Determine number of reference genes
  nrefgenes<-length(refgenefiles)
  
  # Read in gene of interest files
  tempout<-read.qPCRfiles(GOIfiles)
  GOIdata<-tempout$Unknown[order(tempout$Unknown[,"Name"]),]
  n<-nrow(GOIdata)

  # Determine which of Rep.Calc.Conc values are invalid
  GOI.min<-tempout$Min
  GOI.max<-tempout$Max
  valid.ind<-rep(0,times=n)
  valid.ind[GOIdata[,"Rep.Calc.Conc"]<GOI.min]<- -1
  valid.ind[GOIdata[,"Rep.Calc.Conc"]>GOI.max]<- 1
  valid.ind[is.na(GOIdata[,"Rep.Calc.Conc"])]<- 9
  
  # Read in reference gene files and determine validity (i.e. whether the rows correspond with gene of interest set)
  refgenedata<-vector("list",nrefgenes)
  refgenevals<-vector("list",nrefgenes)
  for(j in 1:nrefgenes){
    tempout<-read.qPCRfiles(refgenefiles[[j]])
    refgenedata[[j]]<-tempout$Unknown[order(tempout$Unknown[,"Name"]),]
    for(i in 1:n){
      if(GOIdata[i,"Name"]!=refgenedata[[j]][i,"Name"]){
        if(!(GOIdata[i,"Name"] %in% refgenedata[[j]][,"Name"])){
          stop(paste("\nNo concentration value for ",GOIdata[i,"Name"]," (from gene of interest data) found in Reference gene ",j," file(s)!\n",sep=""))
        }
        if(!(refgenedata[[j]][i,"Name"] %in% GOIdata[,"Name"])){
          stop(paste("\nNo concentration value for ",refgenedata[[j]][i,"Name"]," (from Reference gene ",j," data) found in gene of interest file(s)!\n",sep=""))
        }
        #print(i)
        #print(GOIdata[i,"Name"])
        #print(refgenedata[[j]][i,"Name"])
        #stop(paste("\nWarning: Reference gene ",j," names do not match gene of interest names!\n",sep=""))
      }
    }
    refgenevals[[j]]<-as.numeric(refgenedata[[j]][,"Rep.Calc.Conc"])
  }

  # Normalisation of gene of interest Rep.Calc.Conc values
  refindex<-(Reduce("*",refgenevals))^(1/nrefgenes) # geometric mean of reference gene values
  GOI.normdata<-GOIdata[,"Rep.Calc.Conc"]/refindex
  
  # Select valid gene of interest data only
  GOI.validdata<-GOI.normdata[(valid.ind==0)]
  
  refbaseadd.answer<-"n"
  if(!is.null(refdatabase)){
    
    # Comparison to reference data base
    refdatabase.data<-read.csv(paste(pathname,refdatabase,sep=''),head=TRUE)
    for(i in 1:n){
      if(GOIdata[i,"Name"] %in% refdatabase.data[,"Name"]){
        cat(paste("Warning: GOI name ",GOIdata[i,"Name"]," already found in reference database!\n",sep=""))
      }
    }
    boundvals<-convergence.check(c(refdatabase.data[,"GOI.normalised"],GOI.validdata),graph=graph)
    refdatabase.index<-rep(0,times=nrow(refdatabase.data))
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>=0) & (refdatabase.data[,"GOI.normalised"]<=boundvals[1]))]<-1
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>boundvals[1]) & (refdatabase.data[,"GOI.normalised"]<=boundvals[2]))]<-2
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>boundvals[2]) & (refdatabase.data[,"GOI.normalised"]<=boundvals[3]))]<-3
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>boundvals[3]) & (refdatabase.data[,"GOI.normalised"]<=boundvals[4]))]<-4
    refdatabase.index[refdatabase.data[,"GOI.normalised"]>boundvals[4]]<-5
    refdatabase.data[,"Interval"]<-refdatabase.index
    
    # Harbin test: Reference data base vs. new data
    harbin.out<-harbin.test(x=refdatabase.data[,"GOI.normalised"],y=GOI.validdata,reps=1000)
    cat(paste("\nProportion of labels changing in reference data base: ",round(harbin.out$statistic*100,1),"%\n",sep=""))
    cat("H0: New data originated from same distribution as reference data\n")
    cat("H1: New data and reference data come from different distributions\n")
    cat(paste("Harbin test p-value = ",harbin.out$p.value,"\n",sep=""))
    if(harbin.out$p.value<=0.05){
      cat("\nWARNING: Reference data and new data may not be compatible!\n")
    }
  } else {
    
    # Convergence check using only current normalised data set
    boundvals<-convergence.check(GOI.validdata,graph=graph)
  }  

  # Calculate categories for current data set
  interval.index<-rep(0,times=n)
  interval.index[((GOI.normdata>=0) & (GOI.normdata<=boundvals[1]))]<-1
  interval.index[((GOI.normdata>boundvals[1]) & (GOI.normdata<=boundvals[2]))]<-2
  interval.index[((GOI.normdata>boundvals[2]) & (GOI.normdata<=boundvals[3]))]<-3
  interval.index[((GOI.normdata>boundvals[3]) & (GOI.normdata<=boundvals[4]))]<-4
  interval.index[GOI.normdata>boundvals[4]]<-5
  interval.index[(valid.ind==-1)]<- -999
  interval.index[(valid.ind==1)]<- 999
  resultsmat<-data.frame(Name=GOIdata[,"Name"],
                         GOI.data=GOIdata[,"Rep.Calc.Conc"],
                         Reference.index=refindex,
                         GOI.normalised=GOIdata[,"Rep.Calc.Conc"]/refindex,
                         Interval=interval.index,
                         Date=GOIdata[,"Date"],
                         Time=GOIdata[,"Time"])

  # Write output file for current data set (if applicable)
  if(write.output){
    tkmessageBox(message="Select a name for the output file (for current data)")
    outfile <- tclvalue(tkgetSaveFile())
    #write.csv(resultsmat,paste(pathname,filename,sep=''),row.names=FALSE,quote=FALSE)
    write.csv(resultsmat,outfile,row.names=FALSE,quote=FALSE)
  }
  
  # Add new data to reference data base (if applicable)
  if(!is.null(refdatabase)){
    cat("\nDo you want to add the new data to the reference database (y/n): ")
    refbaseadd.answer <- readline()
    cat("\n")
    refdbsaved.message<-"(NOT saved to reference database)"
    if(refbaseadd.answer %in% c("y","Y","yes","Yes","YES")){
      refdbsaved.message<-"(SAVED to reference database)"
      #write.csv(rbind(refdatabase.data,resultsmat),paste(pathname,refdatabase,sep=''),row.names=FALSE,quote=FALSE)
      write.csv(rbind(refdatabase.data,resultsmat),refdatabase,row.names=FALSE,quote=FALSE)
    }
    cat(paste("Gene of interest (new data) results ",refdbsaved.message,":\n\n",sep=""))
  } else {
    cat("Gene of interest results:\n\n")
  }
  
  # Output current data set results
  return(resultsmat)
}

#**********************************************************************************
