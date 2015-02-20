# Function to read in qPCR files

read.qPCRfiles<-function(files){
  
  # files: list of qPCR files to be read
  
  nfiles<-length(files)
  new.unknown.frame<-NULL
  new.standard.frame<-NULL
  for(j in 1:nfiles){
    datestamp<-scan(files[[j]],what="character",nlines=1,skip=3,sep=",",quiet=TRUE)[2] # scan date stamp on 4th line
    timestamp<-scan(files[[j]],what="character",nlines=1,skip=4,sep=",",quiet=TRUE)[2] # scan time stamp on 5th line
    varnames<-scan(files[[j]],what="character",nlines=1,skip=8,sep=",",quiet=TRUE) # scan variable names on 9th line
    tempdata<-scan(files[[j]],what="character",skip=9,sep=",",quiet=TRUE) # scan all data and store in vector
    nlines<-floor(length(tempdata)/13) # calculate number of lines in the file
    datamat<-NULL
    for(i in 1:nlines){
      temprow<-tempdata[((i-1)*13+1):(i*13)] # convert the data vector to lines again
      datamat<-rbind(datamat,temprow)
    }
    colnames(datamat)<-varnames
    allnames<-unique(datamat[datamat[,"Type"]=="Unknown","Name"])
    datamat<-datamat[datamat[,"Rep. Calc. Conc."]!="",,drop=FALSE] # select rows with non-empty cells in Rep. Calc. Conc. column
    non.empty.names<-unique(datamat[datamat[,"Type"]=="Unknown","Name"])
    names.NA<-allnames[!(allnames %in% non.empty.names)]
    new.unknown<-datamat[datamat[,"Type"]=="Unknown",c("Name","Type","Rep. Calc. Conc."),drop=FALSE] # select only rows of 'Unknown' type
    colnames(new.unknown)<-c("Name","Type","Rep.Calc.Conc.")
    if(length(names.NA)>0){
      new.unknown<-rbind(new.unknown,data.frame(Name=names.NA,Type="Unknown",Rep.Calc.Conc.=NA))
    }
    new.standard<-datamat[datamat[,"Type"]=="Standard",c("Name","Type","Rep. Calc. Conc."),drop=FALSE] # select only rows of 'Standard' type
    unknownvals<-as.numeric(as.character(new.unknown[,"Rep.Calc.Conc."]))
    standardvals<-as.numeric(as.character(new.standard[,"Rep. Calc. Conc."]))
    #newnames<-as.character(newdatamat[newdatamat[,"Rep. Calc. Conc."]!="","Name"])
    if(nrow(new.unknown)>0){
      new.unknown.frame<-rbind(new.unknown.frame,data.frame(Name=new.unknown[,"Name"],Rep.Calc.Conc=unknownvals,Date=datestamp,Time=timestamp))
    }
    if(nrow(new.standard)>0){
      new.standard.frame<-rbind(new.standard.frame,data.frame(Name=new.standard[,"Name"],Rep.Calc.Conc=standardvals))
    }
    #new.unknown.frame<-rbind(new.unknown.frame,as.data.frame(new.unknown))
    #new.standard.frame<-rbind(new.standard.frame,as.data.frame(new.standard))
  }
  new.unknown.frame[,"Name"]<-as.character(new.unknown.frame[,"Name"])
  #new.standard.frame[,"Name"]<-as.character(new.standard.frame[,"Name"])
  #return(list(Unknown=new.unknown.frame,Standard=new.standard.frame))
  return(list(Unknown=new.unknown.frame,Min=min(new.standard.frame[,"Rep.Calc.Conc"]),Max=max(new.standard.frame[,"Rep.Calc.Conc"])))
}