# FUNCTION TO PERFORM HARBIN QUANTIFICATION

harbin.quant<-function()
{
  require(tcltk)
  require(stringr)
  
  happy<-FALSE
  goahead<-FALSE
  
  while(!happy){
    
    # Open gene of interest files
    tkmessageBox(message="Select the gene-of-interest file(s)")
    GOIfiles<-open.files()
    cat("Gene-of-interest file(s)\n")
    for(j in 1:length(GOIfiles)){print(GOIfiles[[j]])}
    cat("\n")
    
    # Specify number of reference genes
    in.interval<-FALSE
    while(!in.interval){
      cat("How many reference genes do you want to use? (Enter an integer between 1 and 6): ")
      refgenenum.answer <- readline()
      if(!(refgenenum.answer %in% c(1:6))){
        cat("\nYou should enter an integer between 1 and 6!")
      }
      else{
        in.interval<-TRUE
        cat("\n")
      }
    }
        
    # Open reference gene files
    refgenefiles<-vector("list",refgenenum.answer)
    for(i in 1:refgenenum.answer){
      tkmessageBox(message=paste("Select the Reference gene ",i," file(s)",sep=""))
      refgenefiles[[i]]<-open.files()
      cat(paste("Reference gene ",i," file(s)\n",sep=""))
      for(j in 1:length(refgenefiles[[i]])){print(refgenefiles[[i]][[j]])}
      cat("\n")
    }
    cat("Go ahead with the selected files? (yes/no/cancel): ")
    goahead.answer <- readline()
    cat("\n")
    if(goahead.answer %in% c("y","Y","yes","Yes","YES")){
      happy<-TRUE
      goahead<-TRUE
    }
    if(goahead.answer %in% c("c","C","cancel","Cancel","CANCEL")){
      happy<-TRUE
      goahead<-FALSE
    }    
  }
  
  if(goahead){
    refdatabase<-NULL
    # Select reference database file (if applicable)
    cat("Do you want to use a reference database file? (yes/no/cancel): ")
    refdatabase.answer <- readline()
    cat("\n")
    if(refdatabase.answer %in% c("y","Y","yes","Yes","YES")){
      tkmessageBox(message="Select the Reference database file")
      refdatabase<-open.files()[[1]]
    }
    if(refdatabase.answer %in% c("c","C","cancel","Cancel","CANCEL")){
      goahead<-FALSE
    }  
  }
    
  if(goahead){    
    # Gene of interest normalisation
    if(is.null(refdatabase)){
      GOI.normalise(GOIfiles,refgenefiles,refdatabase=NULL,write.output=TRUE)
    } else {
      GOI.normalise(GOIfiles,refgenefiles,refdatabase=refdatabase,write.output=TRUE)
    }
  }
}
