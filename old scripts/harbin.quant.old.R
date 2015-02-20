# FUNCTION TO PERFORM HARBIN QUANTIFICATION

harbin.quant<-function()
{
  require(tcltk)
  require(stringr)
  
  # Open gene of interest files
  tkmessageBox(message="Select the gene of interest file(s)")
  GOIfiles<-open.files()
  
  # Open reference gene files
  refgenefiles<-vector("list",3)
  tkmessageBox(message="Select the Reference gene 1 file(s)")
  refgenefiles[[1]]<-open.files()
  tkmessageBox(message="Select the Reference gene 2 file(s)")
  refgenefiles[[2]]<-open.files()
  tkmessageBox(message="Select the Reference gene 3 file(s)")
  refgenefiles[[3]]<-open.files()
  
  # Select reference database file (if applicable)
  tkmessageBox(message="Select the Reference database file")
  refdatabase<-open.files()[[1]]
  
  # Gene of interest normalisation
  if(is.null(refdatabase)){
    GOI.normalise(GOIfiles,refgenefiles,refdatabase=NULL,write.output=TRUE)
  } else {
    GOI.normalise(GOIfiles,refgenefiles,refdatabase=refdatabase,write.output=TRUE)
  }
}