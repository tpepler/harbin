# Function to specify single/multiple files for opening

open.files<-function()
{
  fileName <- tclvalue(tkgetOpenFile(multiple=T))
  if (!nchar(fileName)) {
    tkmessageBox(message = "No file was selected!")
    return(NULL)
  }
  endpoints<-str_locate_all(string=fileName,pattern="\\}")[[1]]
  nfiles<-nrow(endpoints)
  filelist<-vector("list",nfiles)
  for(i in 1:nfiles){
    if(i==1){
      filelist[i]<-substr(fileName,2,endpoints[i,"end"]-1)
    }
    if(i>1){
      filelist[i]<-substr(fileName,endpoints[(i-1),"end"]+3,endpoints[i,"end"]-1)
    }
  }
  return(filelist)
}