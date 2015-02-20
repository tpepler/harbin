require(tcltk)

acceptrefnum<-function(){
  tkdestroy(tt)
  print(as.numeric(tclvalue(SliderValue)))
  #print(tclvalue(GOIfiles))
  
  #numrefgenes<-as.numeric(tclvalue(SliderValue))
  
  #tt2<-tktoplevel()
  #tkgrid(tklabel(tt2,text="Click and browse to select the reference gene files."))
  #if(numrefgenes>=1){
  #  ref1.but <- tkbutton(tt2, text = "Ref1", command = cat("selected ref 1 gene"))
  #  tkgrid(ref1.but)    # Place the button on the window 
  #}
  #if(numrefgenes>=2){
  #  ref2.but <- tkbutton(tt2, text = "Ref2", command = cat("selected ref 2 gene"))
  #  tkgrid(ref2.but)    # Place the button on the window 
  #}
  #if(numrefgenes>=3){
  #  ref3.but <- tkbutton(tt2, text = "Ref3", command = cat("selected ref 3 gene"))
  #  tkgrid(ref3.but)    # Place the button on the window 
  #}
  #tkfocus(tt2)
}

openGOI<-function(){
  #tclvalue(GOIfiles)<-"GOI selected"
  GOIfiles<-open.files()
  write(GOIfiles[[1]],file=paste(GOIfiles[[1]],".ref.txt",sep=""),append=F)
  if(length(GOIfiles)>1){
    for(i in 2:length(GOIfiles)){
      write(GOIfiles[[i]],file=paste(GOIfiles[[1]],".ref.txt",sep=""),append=T)    
    }
  }
}

topfunction1<-function(){
  tt <- tktoplevel()
  GOIfiles<-tclVar(init = "GOI not selected")
  GOI.but <- tkbutton(tt, text = "Select gene of interest files", command = openGOI)
  tkgrid(GOI.but)
  
  SliderValue <- tclVar("3")
  SliderValueLabel <- tklabel(tt,text=as.character(tclvalue(SliderValue)))
  tkgrid(tklabel(tt,text="Select number of reference genes: "),SliderValueLabel)
  tkconfigure(SliderValueLabel,textvariable=SliderValue)
  slider <- tkscale(tt, from=1, to=6,
                    showvalue=F, variable=SliderValue,
                    resolution=1, orient="horizontal")
  tkgrid(slider)
  OK.but <- tkbutton(tt, text = "OK", command = acceptrefnum)
  tkgrid(OK.but)    # Place the button on the window
  tkfocus(tt)
  numrefgenes<-as.numeric(tclvalue(SliderValue))
  return(numrefgenes)
}

topfunction2<-function(numrefgenes){
  tt2<-tktoplevel()
  tkgrid(tklabel(tt2,text="Click and browse to select the reference gene files."))
  if(numrefgenes>=1){
  ref1.but <- tkbutton(tt2, text = "Ref1", command = openGOI)
  tkgrid(ref1.but)    # Place the button on the window 
  }
  if(numrefgenes>=2){
  ref2.but <- tkbutton(tt2, text = "Ref2", command = openGOI)
  tkgrid(ref2.but)    # Place the button on the window 
  }
  if(numrefgenes>=3){
  ref3.but <- tkbutton(tt2, text = "Ref3", command = openGOI)
  tkgrid(ref3.but)    # Place the button on the window 
  }
  tkfocus(tt2)
}

tempval<-topfunction1()
topfunction2(tempval)
  


