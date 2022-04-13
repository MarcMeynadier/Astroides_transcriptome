# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
dataPath<-'../data/net/7_functionnalAnnotation'
wdPath<-paste(scriptPath,dataPath,sep='/')
setwd(wdPath)

file = 