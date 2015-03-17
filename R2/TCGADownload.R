#' @title TCGA Download
#' 
#' 

TCGADownload <-function(fileLocation, earlyStop = 0){
  mainDir <- createDir("data/final")
  if(earlyStop==0) earlyStop <- length(fileLocation)
  for(k in 1:earlyStop){
    dir<-paste(mainDir, strsplit(fileLocation[k], split='/', fixed=TRUE)[[1]][14],sep="/")
    dir.create(dir, showWarnings = F)
    download.file(fileLocation[k],
                  destfile = paste(dir,
                                   strsplit(fileLocation[k], split='/', fixed=TRUE)[[1]][15],
                                   sep="/"),
                  mode="w",
                  method = "curl",
                  quiet = 1)
    if(as.integer(k*100/(length(fileLocation)))%%5==0) 
      print(paste("Downloaded: ",as.integer(k*100/(length(fileLocation))),"%",sep=""))
  }
}