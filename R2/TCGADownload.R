#' @title TCGA Download
#' 
#' 

TCGADownload <-function(fileLocation, earlyStop = 0){
  mainDir <- createDir("data/final")
  if(earlyStop==0) earlyStop <- length(fileLocation)
  for(k in 1:earlyStop){
    dir<-paste(mainDir, strsplit(fileLocation[k], split='/', fixed=TRUE)[[1]][14],sep="/")
    dir.create(dir, showWarnings = F)
    download(fileLocation[k],
                  destfile = paste(dir,
                                   strsplit(fileLocation[k], split='/', fixed=TRUE)[[1]][15],
                                   sep="/"),
                  mode="w",
                  quiet = 1) 
      print(paste("Downloaded ",k," out of ",length(fileLocation),sep=""))
  }
}