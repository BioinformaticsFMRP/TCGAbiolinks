#'@title TCGA Download
#'
#'@description 
#'  Download data previously selected using the TCGAQuery function
#'TCGADownload(dirURL="data/query/", finalDir = "data/final", earlyStop = 0)
#'
#'
#'@param dirURL the location of the query saved URL links.
#'@param finalDir location of the final data saving
#'@param earlyStop If the number of file is too large for you need you can choose a subset
#'
#'
#'@author Davide
#'
#'@seealso TCGAQuery
#'
#'
#'
#'@export
#'@import downloader

TCGADownload <-function(dirURL="data/query/", finalDir = "data/final", earlyStop = 0){
  load(paste(dirURL,"fileURLs.rda",sep=""))
  finalDir <- createDir(finalDir)
  t1 = Sys.time()
  if(earlyStop==0) earlyStop <- length(queryURI)
  for(k in 1:earlyStop){
    dir<-paste(finalDir, strsplit(queryURI[k], split='/', fixed=TRUE)[[1]][14],sep="/")
    dir.create(dir, showWarnings = F)
    download(queryURI[k],
                  destfile = paste(dir,
                                   strsplit(queryURI[k], split='/', fixed=TRUE)[[1]][15],
                                   sep="/"),
                  mode="w",
                  quiet = 1)
      print(paste("Downloaded",k,"out of",length(queryURI),sep=" "))
    if(k%%10==0) print(paste("Estimated time for the end of the process:",
                             (as.numeric(difftime(Sys.time(), t1, units="min"))/k)*earlyStop, "minutes", sep=" "))
  }
  rm(t1)
}
