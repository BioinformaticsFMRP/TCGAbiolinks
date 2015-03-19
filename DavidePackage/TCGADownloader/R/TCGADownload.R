#'@title TCGA Download
#'
#'@description 
#'  Download data previously selected using the TCGAQuery function
#'TCGADownload(fileURL, earlyStop = 0)
#'
#'
#'@param fileURL The character list of URLs referring to the exact location of the data queried
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

TCGADownload <-function(fileURL, earlyStop = 0){
  mainDir <- createDir("data/final")
  if(earlyStop==0) earlyStop <- length(fileURL)
  for(k in 1:earlyStop){
    dir<-paste(mainDir, strsplit(fileURL[k], split='/', fixed=TRUE)[[1]][14],sep="/")
    dir.create(dir, showWarnings = F)
    download(fileURL[k],
                  destfile = paste(dir,
                                   strsplit(fileURL[k], split='/', fixed=TRUE)[[1]][15],
                                   sep="/"),
                  mode="w",
                  quiet = 1)
      print(paste("Downloaded",k,"out of",length(fileURL),sep=" "))
  }
}
