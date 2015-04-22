#' @title .onAttach
#' @description  Load required data into gloval enviroment
#' @keywords internal
.onAttach <- function (libname, pkgname){
  file = system.file("extdata/dataFolders.rda",package="TCGAbiolinks")
  if(file.exists(file)) {
    load(file,.GlobalEnv)
  }
  else{
    message("Please run TCGAUpdate() to obtain the TCGA table")
  }
  metadata = system.file("extdata/barcodes.rdata",package="TCGAbiolinks")
  if(file.exists(metadata)) load(metadata,.GlobalEnv)

}

#' @title Creates a directory
#' @description Internal use only.
#' @param base directory path
#' @keywords internal
createDir <- function(base){
  i="";
  while(file.exists(paste0(base, i))){
    if(i==""){
      i=1;
    }else{
      i=i+1;
    }
  }
  toDir = paste0(base, i)
  dir.create(toDir, showWarnings = F, recursive = T, mode = "0777")
  toDir
}

#' @title DownloadHTML
#' @description DownloadHTML content.
#' @param url url path
#' @keywords internal
#' @import XML downloader
DownloadHTML<- function(url){
    download(url,
             "temp.html",
             mode="wb",
             quiet = 1)
    tmp <- htmlTreeParse("temp.html")
    unlink("temp.html")
    return(capture.output(tmp))
}

#' @title GrepSite
#' @description Grep the requested key in the html character sequence
#' @param siteChar, keyChar
#' @keywords internal
GrepSite <- function(x,Key){
  x <- x[grep(Key, x)]
  x = sapply(strsplit(x, ">"), function(y) y[2])
  x = sapply(strsplit(x, "<"), function(y) y[1])
  x <- x[grep("/", x)]
  if(length(grep("lost",x))!=0) x <- x[-grep("lost", x)]
  return(x)
}