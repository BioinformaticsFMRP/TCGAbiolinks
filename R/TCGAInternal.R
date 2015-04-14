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
#'
#' @description Internal use only.
#'
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


#
# @title DownloadHTML
# @description DownloadHTML content.
# @param url url path
# @keywords internal
# @import RCurl httr
#




#httr STILL WORKING BADLY -----> I ENCAPSULED THE FUNCTIONS but it still slower than the new function when httr is running with many errors
# DownloadHTML<- function(url,secure = F){
#
#     bo2 = T
#     count <- 0
#     handle_find(url)
#     while(bo2){
#       request = try(GET(url, timeout(100)), silent = T)
#       if( class(request) == "try-error"){
#         Sys.sleep(1)
#         bo2 = T
#         count = count + 1
#         handle_find(url)
#         if(count%%15==0) print(paste("Reconnection attempt #",count,sep=""))
#       }else{
#         bo2 = F
#       }
#       if(count>=31){
#         message("Httr failed connecting. If you are downloading very big files (proteins for example) take a look to the documentation for that.")
#         return(DownloadHTMLs(url))
#       }
#     }
#     u<-read.table(textConnection(content(request, as = 'text')), sep = ",", header = T)
#
#     return(deparse(u))
# }

#' @title DownloadHTML
#' @description DownloadHTML content.
#' @param url url path
#' @keywords internal
#' @import RCurl XML downloader


#Slow secure function
.DownloadHTML<- function(url){
  if(RCurl::url.exists(url)){
    download(url,
             "temp.html",
             mode="wb",
             quiet = 1)
    tmp <- htmlTreeParse("temp.html")
    return(capture.output(tmp))
  }else{
    stop("Can't find URL. Please check the web site or the internet connection.")
  }
}


#NEW SOLUTION --- TOO SLOW --- cleaner --- modification needed everywhere not optimal the above is better
# DownloadHTML<- function(url){
#   if(RCurl::url.exists(url)){
#       download(url,
#                "temp.html",
#                mode="wb",
#                quiet = 1)
#       tmp <- htmlTreeParse("temp.html")
#       tmp <- xmlChildren(xmlChildren(xmlChildren(xmlRoot(tmp))$body)$pre)
#       u <- NULL
#       for(i in 1:length(tmp)){
#         if(!is.na(xmlValue(tmp[i]$a)) &&
#              xmlValue(tmp[i]$a) !="Name" &&
#              xmlValue(tmp[i]$a) !="Last modified" &&
#              xmlValue(tmp[i]$a) !="Size" &&
#              xmlValue(tmp[i]$a) !="Parent Directory") u<-c(u, xmlValue(tmp[i]$a))
#       }
#       rm(tmp)
#       unlink("temp.html")
#       return(u)
#   }else{
#     stop("Can't find URL. Please check the web site or the internet connection.")
#   }
# }

GrepSite <- function(x,Key){
  x <- x[grep(Key, x)]
  x = sapply(strsplit(x, ">"), function(y) y[2])
  x = sapply(strsplit(x, "<"), function(y) y[1])
  x <- x[grep("/", x)]
  if(length(grep("lost",x))!=0) x <- x[-grep("lost", x)]
  return(x)
}





# WHAT FOR?
# .DownloadManifest <- function(siteNewLevel){
#   site3 <- paste(siteNewLevel, "MANIFEST.txt",sep="")
#   x <- .DownloadURL(site3)
#   writeLines(x, "x2.txt" )
#   tmp2 <- read.table("x2.txt", quote="\"", stringsAsFactors = F)[,2]
#   tmp2 <- tmp2[ nchar(tmp2) > 20 ]
#   return(tmp2)
# }
#
# .DownloadSdrf <- function(siteNewLevel){
#   x <- .DownloadURL(siteNewLevel)
#   x2 <- x[grep("sdrf",x)]
#   x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
#   x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
#   site3 <- paste(siteNewLevel, x2,sep="" )
#
#   x <- .DownloadURL(site3)
#   writeLines(x, "x2.txt" )
#   tmp2 <- as.data.frame(read.delim("x2.txt",stringsAsFactors = F))
#   return(tmp2)
# }
#
# .DownloadTypeFile <- function(siteNewLevel,keyDown){
#   x <- .DownloadURL(siteNewLevel)
#   x2 <- x[grep(keyDown,x)]
#   x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
#   x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
#   site3 <- paste(siteNewLevel, x2,sep="" )
#
#   x <- .DownloadURL(site3)
#   writeLines(x, "x2.txt" )
#   tmp2 <- as.data.frame(read.delim("x2.txt",stringsAsFactors = F))
#   return(tmp2)
# }
#
# .DownloaDmageTAB <- function(Description,TumorDataList, keySpecies,startK, stopK, typeProtein = F ){
#   Description2 <- paste(Description, keySpecies, sep = "")
#   Description_i_ord <- paste(Description2, "?C=M;O=D", sep = "")
#   x <- .DownloadURL(Description_i_ord)
#   if(length(x)!=10){
#     siteNewLevel <- .FindGrepSite(x,Key="mage-tab",Description2)
#     siteNewLevelSdrf <- DownloadSdrf(siteNewLevel)
#     tmp2 <-  siteNewLevelSdrf$Comment..TCGA.Barcode.
#
#     if(typeProtein==T){
#       siteNewLevelDesign <- DownloadTypeFile(siteNewLevel,"design")
#       tmp2 <- siteNewLevelDesign$Sample.description
#     }
#
#     tmp2 <- tmp2[grep("TCGA",tmp2)]
#     NumberSample <- length(unique(substr(tmp2, startK, stopK)))
#     msgOUT <-  paste(Type, SpecieCurr, CenterCurr, unique(tmp4$Platform) , ".n samples", NumberSample, sep=" ")
#     print(msgOUT)
#     SampleTmp <- unique(substr(tmp2, startK, stopK))
#     idx<- which(names(TumorDataList) == unique(tmp4$Platform))
#     TumorDataList[[idx]] <- SampleTmp
#   }
#   return(TumorDataList)
# }
#
# .DownloaDmageTAB_sdrf <- function(Description,keySpecies,KeyGrep1 = "mage-tab", KeyGrep2 = "sdrf"){
#   #Description2 <- paste(Description, keySpecies, sep = "")
#   Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
#   x <- .DownloadURL(Description_i_ord)
#   print(keySpecies)
#   print(Description)
#
#   if(length(x)!=10){
#     siteNewLevel <- .FindGrepSite(x,Key=KeyGrep1,Description)
#     print(siteNewLevel)
#     x <- .DownloadURL(x)
#     #print(x)
#     #print(KeyGrep2)
#     x2 <- x[grep(KeyGrep2,x)]
#     #print(x2)
#     x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
#     x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
#     site3 <- paste(siteNewLevel, x2,sep="" )
#     #site4 <- paste(keySpecies,unlist(strsplit(site3,keySpecies))[2],sep="")
#     site4 <- unlist(strsplit(site3,keySpecies))[2]
#
#     print(site4)
#     return(site4)
#   } else{return("")}
#
# }
#
# #' @title DownloadData_fromURL
# #' @description Download HTML content.
# #' @param url url path
# #' @keywords internal
# #' @import httr
# .DownloadData_fromURL <- function(url, sep = ",", header = TRUE){
#   handle_find(url)
#   handle_reset(url)
#   request <- GET(url)
#   stop_for_status(request)
#   handle <- textConnection(content(request, as = 'text'))
#   on.exit(close(handle))
#   read.table(handle, sep = sep, header = header)
# }
