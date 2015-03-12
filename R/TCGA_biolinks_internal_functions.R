.FindGrepSite <- function(x,Key,Description){
  x2 <- x[grep(Key, x)]
  
  if( Key != "sdrf"){ x2 <- x2 [- grep("tar.gz", x2)][1] }
  
  
  x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
  x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
  site2 <- paste(Description, x2,sep="" )
  return(site2)
}

.DownloadURL <- function(Site){
  opts = curlOptions(ftp.use.epsv = T,
                     dirlistonly = T,
                     ssl.verifypeer = F,
                     timeout = 50)
  handle = getCurlHandle(.opts = opts)
  Site <- URLencode(Site)
  
  bo1 = T
  count1 = 0
  while(bo1){
    bo1 = !url.exists(Site, .opts = opts)
    count1 = count1 + 1
    if(count1>1) {
      print(paste("Reonnection to the url #",count1,sep=""))
      Sys.sleep(1)
    }
    if(count1 == 20) stop(paste("The url (",Site,") is not existing or not available.
                                Please check your internet connection or contact the mantainers for an update.",sep=""))
  }
  
  if(interactive() && ("ssl" %in% names(curlVersion()$features))) {
    bo2 = T
    count <- 0
    while(bo2){
      x = try(getURLContent(Site, verbose = F, curl = handle))
      if( class(x) == "try-error"){
        Sys.sleep(1)
        bo2 = T
        count = count + 1
        if(count>=1) print(paste("Reconnection attempt #",count,sep=""))
      }else if(count==20){
        stop("Connetion limit exceded. Check your internet connection and your proxy settings.")
      }else{
        bo2 = F
      }
    }
  }else{
    stop("Curl is not configured for ssl.")
  }
  
  x <- unlist(strsplit(x, "\n"))
  return(x)
  }

.DownloadManifest <- function(siteNewLevel){
  site3 <- paste(siteNewLevel, "MANIFEST.txt",sep="")
  x <- .DownloadURL(site3)
  writeLines(x, "x2.txt" )
  tmp2 <- read.table("x2.txt", quote="\"", stringsAsFactors = F)[,2]
  tmp2 <- tmp2[ nchar(tmp2) > 20 ]
  return(tmp2)
}

.DownloadSdrf <- function(siteNewLevel){
  x <- .DownloadURL(siteNewLevel)
  x2 <- x[grep("sdrf",x)]
  x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
  x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
  site3 <- paste(siteNewLevel, x2,sep="" )
  
  x <- .DownloadURL(site3)
  writeLines(x, "x2.txt" )
  tmp2 <- as.data.frame(read.delim("x2.txt",stringsAsFactors = F))
  return(tmp2)
}

.DownloadTypeFile <- function(siteNewLevel,keyDown){
  x <- .DownloadURL(siteNewLevel)
  x2 <- x[grep(keyDown,x)]
  x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
  x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
  site3 <- paste(siteNewLevel, x2,sep="" )
  
  x <- .DownloadURL(site3)
  writeLines(x, "x2.txt" )
  tmp2 <- as.data.frame(read.delim("x2.txt",stringsAsFactors = F))
  return(tmp2)
}

.DownloaDmageTAB <- function(Description,TumorDataList, keySpecies,startK, stopK, typeProtein = F ){
  Description2 <- paste(Description, keySpecies, sep = "")
  Description_i_ord <- paste(Description2, "?C=M;O=D", sep = "")
  x <- .DownloadURL(Description_i_ord)
  if(length(x)!=10){
    siteNewLevel <- .FindGrepSite(x,Key="mage-tab",Description2)
    siteNewLevelSdrf <- DownloadSdrf(siteNewLevel)
    tmp2 <-  siteNewLevelSdrf$Comment..TCGA.Barcode.
    
    if(typeProtein==T){
      siteNewLevelDesign <- DownloadTypeFile(siteNewLevel,"design")
      tmp2 <- siteNewLevelDesign$Sample.description
    }
    
    tmp2 <- tmp2[grep("TCGA",tmp2)]
    NumberSample <- length(unique(substr(tmp2, startK, stopK)))
    msgOUT <-  paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,  " .n samples ", NumberSample, sep="")
    print(msgOUT)
    SampleTmp <- unique(substr(tmp2, startK, stopK))
    idx<- which(names(TumorDataList) == unique(tmp4$Platform))
    TumorDataList[[idx]] <- SampleTmp
  }
  return(TumorDataList)
}

.DownloaDmageTAB_sdrf <- function(Description,keySpecies,KeyGrep1 = "mage-tab", KeyGrep2 = "sdrf"){
  Description2 <- paste(Description, keySpecies, sep = "")
  Description_i_ord <- paste(Description2, "?C=M;O=D", sep = "")
  x <- .DownloadURL(Description_i_ord)
  if(length(x)!=10){
    siteNewLevel <- .FindGrepSite(x,Key=KeyGrep1,Description2)
    x <- .DownloadURL(siteNewLevel)
    x2 <- x[grep(KeyGrep2,x)]
    
    x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
    x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
    site3 <- paste(siteNewLevel, x2,sep="" )
    #site4 <- paste(keySpecies,unlist(strsplit(site3,keySpecies))[2],sep="")
    site4 <- unlist(strsplit(site3,keySpecies))[2]
    
    print(site4)
    return(site4)
  } else{return("")}
  
}

#' @title Creates a directory
#'
#' @description Internal use only.
#'
#' @param base directory path
#' @keywords internal

.createDirectory <- function(base){
  i="";
  while(file.exists(paste(base, i, sep=""))){
    if(i==""){
      i=1;
    }else{
      i=i+1;
    }
  }
  toDir = paste(base, i, sep="")
  dir.create(toDir)
  #dir.create(base, showWarnings = FALSE, recursive = FALSE, mode = "0777")
  
  toDir
}

.DownloadData_fromURL <- function(url, sep = ",", header = TRUE){
  require(httr)
  handle_find(url)
  handle_reset(url)
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}
