#' @title TCGA query
#'
#' @description  TCGA Query
#'
#' @param Tumor -  tumor code
#' @param centerType 
#' @param center
#' @param platform
#' @param level
#' @param version
#' @param i - interactive
#' @param file - link reference data matrix
#' @export

TCGAQuery <- function(tumor = "all", 
                      centerType = "all", 
                      center = "all", 
                      platform = "all",
                      level = "all",
                      version = "all",
                      i = F,
                      file = "data/PlatformMat.rda"){
  load(file = file) #please add a way to access the package data storing
  if(!i){
    ifelse(tumor != "all",x<-subset(PlatformMat, PlatformMat[,"Tumor"] == tolower(tumor)),x<-PlatformMat)
    
    if(centerType != "all" && is.null(nrow(x))) x<-subset(x, x["CenterType"] == tolower(centerType))
    if(centerType != "all" && !is.null(nrow(x))) x<-subset(x, x[,"CenterType"] == tolower(centerType))
   
    
    if(platform != "all" && is.null(nrow(x))) x<-subset(x, x["Platform"] == tolower(platform))
    if(platform != "all" && !is.null(nrow(x))) x<-subset(x, x[,"Platform"] == tolower(platform))
    
    
    if(level != "all"){
      if(is.null(nrow(x))){
        l <- grep(paste("Level_",as.character(level),sep=""),x["Folder"])
        if(length(l)==0) x<-NULL
      }else{
        l <- grep(paste("Level_",as.character(level),sep=""),x[,"Folder"])
        x<-x[l,]
      }
    }
#     if(version != "all"){ #something smarter could be done (maybe "lastUp" as keyword)
#       if(is.null(nrow(x)))l <- grep(paste("Level_",as.character(version),sep=""),x["Folder"])
#       if(!is.null(nrow(x)))l <- grep(paste("Level_",as.character(version),sep=""),x[,"Folder"])
#       if(length(l)>1) x<-x[l,]
#       if(length(l)==0)x<-NULL
#     }
    
    if(is.null(x)){
      return("Nothing found. Check the proper spelling in the documentation.")
    }else if(is.null(nrow(x))) {
      print("Found: 1 folder. Start downloading filenames:")
    }else{
      print(paste("Found:", length(x[,1]), "folders. Start downloading filenames:",sep=" "))
    }
  
    ret = NULL
    dataDir<-createDir("data/query")
    if(is.null(nrow(x))){
      download(x["Manifest"],
                    destfile = paste(dataDir,"/filenames.txt",sep=""),
                    mode="a",
                    quiet = 0)
      ret <- paste(unlist(strsplit(x["Manifest"], split='MANIFEST.txt', fixed=TRUE)),
                   as.character(read.table(file = paste(dataDir,"/filenames.txt",sep=""))[2]$V2),sep="")
      print("Donwloaded.")
    }else{
      for(j in 1:length(x[,"Tumor"])){
        download(x[,"Manifest"][j],
                      destfile = paste(dataDir,"/filenames.txt",sep=""),
                      mode="a",
                      quiet = 1) #character. The mode with which to write the file. 
                                 #Useful values are "w", "wb" (binary), "a" (append) and "ab". 
                                 #Only used for the "internal" method.
                      #APPEND IS NOT WORKING
        if(as.integer(j*100/(length(x[,"Tumor"])))%%5 == 0) 
                      print(paste("Downloaded: ",as.integer(j*100/(length(x[,"Tumor"]))),"%",sep=""))
        
        ret<- c(ret,paste(unlist(strsplit(x[,"Manifest"][j], split='MANIFEST.txt', fixed=TRUE)),
                     as.character(read.table(file = paste(dataDir,"/filenames.txt",sep=""))[2]$V2),sep=""))
      }
    }
  }
  return(ret)
}

