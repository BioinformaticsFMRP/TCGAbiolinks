#' @title TCGA query
#'
#' @description  Crossfinding of file locations for downloading (TCGADownload)
#' TCGAQuery(tumor = "all",centerType = "all",center = "all",
#' platform = "all",level = "all",version = "all",i = F,file = "data/dataFolders.rda",
#' qOutput = "data/query/")
#' @param tumor tumor code
#' @param centerType type code
#' @param center center code
#' @param platform platform code
#' @param level level 1 2 3
#' @param version version code -TO DO-
#' @param i - interactive -TO DO-
#' @param file - link reference data matrix
#' @param qOutput place where the query is saved to be downloaded automatically. 
#'        The folder can be specified in both TCGAQuery and TCGADownload
#' 
#' @author Davide
#' 
#' @seealso TCGADownload
#' @export
#' @import downloader

TCGAQuery <- function(tumor = "all",
                      centerType = "all",
                      center = "all",
                      platform = "all",
                      level = "all",
                      version = "all",
                      i = F,
                      file = "data/dataFolders.rda",
                      qOutput = "data/query/"){
  load(file = file) #please add a way to access the package data storing
  if(!i){
    ifelse(tumor != "all",x<-subset(dataFolders, dataFolders[,"Tumor"] == tolower(tumor)),x<-dataFolders)

    if(centerType != "all" && is.null(nrow(x))) x<-subset(x, x["CenterType"] == tolower(centerType))
    if(centerType != "all" && !is.null(nrow(x))) x<-subset(x, x[,"CenterType"] == tolower(centerType))

    if(center != "all" && is.null(nrow(x))) x<-subset(x, x["Center"] == tolower(center))
    if(center != "all" && !is.null(nrow(x))) x<-subset(x, x[,"Center"] == tolower(center))
    

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

    queryURI = NULL
    dir.create(path = qOutput, showWarnings = F)
    if(is.null(nrow(x))){
      download(x["Manifest"],
                    destfile = paste(qOutput,"/filenames.txt",sep=""),
                    mode="w",
                    quiet = 1)
      queryURI <- paste(unlist(strsplit(x["Manifest"], split='MANIFEST.txt', fixed=TRUE)),
                   as.character(read.table(file = paste(qOutput,"/filenames.txt",sep=""))[2]$V2),sep="")
      print("Donwloaded.")
    }else{
      for(j in 1:length(x[,"Tumor"])){
        download(x[,"Manifest"][j],
                      destfile = paste(qOutput,"/filenames.txt",sep=""),
                      mode="w",
                      quiet = 1) #character. The mode with which to write the file.
                                 #Useful values are "w", "wb" (binary), "a" (append) and "ab".
                                 #Only used for the "internal" method.
                      #APPEND IS NOT WORKING

        print(paste("Downloaded:",j,"out of",length(x[,"Tumor"]),sep=" "))
        queryURI<-c(queryURI,paste(unlist(strsplit(x[,"Manifest"][j], split='MANIFEST.txt', fixed=TRUE)),
              as.character(read.table(file = paste(qOutput,"/filenames.txt",sep=""))[2]$V2),sep=""))
        unlink(paste(qOutput,"/filenames.txt",sep=""))
      }
    }
  }
  print(paste("We found",length(queryURI),"files",sep=" "))
  save(queryURI, file = paste(qOutput,"/fileURLs.rda",sep=""))
  #todo - add the showing of the result in a human readable way
  #return(queryURI)
}

