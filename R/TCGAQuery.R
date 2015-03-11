#' @title TCGA Query
#'
#' @description  TCGA Query
#'
#' @param Tumor -  tumor code
#' @param siteTCGA - Plataform
#' @export


TCGAQuery <- function(Tumor,sdrfFolder,PlatformAndAssociatedData){
  message(Tumor)
  tmplf <- PlatformAndAssociatedData[PlatformAndAssociatedData$Tumor %in% tolower(Tumor),]
  TumorDataList <- vector("list", length(unique(tmplf$Platform)))
  names(TumorDataList) <- unique(tmplf$Platform)
  FolderWd <- getwd()
  setwd(sdrfFolder)
  
  for( i in 1: nrow(tmplf)){
    startK <- 1
    stopK <- 16
    if(tmplf[i,"FileName"] == "") next
    currentTFile <- gsub("/","_",tmplf[i,"FileName"])
    
    
    currentTFile<-paste(sdrfFolder,currentTFile,sep="")
    tmp2 <- read.delim(currentTFile,stringsAsFactors = F,skip =  tmplf[i,"KeySkip"], header =T )
    
    if( length(intersect(colnames(tmp2), tmplf[i,"nameCol"]))==1){ tmp2 <- tmp2[, tmplf[i,"nameCol"]  ] }
    if( length(intersect(colnames(tmp2), "Sample.description"))==1){ tmp2 <- tmp2[, "Sample.description"] }
    if( length(intersect(colnames(tmp2), "Biospecimen.Barcode"))==1){ tmp2 <- tmp2[, "Biospecimen.Barcode"] }
    if( length(intersect(colnames(tmp2), "Comment..TCGA.Barcode."))==1){ tmp2 <- tmp2[, "Comment..TCGA.Barcode."] }
    tmp2 <- as.matrix(tmp2)
    if( ncol(tmp2)!=1){
      if( length( grep("nationwidechildrens", tmplf[i,"FileName"]))==1 || length( grep("biospecimen", tmplf[i,"FileName"]))==1 || length( grep("bio.Level", tmplf[i,"FileName"]))==1   ){
        tmp2 <- read.table(currentTFile, quote="\"",stringsAsFactors = F )[,2]
        tmp2 <- tmp2[grep("TCGA",tmp2)]
        if( length( grep("nationwidechildrens",tmp2 ))==1 || length( grep("biospecimen", tmp2))==1 ){
          if ( length(tmp2)!=0){ tmp2 <- as.matrix(unlist(strsplit(tmp2, ".",fixed = T))) }
        }
      }
    }
    
    tmp2 <- tmp2[grep("TCGA",tmp2)]
    SampleTmp <- unique(substr(tmp2, startK, stopK))
    SampleTmp <- SampleTmp[grep("TCGA",SampleTmp)]
    msgOUT <-  paste(i," ", Tumor," ", tmplf[i,"Type"], " ", tmplf[i,"Species"], " ",   tmplf[i,"Center"], " ", tmplf[i,"Platform"] , " " ,  " .n samples ", length(SampleTmp), sep="")
    print(msgOUT)
    
    if ( length(SampleTmp) < 2) { message("Maybe.....ERROR....Check IT")}
    
    idx<- which(names(TumorDataList) == tmplf[i,"Platform"])
    TumorDataList[[idx]] <- SampleTmp
    
  }
  #Tumor_completeLIST[[k]]<-TumorDataList
  
  setwd(FolderWd)
  return(TumorDataList)
  
}
