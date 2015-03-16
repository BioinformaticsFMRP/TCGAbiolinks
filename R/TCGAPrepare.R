TCGAPrepare<- function(Tumor, PlatformAndAssociatedData, sdrfFolder = "", downloadFolder = "",PlatformType,nsample=0, listSample=0, newsample =FALSE){
  require(RCurl)
  
  
  
  lf <- list.files(downloadFolder)
  lf1<-lf
  lf <- paste(downloadFolder, lf, sep = "/")
  
  tmpData <- read.csv(lf[1], sep = ",", stringsAsFactors = FALSE)
  geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
  tmpData <- tmpData[,-1]
  
  if(PlatformType == "agilentg4502a_07_3" | PlatformType == "agilentg4502a_07_2" | PlatformType == "agilentg4502a_07_1"){
    tmpData <- tmpData[-1, ]
    geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
    rownames(geData) <- tmpData$Hybridization.REF
    colNames <- rep("", length(lf))
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      tmpData <- tmpData[,-1]
      colNames[i] <- colnames(tmpData)[2]
      tmpData <- tmpData[-1, 2]
      tmpData[which(tmpData == "null")] <- NA
      geData[, i] <- as.numeric(tmpData)
      print(i)
    }
    colnames(geData) <- substr(gsub("\\.", "-", colNames),1,16)
  }
  
  if(PlatformType == "illuminahiseq_rnaseq"){
    rownames(geData) <- tmpData$gene
    colnames(geData) <- substr(paste("TCGA", sapply(strsplit(lf, "TCGA"), function(y) y[2]), sep = ""), 1, 28)
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      tmpData <- tmpData[,-1]
      expr <- tmpData$raw_counts
      geData[, i] <- expr
      print(i)
    }
  }
  
  if(PlatformType == "illuminahiseq_rnaseqv2" || PlatformType == "illuminahiseq_totalrnaseqv2"){
    rownames(geData) <- tmpData$gene_id
    path2 <- paste(Description, key2a,plt2$FileName,sep="")
    xpath <- .DownloadURL(path2)
    nc<-length(unlist(strsplit(xpath[1], "\t")))
    xmat <- matrix( unlist(strsplit(xpath, "\t") ), ncol=nc,byrow=T)
    colnames(xmat)<-gsub(" ","_",xmat[1,])
    
    xmat<-xmat[-1,]
    xmat <- as.data.frame(xmat,stringsAsFactors=F)
    xmat2 <- xmat[grep("rsem.genes.results",xmat$Derived_Data_File),]
    xmat2 <- xmat[sapply(lf1,function(y) grep(y, xmat$Derived_Data_File)),]
    
    colnames(geData)<-xmat2$'Comment_[TCGA_Barcode]'
    
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      expr <- tmpData$raw_count
      geData[, i] <- expr
      print(i)
    }
  }
  
  
  
  if(PlatformType == "humanmethylation27" || PlatformType == "humanmethylation450" ){
    tmpData <- tmpData[-1, ]
    geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
    rownames(geData) <- tmpData$Hybridization.REF
    colNames <- rep("", length(lf))
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      colNames[i] <- colnames(tmpData)[2]
      tmpData <- tmpData[-1, 2]
      tmpData[which(tmpData == "null")] <- NA
      geData[, i] <- as.numeric(tmpData)
      print(i)
    }
    colnames(geData) <- substr(gsub("\\.", "-", colNames),1,28)
  }
  
  
  
  if(PlatformType == "illuminaga_mirnaseq" || PlatformType == "illuminahiseq_mirnaseq"){
    rownames(geData) <- tmpData$miRNA_ID
    colnames(geData) <- substr(paste("TCGA", sapply(strsplit(lf, "TCGA"), function(y) y[2]), sep = ""), 1, 28)
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      expr <- tmpData$read_count
      geData[, i] <- expr
      print(i)
    }
  }
  
  if(PlatformType == "mda_rppa_core"){
    tmpData <- tmpData[-1, ]
    geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
    rownames(geData) <- tmpData$Sample.REF
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      tmpData <- tmpData[-1, 2]
      tmpData[which(tmpData == "null")] <- NA
      geData[, i] <- as.numeric(tmpData)
      print(i)
    }
    colNames <- gsub(".txt", "", sapply(strsplit(lf1, "Level_3."), function(x) x[2]))
    
    toDdl <- .DownloaDmageTAB_sdrf(Description, key2a, KeyGrep1 = "mage-tab", KeyGrep2 = "array_design.txt")
    toDdl <- paste(Description, key2a, toDdl, sep = "")
    x <- .DownloadURL(toDdl)
    x <- strsplit(x, "\t")
    x <- matrix(unlist(x), nrow = length(x), byrow = T)
    x <- gsub("\r", "", x)
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    rownames(x) <- x[, "UUID"]
    
    colnames(geData) <- x[colNames, "Sample.barcode"]
  }
  
  if(PlatformType == "illuminahiseq_dnaseqc" | PlatformType == "hg-cgh-415k_g4124a" | PlatformType == "hg-cgh-244a"){
    geData <- vector("list", length(x))
    for(i in 1:length(lf)) geData[[i]] <- read.csv(lf[i], sep = ",", stringsAsFactors = FALSE)
    names(geData) <- substr(paste("TCGA", sapply(strsplit(lf, "TCGA"), function(y) y[2]), sep = ""), 1, 28)
  }
  
  if(PlatformType == "genome_wide_snp_6"){
    
    path2 <- paste(Description, key2a,plt2$FileName,sep="")
    xpath <- .DownloadURL(path2)
    nc<-length(unlist(strsplit(xpath[1], "\t")))
    xmat <- matrix( unlist(strsplit(xpath, "\t") ), ncol=nc,byrow=T)
    colnames(xmat)<-gsub(" ","_",xmat[1,])
    
    xmat<-xmat[-1,]
    xmat <- as.data.frame(xmat,stringsAsFactors=F)
    lf2<-gsub("hg19","hg18",lf1)
    lf3 <- intersect(lf2, xmat$Derived_Array_Data_File)
    
    xmat2 <- xmat[sapply(lf3,function(y) grep(y, xmat$Derived_Array_Data_File)),]
    #grep(lf2[5], xmat$Derived_Array_Data_File)
    
    geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf3))
    colnames(geData)<-xmat2$'Comment_[TCGA_Barcode]'
    
    # rownames(geData) <- tmpData$gene
    # colnames(geData) <- substr(paste("TCGA", sapply(strsplit(lf, "TCGA"), function(y) y[2]), sep = ""), 1, 16)
    for(i in 1:length(lf)){
      tmpData <- read.delim(lf[i], stringsAsFactors = FALSE, sep = ",")
      expr <- tmpData$Segment_Mean
      geData[, i] <- expr
      print(i)
    }
  }
  
  if(PlatformType == "huex-1_0-st-v2"){
    tmpData <- tmpData[-1, ]
    rownames(tmpData) <- tmpData$Hybridization.REF
    tmpData <- tmpData[, -1]
    geData <- apply(tmpData, 2, as.numeric)
    rownames(geData) <- rownames(tmpData)
    colNames <- sub("X", "", colnames(geData))
    
    
    toDdl <- .DownloaDmageTAB_sdrf(Description, key2a, KeyGrep1 = "mage-tab", KeyGrep2 = "sdrf.txt")
    toDdl <- paste(Description, key2a, toDdl, sep = "")
    
    x <- .DownloadURL(toDdl)
    x <- strsplit(x, "\t")
    x <- matrix(unlist(x), nrow = length(x), byrow = T)
    x <- gsub("\r", "", x)
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    rownames(x) <- x[, "Hybridization Name"]
    
    colnames(geData) <- sub("_LE", "", x[colNames, "Labeled Extract Name"])
  }
  
  if(PlatformType == "human1mduo" | PlatformType == "humanhap550"){
    geData <- tmpData
  }
  
  if(PlatformType == "ht_hg-u133a"){
    tmpData <- tmpData[-1, ]
    
    geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
    rownames(geData) <- tmpData$Hybridization.REF
    colNames <- rep("", ncol(geData))
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      tmpData <- tmpData[-1, -1]
      expr <- tmpData[, 2]
      expr[expr == "N/A"] <- NA
      geData[, i] <- as.numeric(expr)
      colNames[i] <- gsub(".", "-", colnames(tmpData)[2], fixed = T)
      print(i)
    }
    
    toDdl <- .DownloaDmageTAB_sdrf(Description, key2a, KeyGrep1 = "mage-tab", KeyGrep2 = "sdrf.txt")
    toDdl <- paste(Description, key2a, toDdl, sep = "")
    
    x <- .DownloadURL(toDdl)
    x <- strsplit(x, "\t")
    x <- matrix(unlist(x), nrow = length(x), byrow = T)
    x <- gsub("\r", "", x)
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    rownames(x) <- x[, "Hybridization Name"]
    
    colnames(geData) <- x[colNames, "Comment [TCGA Barcode]"]
  }
  
  if(PlatformType == "human1mduo" | PlatformType == "humanhap550"){
    geData <- tmpData
  }
  
  
  if(PlatformType == "illuminadnamethylation_oma003_cpi" | PlatformType == "illuminadnamethylation_oma002_cpi" | PlatformType == "hg-u133_plus_2" | PlatformType == "h-mirna_8x15kv2" | PlatformType == "h-mirna_8x15k"){
    tmpData <- tmpData[-1, ]
    geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
    rownames(geData) <- tmpData$Hybridization.REF
    colNames <- rep("", ncol(geData))
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      tmpData <- tmpData[-1, -1]
      expr <- tmpData[, 2]
      expr[expr == "N/A"] <- NA
      geData[, i] <- as.numeric(expr)
      colNames[i] <- gsub(".", "-", colnames(tmpData)[2], fixed = T)
      print(i)
    }
    colnames(geData) <- colNames
  }
  
  
  setwd(FolderWd)
  return(geData)
  

  
}
