TCGADownload <-
function(Tumor, PlatformAndAssociatedData,PlatformType,nsample=0, listSample=0, newsample =FALSE){
  require(RCurl)
  FolderWd <- getwd()
  #setwd(downloadFolder)
  
  
  

  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  tmp <- PlatformAndAssociatedData[toupper(PlatformAndAssociatedData$Tumor) == toupper(Tumor)
                                   & toupper(PlatformAndAssociatedData$Platform) == toupper(PlatformType), ]
  
  key1a <- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
  Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
  key2a <- paste("/",tmp$Folder,"/",sep="")
  
  if(PlatformType == "illuminaga_dnaseq"){
    toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_2", KeyGrep2 = "somatic.maf")
    toDdl <- paste(Description, key2a, toDdl, sep = "")
    
    x <- .DownloadURL(toDdl)
    x <- strsplit(x, "\t")
    x <- x[-1]
    x <- matrix(unlist(x), nrow = length(x), byrow = T)
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    
    return(x)
  }
  
  
  lstFileSdrf <- list.files(file.path( FolderWd))
  lstFileSdrf_plt <- lstFileSdrf[grep(tolower(PlatformType), tolower(lstFileSdrf))]
  listSample_fromSdrf <- read.delim(paste( FolderWd,"/",  lstFileSdrf_plt,sep=""))
  
  
  toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_3", KeyGrep2 = "MANIFEST.txt") 
  toDdl <- paste(Description, key2a, toDdl, sep = "")
  
  x <- .DownloadURL(toDdl)
  x <- sapply(strsplit(x, "  "), function(y) y[2])
  
  if(PlatformType == "illuminahiseq_rnaseq"){  x <- x[grep("gene.quantification", x)] }
  if(PlatformType == "agilentg4502a_07_3"){    x <- x[grep("tcga_level3", x)]}
  if(PlatformType == "illuminahiseq_rnaseqv2"){ x <- x[grep("rsem.genes.results", x)] }
  if(PlatformType == "humanmethylation27"){ x <- x[grep("HumanMethylation27", x)] }
  if(PlatformType == "humanmethylation450"){ x <- x[grep("HumanMethylation450", x)] }
  if(PlatformType == "illuminaga_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }
  if(PlatformType == "illuminahiseq_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }
  
  
  if(PlatformType == "genome_wide_snp_6"){ x <- x[grep("hg19.seg", x)]
                                           x <- x[-grep("nocnv", x)]}
  
  if(PlatformType == "mda_rppa_core"){  x <- x[grep("protein_expression", x)] }
  
  
  
  if(length(listSample)!=0){
  if(PlatformType == "illuminahiseq_rnaseqv2"){
    xSplit <- gsub(".rsem.genes.results", "", sapply(strsplit(x, "unc.edu"), function(x) x[2]))
    xSplit2 <- substr(xSplit, 2,37)
    xuuid <- as.data.frame(cbind( list = as.character(x), uuid = as.character(xSplit2), barcode = as.character(xSplit2)))
    xuuid$list <- as.character(xuuid$list)
    xuuid$uuid <- as.character(xuuid$uuid)
    xuuid$barcode <- as.character(xuuid$barcode)
    print(paste("Finding uuid for ", length(listSample), " samples with TCGA barcode selected",sep=""))
    for( g in 1 : nrow(xuuid) ){
      
      tmpsdrf <-  listSample_fromSdrf[ as.character(listSample_fromSdrf$Extract.Name) ==  as.character(xuuid$uuid[g]),]
      tmpbarcode <- as.character(tmpsdrf$Comment..TCGA.Barcode.)[1]
      if(length(tmpbarcode) == 0 ) next
    xuuid$barcode[g] <- as.character(tmpsdrf$Comment..TCGA.Barcode.)[1]
  }
  xuuid_complete <- xuuid[!(is.na(xuuid$barcode)),]
  xuuid_complete_todown <- xuuid_complete[substr(xuuid_complete$barcode,1,nchar(listSample[1])) %in% listSample,]
   x <- x [x %in% xuuid_complete_todown$list]
  }
  }
   
  plt <- PlatformAndAssociatedData
  plt2 <- plt[toupper(plt$Tumor) == toupper(Tumor) & plt$Platform == PlatformType,]
  
  version <- gsub(Description,"",toDdl)
  version <- gsub(key2a,"",version)
  version <- gsub("MANIFEST.txt","",version)
  
  print(paste("Found ", length(x), " samples of ",Tumor, "version: ", version, "..." ))
   
  setwd("..")
  .createDirectory(PlatformType)
  
  downloadFolder<-paste(getwd(),PlatformType,sep="/")
  
  
 print(downloadFolder)
 setwd(downloadFolder)
  
 
 if (newsample == TRUE){
   lstDownAlready <-  list.files(downloadFolder)
   print(paste("Found ", length(lstDownAlready), " samples downloaded"),sep="")
   sampleNew <- setdiff(x,lstDownAlready) 
   print(paste("Found ", length(sampleNew), " new samples to download"),sep="")
   x <- sampleNew
 }
 
 samplesList <- paste(gsub("MANIFEST.txt", "", toDdl), x, sep = "")
 if (nsample!=0){
   samplesList <- samplesList[1:nsample]
 }
 
 
 
 
for(i in 1:length(samplesList)){
  tmp2 <- .DownloadData_fromURL(url =samplesList[i] ,sep="", header = T)
  #filename <- paste(downloadFolder, x[i], sep = "")
  filename <- x[i]
  write.csv(tmp2, filename)
  print(paste(x[i], " ... sample n. ", i, " of ", length(samplesList), sep = ""))
}

  
  lf <- list.files(downloadFolder)
  lf1<-lf
  lf <- paste(downloadFolder, lf, sep = "/")
  
  tmpData <- read.csv(lf[1], sep = ",", stringsAsFactors = FALSE)
  geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
  tmpData <- tmpData[,-1]
  
  if(PlatformType == "agilentg4502a_07_3"){
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
  
  if(PlatformType == "illuminahiseq_rnaseqv2"){
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
  
  
  setwd(FolderWd)
  return(geData)
  
}
