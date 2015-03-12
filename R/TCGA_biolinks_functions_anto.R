TCGADownload <- function(Tumor, PlatformAndAssociatedData, sdrfFolder = "", downloadFolder = "",PlatformType,nsample=0, listSample=0, newsample =FALSE){
  require(RCurl)
  FolderWd <- getwd()
  setwd(downloadFolder)
  .createDirectory(PlatformType)
  downloadFolder<-paste(downloadFolder,PlatformType,sep="/")

  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  tmp <- PlatformAndAssociatedData[toupper(PlatformAndAssociatedData$Tumor) == toupper(Tumor)
                                   & toupper(PlatformAndAssociatedData$Platform) == toupper(PlatformType), ]

  key1a <- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
  Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
  key2a <- paste("/",tmp$Folder,"/",sep="")
  
  if(PlatformType == "illuminaga_dnaseq" | PlatformType == "solid_dnaseq" | PlatformType == "solid_dnaseq_curated" | PlatformType == "mixed_dnaseq_curated" | PlatformType == "mixed_dnaseq" | 
       PlatformType == "illuminahiseq_dnaseq_automated" | PlatformType == "illuminaga_dnaseq_curated" | PlatformType == "illuminaga_dnaseq_automated" | PlatformType == "illuminaga_dnaseq"){

    toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_2", KeyGrep2 = "somatic.maf")
    toDdl <- paste(Description, key2a, toDdl, sep = "")

    x <- .DownloadURL(toDdl)
    x <- strsplit(x, "\t")
    if(PlatformType == "illuminaga_dnaseq_curated" | PlatformType == "illuminaga_dnaseq") x <- x[-c(1:4)] else x <- x[-1]
    x <- matrix(unlist(x), nrow = length(x), byrow = T)
    colnames(x) <- x[1, ]
    x <- x[-1, ]

    return(x)
  }


  lstFileSdrf <- list.files(file.path(sdrfFolder))
  lstFileSdrf_plt <- lstFileSdrf[grep(tolower(PlatformType), tolower(lstFileSdrf))]
  listSample_fromSdrf <- read.delim(paste( sdrfFolder,  lstFileSdrf_plt,sep=""))

  
  #if(PlatformType == "illuminadnamethylation_oma003_cpi" | PlatformType == "illuminadnamethylation_oma002_cpi"){
  #  toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_2", KeyGrep2 = "MANIFEST.txt") 
  #  toDdl <- paste(Description, key2a, toDdl, sep = "")
  #}else{
    toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_3", KeyGrep2 = "MANIFEST.txt") 
    toDdl <- paste(Description, key2a, toDdl, sep = "")
  #}
  

  x <- .DownloadURL(toDdl)
  x <- sapply(strsplit(x, "  "), function(y) y[2])

  if(PlatformType == "illuminahiseq_rnaseq"){  x <- x[grep("gene.quantification", x)] }
  if(PlatformType == "agilentg4502a_07_3"){    x <- x[grep("tcga_level3", x)]}
  if(PlatformType == "illuminahiseq_rnaseqv2" || PlatformType == "illuminahiseq_totalrnaseqv2"){ x <- x[grep("rsem.genes.results", x)] }
  if(PlatformType == "humanmethylation27"){ x <- x[grep("HumanMethylation27", x)] }
  if(PlatformType == "humanmethylation450"){ x <- x[grep("HumanMethylation450", x)] }
  if(PlatformType == "illuminadnamethylation_oma003_cpi"){ x <- x[grep("IlluminaDNAMethylation_OMA003_CPI", x)] }
  if(PlatformType == "illuminaga_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }
  if(PlatformType == "illuminahiseq_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }

  if(PlatformType == "genome_wide_snp_6"){ x <- x[grep("hg19.seg", x)]
                                           x <- x[-grep("nocnv", x)]}
  if(PlatformType == "illuminahiseq_dnaseqc"){ x <- x[grep("Segment", x)] }

  #if(PlatformType == "humanhap550"){ x <- x[grep("seg.txt", x)] }
  #if(PlatformType == "human1mduo"){ x <- x[grep("seg.txt", x)] }
  

  if(PlatformType == "mda_rppa_core"){  x <- x[grep("protein_expression", x)] }



  if(length(listSample)!=0){
  if(PlatformType == "illuminahiseq_rnaseqv2" || PlatformType == "illuminahiseq_totalrnaseqv2"){
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


  # .DownloadData_fromURL("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/barcode/TCGA-BH-A0H9-01A")



  plt <- PlatformAndAssociatedData
  plt2 <- plt[toupper(plt$Tumor) == toupper(Tumor) & plt$Platform == PlatformType,]

  version <- gsub(Description,"",toDdl)
  version <- gsub(key2a,"",version)
  version <- gsub("MANIFEST.txt","",version)

  print(paste("Found ", length(x), " samples of ",Tumor, "version: ", version, "..." ))

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


# start create matrix

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

  if(PlatformType == "illuminahiseq_dnaseqc"){
    geData <- vector("list", length(x))
    for(i in 1:length(lf)) geData[[i]] <- read.csv(lf[i], sep = ",", stringsAsFactors = FALSE)
    names(geData) <- substr(x, 1, 28)
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



TCGAVersion <- function(Tumor, PlatformType,PlatformAndAssociatedData){

  require(RCurl)
  #downloadFolder<-paste(downloadFolder,PlatformType,"/",sep="")
  #.createDirectory(PlatformType)
  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  tmp <- PlatformAndAssociatedData[toupper(PlatformAndAssociatedData$Tumor) == toupper(Tumor)
                                   & toupper(PlatformAndAssociatedData$Platform) == toupper(PlatformType), ]

  key1a <- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
  Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
  key2a <- paste("/",tmp$Folder,"/",sep="")

  #toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_3", KeyGrep2 = "MANIFEST.txt")
  #toDdl <- paste(Description, key2a, toDdl, sep = "")
  #x <- .DownloadURL(toDdl)
  #x <- sapply(strsplit(x, "  "), function(y) y[2])

  toDdl <- paste(Description, key2a, sep = "")
  x <- .DownloadURL(toDdl)
  xver <- x[grep("Level_3",x)]
  xver <- as.matrix(xver[-grep("tar.gz",xver)])
  xverMat <- as.data.frame(matrix(0,nrow(xver),2))
  colnames(xverMat)<-c("Version","Date")

  for( i in 1: nrow(xverMat)){
    xtmp1 <- xver[i]
    xver2 <- as.matrix(unlist(strsplit(xtmp1, "  ")))
    timeVer <- xver2[grep(":",xver2)]
    Vers <- xver2[grep("Level_3",xver2)]
    Vers  <- as.matrix(sapply(strsplit(Vers, ">"), function(y) y[2]))
    Vers  <- as.matrix(sapply(strsplit(Vers, "<"), function(y) y[1]))
    xverMat$Version[i]<- Vers
    xverMat$Date[i]<-timeVer
  }


  xverMat <- cbind(xverMat, Samples = matrix(0, nrow(xverMat),1), SizeMbyte = matrix(0, nrow(xverMat),1))
  print(paste("Found ", nrow(xverMat), " Version of ", PlatformType,sep=""))

  for( i in 1: nrow(xverMat)){

    todown1<- paste(Description,key2a,xverMat$Version[i],sep="")
    print(paste("Version ", i , " of ", nrow(xverMat), " ", xverMat$Version[i], " ...done",sep=""))
    x <- .DownloadURL(todown1)

    if(PlatformType == "illuminahiseq_rnaseq"){  x <- x[grep("gene.quantification", x)] }
    if(PlatformType == "agilentg4502a_07_3"){    x <- x[grep("tcga_level3", x)]}
    if(PlatformType == "illuminahiseq_rnaseqv2"){ x <- x[grep("rsem.genes.results", x)] }
    if(PlatformType == "humanmethylation27"){ x <- x[grep("HumanMethylation27", x)] }
    if(PlatformType == "humanmethylation450"){ x <- x[grep("HumanMethylation450", x)] }
    if(PlatformType == "illuminaga_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }
    if(PlatformType == "genome_wide_snp_6"){ x <- x[grep("hg19.seg", x)]}

    x2<- sapply(strsplit(x, ":"), function(y) y[2])
    x3<- sapply(strsplit(x2, " "), function(y) y[3])
    sizeK <- x3[grep("K",x3)]
    sizeM <- x3[grep("M",x3)]
    sizeK_1 <- as.numeric(gsub("K","",sizeK))
    sizeM_1 <- as.numeric(gsub("M","",sizeM))
    sizeK_2<- round(sum(sizeK_1)/1000)
    sizeM_2<- sum(sizeM_1)
    sizeTot<- sizeK_2+sizeM_2
    xverMat$SizeMbyte[i]<-sizeTot
    xverMat$Samples[i]<-length(x3)

  }
  return(xverMat)
}

TCGAmanifest <- function(Tumor){
  FolderWd <- getwd()
  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  .createDirectory("TCGAsdrf")
  setwd("./TCGAsdrf")
  tmplf <- PlatformAndAssociatedData[PlatformAndAssociatedData$Tumor %in% tolower(Tumor),]
  token <- 1
  for ( w in token: nrow( tmplf)){
    tmp <- tmplf[w,]
    if(tmp["FileName"] == "") next
    key1a<- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
    Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
    DescriptionSite <- paste(Description, tmp$Folder, tmplf[w,"FileName"], sep="/")
    x <- .DownloadURL(DescriptionSite)
    #x <- .DownloadData_fromURL(DescriptionSite)
    writeLines(x, gsub("/","_",tmplf[w,"FileName"]))

    print(paste(w, " of ", nrow(tmplf ), " ", tmplf[w,"FileName"]),sep="")
  }
  setwd(FolderWd)
}

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




.TCGAQuery2 <- function(Tumor,siteTCGA){

  TumorData <- matrix(0, 1, 8)
  colnames(TumorData) <- c("Total", "Exome", "SNP", "Methylation", "mRNA", "miRNA", "Clinical","Protein")
  rownames(TumorData) <- Tumor
  LevelsPlatforms <- c("Level_1", "Level_2", "Level_3")

  #for ( i in 2: ncol(TumorData)){


  i<-5
  Type <-  colnames(TumorData)[i]
  tmp <- PlatformAndAssociatedData[PlatformAndAssociatedData$Type%in% Type,]

  tmp2a <- tmp[ grep(tolower(Tumor),tmp$Tumor),]

  if( nrow(tmp2a)!=0){

    Species <- unique(tmp2a$Species)



    for( j in 1:length(Species)){
      SpecieCurr <- Species[j]
      tmp3 <- tmp2a[tmp2a$Species %in% SpecieCurr, ]
      Centers <- unique(tmp3$Center)

      for( k in 1:length(Centers)){
        CenterCurr <- Centers[k]
        tmp3b <- tmp3[tmp3$Center %in% CenterCurr,]

        Platforms <- unique(tmp3b$Platform)

        for( q in 1:length(Platforms)){

          PlatformCurr <- Platforms[q]
          tmp4 <- tmp3b[tmp3b$Platform %in% PlatformCurr, ]

          key1<- paste(unique(tmp4$CenterType), unique(tmp4$Center), unique(tmp4$Platform), sep="/")
          Description <- paste(siteTCGA, tolower(Tumor), "/",key1, sep="")

          if( length(grep("agilentg4502a_07",unique(tmp4$Platform))) > 0   || unique(tmp4$Platform) ==  "ht_hg-u133a" || unique(tmp4$Platform) ==  "hg-u133_plus_2" || unique(tmp4$Platform) ==  "illuminaga_mrna_dge" ){
            Description <- paste(Description, "/transcriptome/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)



            for( w in 1:length(LevelsPlatforms)){

              siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[w],Description)
              tmp2 <- DownloadManifest(siteNewLevel)


              if(length(grep("agilentg4502a_07",unique(tmp4$Platform))) > 0){
                NumberSample <- length(unique(substr(tmp2, 1, 23)))
              }

              if(unique(tmp4$Platform) ==  "ht_hg-u133a"){
                NumberSample <- length(unique(substr(tmp2, 44, 58)))
              }

              if(unique(tmp4$Platform) ==  "hg-u133_plus_2" || unique(tmp4$Platform) ==  "illuminaga_mrna_dge"){
                NumberSample <- length(unique(substr(tmp2, 1, 16)))
              }

              message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[w] ,  " .n samples ", NumberSample, sep="")
              print(message)

            } #end for LevelsPlatform
          } #end platform Exp-Gene


          if(unique(tmp4$Platform) ==  "huex-1_0-st-v2"){
            Description <- paste(Description, "/exon/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key="mage-tab",Description)
            siteNewLevelSdrf <- DownloadSdrf(siteNewLevel)
            tmp2 <-  siteNewLevelSdrf$Comment..TCGA.Barcode.
            tmp2 <- tmp2[grep("TCGA",tmp2)]
            NumberSample <- length(unique(substr(tmp2, 1, 16)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " , " .n samples ", NumberSample, sep="")
            print(message)
          }



          if(unique(tmp4$Platform) ==  "illuminahiseq_rnaseq" || unique(tmp4$Platform) ==  "illuminaga_rnaseq"){
            Description <- paste(Description, "/rnaseq/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[3],Description)
            tmp2 <- DownloadManifest(siteNewLevel)
            NumberSample <- length(unique(substr(tmp2, 13, 29)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[3] ,  " .n samples ", NumberSample, sep="")
            print(message)
          }



          if(unique(tmp4$Platform) ==  "illuminahiseq_rnaseqv2" || unique(tmp4$Platform) ==  "illuminaga_rnaseqv2" ){
            Description <- paste(Description, "/rnaseqv2/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[3],Description)
            tmp2 <- DownloadManifest(siteNewLevel)

            NumberSample <- length(unique(substr(tmp2, 13, 29)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[3] ,  " .n samples ", NumberSample, sep="")
            print(message)
          }


          if(unique(tmp4$Platform) ==  "illuminahiseq_totalrnaseqv2"){
            Description <- paste(Description, "/totalrnaseqv2/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[3],Description)
            tmp2 <- DownloadManifest(siteNewLevel)

            NumberSample <- length(unique(substr(tmp2, 9, 44)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[3] ,  " .n samples ", NumberSample, sep="")
            print(message)
          }



        } # end platform
      } #end for Centers
    } # end for species



    #}
  }

}

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
      print(paste("Reconnection to the url #",count1,sep=""))
      Sys.sleep(1)
    }
    if(count1 == 20) stop(paste("The url (",Site,") is not existing or not available.
                          Please check your internet connection or contact the mantainers for an update.",sep=""))
  }

  if(interactive() && ("ssl" %in% names(curlVersion()$features))) {
    bo2 = T
    count <- 0
    while(bo2){
      x = try(getURLContent(Site, verbose = F, curl = handle),silent = T)
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
