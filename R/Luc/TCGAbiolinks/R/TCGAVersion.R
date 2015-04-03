TCGAVersion <-
function(Tumor, PlatformType,PlatformAndAssociatedData){
  
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
