TCGAVersionDetailed <- function(Tumor,
                                PlatformType,
                                listSample=0,
                                sdrfFolder,
                                PlatformAndAssociatedData){
  #downloadFolder<-paste(downloadFolder,PlatformType,"/",sep="")
  #.createDirectory(PlatformType)
  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  tmp <- PlatformAndAssociatedData[toupper(PlatformAndAssociatedData$Tumor) == toupper(Tumor)
                                   & toupper(PlatformAndAssociatedData$Platform) == toupper(PlatformType), ]

  key1a       <- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
  Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
  key2a       <- paste0("/",tmp$Folder,"/")
  link        <- paste0(Description, key2a)
  linkInfo    <- .DownloadURL(link)
  version     <- linkInfo[grep("Level_3",linkInfo)]
  version     <- as.matrix(version[-grep("tar.gz",version)])
  versionMat  <- as.data.frame(matrix(0,nrow(version),2))
  colnames(versionMat) <- c("Version","Date")

  aux  <- as.matrix(unlist(strsplit(version, "  ")))
  time <- aux[grep(":",aux)]
  vers <- aux[grep("Level_3",aux)]
  vers <- as.matrix(sapply(strsplit(vers, ">"), function(y) y[2]))
  vers <- as.matrix(sapply(strsplit(vers, "<"), function(y) y[1]))
  versionMat$Version <- vers
  versionMat$Date    <- time

  versionMat <- cbind(versionMat,
                      Samples = matrix(0, nrow(versionMat),1),
                      SizeMbyte = matrix(0, nrow(versionMat),1)
                      )
  print(paste("Found", nrow(versionMat), "Version of", PlatformType, sep = " "))

  for (i in 1: nrow(versionMat)){

    linkToVersion <- paste0(Description,key2a,versionMat$Version[i])
    samplesInfo   <- findSampleVersions(listSample,
                                        linkToVersion,
                                        versionMat$Version[i]
                                        )

    print(paste("Version", i , "of", nrow(versionMat),
                versionMat$Version[i], "...done",
                sep=" ")
          )
    x <- .DownloadURL(linkToVersion)

    if(PlatformType == "illuminahiseq_rnaseq"  ){ x <- x[grep("gene.quantification" , x)]}
    if(PlatformType == "agilentg4502a_07_3"    ){ x <- x[grep("tcga_level3"         , x)]}
    if(PlatformType == "illuminahiseq_rnaseqv2"){ x <- x[grep("rsem.genes.results"  , x)]}
    if(PlatformType == "humanmethylation27"    ){ x <- x[grep("HumanMethylation27"  , x)]}
    if(PlatformType == "humanmethylation450"   ){ x <- x[grep("HumanMethylation450" , x)]}
    if(PlatformType == "illuminaga_mirnaseq"   ){ x <- x[grep("mirna.quantification", x)]}
    if(PlatformType == "genome_wide_snp_6"     ){ x <- x[grep("hg19.seg"            , x)]}

    sizeList <- sapply(strsplit(x, ":"),
                    function(y) {
                      sapply(strsplit(y[2], " "),
                             function(z) z[3])
                    }

             )

    versionMat$SizeMbyte[i] <- getTotalSize(sizeList)
    versionMat$Samples[i]   <- length(sizeList)
    versionMat$Version[i]   <- findSampleVersions(listSample,)

  }
  return(versionMat)
}

# get a detail matrix file, uuid, barcode, size, date
# input: listSample list of bar code to filter output
#        linkToVersion link to version
#        version: father folder?
#        sdrfFolder: folder created by TCGAmanifest
findSampleVersions <- function (listSample,
                                linkToVersion,
                                sdrfFolder="/Users/tiago/Downloads/trash/TCGAsdrf/"){

  lstFileSdrf <- list.files(file.path(sdrfFolder))
  lstFileSdrf_plt <- lstFileSdrf[grep(tolower(PlatformType), tolower(lstFileSdrf))]
  listSample_fromSdrf <- read.delim(paste( sdrfFolder,  lstFileSdrf_plt,sep=""))

  ftpContent <- .DownloadURL(linkToVersion)
  files <-  ftpContent[grep("[a-z0-9]{8}-[a-z0-9]{4}", ftpContent)]

  # creating output matrix
  versionMat  <- as.data.frame(matrix(0,length(files),5))
  colnames(versionMat) <- c("file","uuid","barcode","Date","Size")

  # inserting data into matrix (uuid, data, size, barcode)
  aux  <- as.matrix(unlist(strsplit(files, "  ")))
  filesName <- aux[grep("[a-z0-9]{8}-[a-z0-9]{4}",aux)]
  fileAux <- unlist(strsplit(unlist(strsplit(filesName,">")),"="))
  versionMat$file <- fileAux[-grep("<",fileAux)]
  versionMat$uuid <- substr(filesName, 17,52)
  versionMat$Date <- aux[grep(":",aux)]
  versionMat$Size <- aux[grep("[0-9]K|[0-9]M",aux)]

  if(length(listSample)!=0){
    if(PlatformType == "illuminahiseq_rnaseqv2" || PlatformType == "illuminahiseq_totalrnaseqv2"){
      print(paste("Finding uuid for", length(listSample), "samples with TCGA barcode selected",sep=" "))

      for( g in 1:length(files) ){
        tmpsdrf <-  listSample_fromSdrf[ as.character(listSample_fromSdrf$Extract.Name) ==  as.character(versionMat$uuid[g]),]
        tmpbarcode <- as.character(tmpsdrf$Comment..TCGA.Barcode.)[1]
        if(length(tmpbarcode) == 0 ) next
        versionMat$barcode[g] <- as.character(tmpsdrf$Comment..TCGA.Barcode.)[1]
      }
      versionMat <- versionMat[!(is.na(versionMat$barcode)),]
      versionMat <- versionMat[substr(versionMat$barcode,1,nchar(listSample[1])) %in% listSample,]
    }
  }
  return (versionMat)
}
# param sizeList - list of sizes in KB, MB
# output total size in MB
getTotalSize <- function(sizeList){
  sizeK  <- sizeList[grep("K",sizeList)]
  sizeM  <- sizeList[grep("M",sizeList)]
  totalK <- round(sum(as.numeric(gsub("K","",sizeK)))/1000)
  totalM <- sum(as.numeric(gsub("M","",sizeM)))
  return (totalM + totalK)
}

#' @title TCGA Version
#'
#' @description  TCGA Version
#'
#' @param Tumor  a character string indicating the cancer type for
#'        which to download data. Options include ACC, BLCA, BRCA,
#'        CESC, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML,
#'        LGG, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, SARC, SKCM, STAD,
#'        THCA, UCEC, UCS. Look at https://tcga-data.nci.nih.gov/tcga/
#'        for Available Cancer Types.
#' @param PlatformType "illuminahiseq_rnaseq"
#' @param PlatformAndAssociatedData data frame 615 observations of 12 variables,
#'        indicating the different characteristics of the data
#'        e.g. tumour, type, species.
#' @import RCurl httr bitops
#' @export

TCGAVersion <- function(Tumor, PlatformType,PlatformAndAssociatedData){
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
