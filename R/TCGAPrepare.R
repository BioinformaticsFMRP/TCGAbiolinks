#' @title TCGA Prepare
#' @description Prepare data matrices for downstream analysis, ready to use also with other R packages using data previously downloaded using TCGAdownload
#' @param sdrfFolder TCGAQuery location of the output data saving
#' @param downloadFolder TCGADownload location of the output data saving
#' @param PlatformType Platform as output of TCGAQuery
#' @seealso TCGAQuery
#' @examples
#' \dontrun{
#'    TCGAPrepare(sdrfFolder = 'folder1', downloadFolder= 'folder2')
#' }
#' @export
#' @importFrom downloader download
#' @return Data matrix
TCGAPrepare <- function(sdrfFolder = "", downloadFolder = "", PlatformType =""){

  Description <- "To fix with new version of TCGAQuery"
  key2 <- "To fix with new version of TCGAQuery"
  key2a <- "To fix with new version of TCGAQuery"
  plt2 <-  "To fix with new version of TCGAQuery"
  DownloadHTML <- function(out ="To fix with new version of TCGAQuery") {return(out)}
  .DownloaDmageTAB_sdrf <- function(out ="To fix with new version of TCGAQuery") {return(out)}

  FolderWd <- getwd()
  lf <- list.files(downloadFolder)
  lf1<-lf
  lf <- paste(downloadFolder, lf, sep = "/")

  tmpData <- read.csv(lf[1], sep = ",", stringsAsFactors = FALSE)
  geData <- matrix(0, nrow = nrow(tmpData), ncol = length(lf))
  tmpData <- tmpData[,-1]

  if(PlatformType == "agilentg4502a_07_3" |
       PlatformType == "agilentg4502a_07_2" |
       PlatformType == "agilentg4502a_07_1"){
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
    colnames(geData) <- substr(paste("TCGA",
                                     sapply(strsplit(lf, "TCGA"),
                                            function(y) y[2]), sep = ""), 1, 28)
    for(i in 1:length(lf)){
      tmpData <- read.csv(lf[i], stringsAsFactors = FALSE, sep = ",")
      tmpData <- tmpData[,-1]
      expr <- tmpData$raw_counts
      geData[, i] <- expr
      print(i)
    }
  }

  if(PlatformType == "illuminahiseq_rnaseqv2" ||
       PlatformType == "illuminahiseq_totalrnaseqv2"){
    rownames(geData) <- tmpData$gene_id
    path2 <- paste(Description, key2a,plt2$FileName,sep="")
    xpath <- DownloadHTML(path2)
    nc<-length(unlist(strsplit(xpath[1], "\t")))
    xmat <- matrix( unlist(strsplit(xpath, "\t") ), ncol=nc,byrow=TRUE)
    colnames(xmat)<-gsub(" ","_",xmat[1,])

    xmat<-xmat[-1,]
    xmat <- as.data.frame(xmat,stringsAsFactors=FALSE)
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



  if(PlatformType == "humanmethylation27" ||
       PlatformType == "humanmethylation450" ){
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

  if(PlatformType == "illuminaga_mirnaseq" ||
       PlatformType == "illuminahiseq_mirnaseq"){
    rownames(geData) <- tmpData$miRNA_ID
    colnames(geData) <- substr(paste("TCGA",
                                     sapply(strsplit(lf, "TCGA"),
                                            function(y) y[2]), sep = ""), 1, 28)
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
    colNames <- gsub(".txt", "", sapply(strsplit(lf1, "Level_3."),
                                        function(x) x[2]))

    toDdl <- .DownloaDmageTAB_sdrf(Description, key2a,
                                   KeyGrep1 = "mage-tab",
                                   KeyGrep2 = "array_design.txt")
    toDdl <- paste(Description, key2a, toDdl, sep = "")
    x <- DownloadHTML(toDdl)
    x <- strsplit(x, "\t")
    x <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
    x <- gsub("\r", "", x)
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    rownames(x) <- x[, "UUID"]

    colnames(geData) <- x[colNames, "Sample.barcode"]
  }

  if(PlatformType == "illuminahiseq_dnaseqc" |
       PlatformType == "hg-cgh-415k_g4124a" |
       PlatformType == "hg-cgh-244a"){
    geData <- vector("list", length(x))
    for(i in 1:length(lf)) geData[[i]] <- read.csv(lf[i], sep = ",",
                                                   stringsAsFactors = FALSE)
    names(geData) <- substr(paste("TCGA",
                                  sapply(strsplit(lf, "TCGA"),
                                         function(y) y[2]), sep = ""), 1, 28)
  }

  if(PlatformType == "genome_wide_snp_6"){

    path2 <- paste(Description, key2a,plt2$FileName,sep="")
    xpath <- DownloadHTML(path2)
    nc<-length(unlist(strsplit(xpath[1], "\t")))
    xmat <- matrix( unlist(strsplit(xpath, "\t") ), ncol=nc,byrow=TRUE)
    colnames(xmat)<-gsub(" ","_",xmat[1,])

    xmat<-xmat[-1,]
    xmat <- as.data.frame(xmat,stringsAsFactors=FALSE)
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

    x <- DownloadHTML(toDdl)
    x <- strsplit(x, "\t")
    x <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
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
      colNames[i] <- gsub(".", "-", colnames(tmpData)[2], fixed = TRUE)
      print(i)
    }

    toDdl <- .DownloaDmageTAB_sdrf(Description, key2a, KeyGrep1 = "mage-tab", KeyGrep2 = "sdrf.txt")
    toDdl <- paste(Description, key2a, toDdl, sep = "")

    x <- DownloadHTML(toDdl)
    x <- strsplit(x, "\t")
    x <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
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
      colNames[i] <- gsub(".", "-", colnames(tmpData)[2], fixed = TRUE)
      print(i)
    }
    colnames(geData) <- colNames
  }


  setwd(FolderWd)
  return(geData)



}

#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapuuidbarcode <- function(uuids){
  # Using tcha api: https://goo.gl/M1uQLR
  uuids <- paste0(uuids, collapse = ",")
  ans <- fromJSON(postForm("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch",
                           .opts = list(postfields = uuids,
                                        httpheader = c('Content-Type' = 'text/plain'),
                                        ssl.verifypeer = FALSE)))

  # transform to dataframe
  x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame))[[1]])
  barcodes <- seq(1,length(x),2)
  uuid <- seq(2,length(x),2)
  x <- data.frame(x[uuid],as.character(x[barcodes]),
                  stringsAsFactors = FALSE)
  colnames(x) <- c("uuid","barcode")
  return(x)
}

# Get sdrf file/array_design of a line
# example
# query <- TCGAQuery(tumor = "BRCA")
# getMagecontent(query[1,])
# Obs: delete the file after reading
#      is it better to save it?
getMage <- function(line){
  tcga.db <- get("tcga.db")
  root <- "https://tcga-data.nci.nih.gov"
  mages <-  tcga.db[grep("mage-tab",tcga.db$name),]
  # get the mage file for the entry
  mage <- subset(mages, mages$Disease == line$Disease &
                   mages$Platform == line$Platform &
                   mages$Center == line$Center)

  if (dim(mage)[1] != 0) {
    file <- basename(mage$deployLocation)
    if ( !file.exists(file)) {
      download(paste0(root,mage$deployLocation), file, quiet = TRUE)
    }
    folder <- gsub(".tar.gz","",file)
    if ( !file.exists(folder)) {
      untar(file)
    }
    files <- list.files(folder)
    if (line$Platform == "MDA_RPPA_Core") {
      sdrf <- files[grep("array_design",files)]
    } else {
      # Platform is not MDA_RPPA_Core
      sdrf <- files[grep("sdrf",files)]
    }
    if (length(sdrf) > 1) {
      sdrf <- sdrf[1]
    }
    df <- read.delim(file = file.path(folder,sdrf),
                     sep = "\t",
                     stringsAsFactors = FALSE,
                     fileEncoding = "latin1")
  }
  unlink(file)
  unlink(folder, recursive = TRUE)
  return(df)
}

#' @title TCGA TCGAPrepare
#' @description
#'  Prepare data matrices for downstream analysis,
#'  ready to use also with other R packages
#'  DNA Methylation
#'   Row: matrix with genes/loci
#'   Cols: samples in columns
#'
#' @export
TCGAPrepare2 <- function(query, dir = NULL, type = NULL){

  if (is.null(dir)) {
    message("Argument dir is NULL. Plese provide the directory
            with the folders to be prepared. ")
    return(NULL)
  }
  if (length(unique(query$Platform)) > 1 | length(unique(query$Center)) > 2) {
    message("Sorry! But, for the moment, we can only prepare on type of
            platform per call")
    return(NULL)
  } else {
    platform <- unique(query$Platform)
  }

  files <- NULL
  dirs <- gsub(".tar.gz","",basename(query$deployLocation))
  for (i in seq_along(dirs)) {
    aux <- list.files(file.path(dir,dirs[i]), full.names = TRUE, recursive = TRUE)
    files <- c(files, aux )
  }
  files <- files[-grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE",files)]

  df <- NULL
  if (grepl("humanmethylation",tolower(platform))) {
    for (i in seq_along(files)) {
      data <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      sample <- gsub("\\.", "-", colnames(data)[2])
      colnames(data) <- data[1,]
      data <- data[-1,] # removing Composite Element REF
      colnames(data)[2] <- sample
      if (i == 1) {
        df <- data[, c(1, 3:5, 2)]
      } else {
        df <- merge(df, data[, c(1, 2)],
                    by = "Composite Element REF")
      }
    }
    rownames(df) <- df$Composite.Element.REF
    message("Removing X Y chromossomes")
    df <- df[df$Chromosome != "X" & df$Chromosome != "Y", ]
    # methylation$Chromosome <- NULL

    # remove NA lines
    message("Removing NA Lines")
    df <- na.omit(df)

  }

  if (grepl("mda_rppa_core",tolower(platform))) {
    for (i in seq_along(files)) {
      data <- read.table(files[i], header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
      sample <- gsub("\\.", "-", colnames(data)[2])
      colnames(data) <- data[1,]
      data <- data[-1,] # removing Composite Element REF
      colnames(data)[2] <- sample
      print(files[i])

      print(sample)
      if (i == 1) {
        df <- data
      } else {
        df <- merge(df, data,by = "Composite Element REF")
      }
    }

    # get array_design.txt from mage folder
    # and change uuid by Barcode
    uuid <- colnames(df)
    idx <- grep("Sample|Control",uuid)
    uuid <- uuid[-idx]
    map <- mapuuidbarcode(uuid)
    idx <- which(colnames(df) %in% map$uuid)
    colnames(df)[idx] <- map$barcode

  }
  # case: header has barcode
  # Line 2 is useless
  if (grepl("agilent",tolower(platform))) {
    for (i in seq_along(files)) {
      data <- read.table(files[i], header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
      data <- data[-1,] # removing Composite Element REF
      if (i == 1) {
        df <- data
      } else {
        df <- merge(df, data,by = "Hybridization REF")
      }
    }
  }

  if (grepl("illuminaga",tolower(platform))) {
    for (i in seq_along(files)) {
      data <- read.table(files[i], header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE,
                         comment.char = "#",fill = TRUE)
      data <- data[-1,] # removing Composite Element REF
      if (i == 1) {
        df <- data
      } else {
        df <- merge(df, data,by = "Hybridization REF")
      }
    }
  }
  return(df)
}
