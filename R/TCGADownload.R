#' @title TCGA Download
#'
#' @description  
#'   Download data of samples from TCGA specified by tumor,type,
#'   species and platform type.
#'
#' @details
#'   TCGA download retrieves and stores data of samples belonging 
#'   to the specified cancer type and measured by the specified assay platform.      
#'
#' @author 
#' Antonio Colaprico \email{antonio.colaprico@@gmail.com},
#' Luciano Garofano \email{lucianogarofano88@@gmail.com}, 
#' Claudia Cava \email{claud.cava@@gmail.com},
#' Gianluca Bontempi,
#' Michele Ceccarelli
#'
#' Maintainer: Antonio Colaprico \email{antonio.colaprico@@gmail.com}
#'
#' @param Tumor  a character string indicating the cancer type for 
#'        which to download data. Options include ACC, BLCA, BRCA, 
#'        CESC, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, 
#'        LGG, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, SARC, SKCM, STAD, 
#'        THCA, UCEC, UCS. Look at https://tcga-data.nci.nih.gov/tcga/ 
#'        for Available Cancer Types. 
#' @param Type "mRNA"
#' @param Species "RNASeq" 
#' @param PlatformAndAssociatedData data frame 615 observations of 12 variables,
#'        indicating the different characteristics of the data
#'        e.g. tumour, type, species.  
#' @param downloadFolder = ""
#' @param PlatformType "illuminahiseq_rnaseq"
#' @param nsample number of samples to be analyzed, default all samples
#' @param listSample  a character that indicates TCGA barcodes of the patients/samples to analyze
#' @return
#' TCGADownload returns a matrix with Samples in columns with barcode and probeID in rows.
#' ProbeID and values changes with platforms.
#' \describe{
#' \item{humanmethylation27}{If platform is humanmethylation27 returns values as Hybridization.RE}
#' \item{illuminahiseq_rnaseq}{If platform is illuminahiseq_rnaseq returns values as raw_counts}
#' \item{agilentg4502a_07_3}{If platform is agilentg4502a_07_3 returns values as Hybridization.REF}
#' \item{illuminahiseq_rnaseqv2}{If platform is illuminahiseq_rnaseqv2 returns values as raw_counts}
#' \item{humanmethylation450}{If platform is humanmethylation450 returns values as Hybridization.REF}
#' \item{illuminaga_mirnaseq}{If platform is illuminaga_mirnaseq returns values as read_count}
#' \item{illuminahiseq_mirnaseq}{If platform is illuminahiseq_mirnaseq returns values as read_count}
#' \item{genome_wide_snp_6}{If platform is genome_wide_snp_6 returns values as Segment_Mean}
#' \item{mda_rppa_core}{If platform is mda_rppa_core returns values as Protein Expression}
#' }
#' @example inst/examples/example_TCGADownload.R
#' @import RCurl httr bitops
#' @export


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
    x <- x[-1]
    x <- matrix(unlist(x), nrow = length(x), byrow = T)
    colnames(x) <- x[1, ]
    x <- x[-1, ]
    
    return(x)
  }
  
  
  lstFileSdrf <- list.files(file.path(sdrfFolder))
  lstFileSdrf_plt <- lstFileSdrf[grep(tolower(PlatformType), tolower(lstFileSdrf))]
  listSample_fromSdrf <- read.delim(paste( sdrfFolder,  lstFileSdrf_plt,sep=""))
  
  
  if(PlatformType == "illuminadnamethylation_oma003_cpi" | PlatformType == "illuminadnamethylation_oma002_cpi"){
    toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_2", KeyGrep2 = "MANIFEST.txt") 
    toDdl <- paste(Description, key2a, toDdl, sep = "")
  }else{
    toDdl <- .DownloaDmageTAB_sdrf(Description, keySpecies = key2a, KeyGrep1 = "Level_3", KeyGrep2 = "MANIFEST.txt")
    toDdl <- paste(Description, key2a, toDdl, sep = "")
  }
  
  x <- .DownloadURL(toDdl)
  x <- sapply(strsplit(x, "  "), function(y) y[2])
  
  if(PlatformType == "illuminahiseq_rnaseq"){  x <- x[grep("gene.quantification", x)] }
  if(PlatformType == "agilentg4502a_07_3" | PlatformType == "agilentg4502a_07_2" | PlatformType == "agilentg4502a_07_1"){    x <- x[grep("tcga_level3", x)]}
  if(PlatformType == "ht_hg-u133a" | PlatformType == "h-mirna_8x15kv2" | PlatformType == "h-mirna_8x15k"){    x <- x[grep("level3", x)]}
  if(PlatformType == "hg-u133_plus_2"){    x <- x[grep("pergene", x)]}
  if(PlatformType == "illuminahiseq_rnaseqv2" || PlatformType == "illuminahiseq_totalrnaseqv2"){ x <- x[grep("rsem.genes.results", x)] }
  if(PlatformType == "humanmethylation27"){ x <- x[grep("HumanMethylation27", x)] }
  if(PlatformType == "humanmethylation450"){ x <- x[grep("HumanMethylation450", x)] }
  if(PlatformType == "illuminadnamethylation_oma003_cpi"){ x <- x[grep("IlluminaDNAMethylation_OMA003_CPI", x)] }
  if(PlatformType == "illuminadnamethylation_oma002_cpi"){ x <- x[grep("IlluminaDNAMethylation_OMA002_CPI", x)] }
  if(PlatformType == "illuminaga_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }
  if(PlatformType == "illuminahiseq_mirnaseq"){ x <- x[grep("mirna.quantification", x)] }
  
  if(PlatformType == "genome_wide_snp_6"){ x <- x[grep("hg19.seg", x)]
                                           x <- x[-grep("nocnv", x)]}
  if(PlatformType == "illuminahiseq_dnaseqc" | PlatformType == "hg-cgh-415k_g4124a" | PlatformType == "hg-cgh-244a"){ x <- x[grep("Segment", x)] }
  if(PlatformType == "humanhap550"){ x <- x[grep("seg.txt", x)] }
  if(PlatformType == "human1mduo"){ x <- x[grep("seg.txt", x)] }
  
  if(PlatformType == "huex-1_0-st-v2"){  x <- x[grep("gene", x)] }
  
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
    if(PlatformType == "huex-1_0-st-v2" | PlatformType == "ht_hg-u133a" | PlatformType == "hg-u133_plus_2" | PlatformType == "h-mirna_8x15kv2" | PlatformType == "h-mirna_8x15k" |
         PlatformType == "human1mduo" | PlatformType == "humanhap550" | PlatformType == "illuminadnamethylation_oma003_cpi" | PlatformType == "illuminadnamethylation_oma002_cpi"){
      tmp2 <- .DownloadData_fromURL(url =samplesList[i] ,sep="\t", header = T)
    }else{
      tmp2 <- .DownloadData_fromURL(url =samplesList[i] ,sep="", header = T)
    }
    #filename <- paste(downloadFolder, x[i], sep = "")
    filename <- x[i]
    write.csv(tmp2, filename)
    print(paste(x[i], " ... sample n. ", i, " of ", length(samplesList), sep = ""))
  }

  
}