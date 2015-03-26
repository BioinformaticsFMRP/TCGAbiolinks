#' @title TCGA Version Detailed
#'
#' @description  TCGA Version Detailed
#'
#' @param tumor tumor code between "acc"  "blca" "brca" "cesc" "chol" "cntl" "coad" "dlbc" "esca" "fppp" "gbm"
#'                                 "hnsc" "kich" "kirc" "kirp" "laml" "lcml" "lgg"  "lihc" "lnnh" "luad" "lusc"
#'                                 "meso" "misc" "ov"   "paad" "pcpg" "prad" "read" "sarc" "skcm" "stad" "tgct"
#'                                 "thca" "thym" "ucec" "ucs"  "uvm"
#'
#' @param centerType type code between "bcr"  "cgcc" "gsc"
#'
#' @param center center code between "biotab"                  "nationwidechildrens.org"
#'                                   "bcgsc.ca"                "broad.mit.edu"
#'                                   "jhu-usc.edu"             "mdanderson.org"
#'                                   "unc.edu"                 "hgsc.bcm.edu"
#'                                   "hms.harvard.edu"         "genome.wustl.edu"
#'                                   "ucsc.edu"                "intgen.org"
#'                                   "hudsonalpha.org"         "lbl.gov"
#'                                   "mskcc.org"               "supplemental"
#'
#' @param platform platform code between  "clin"                                "bio"
#'                                        "biotab"                              "diagnostic_images"
#'                                        "pathology_reports"                   "tissue_images"
#'                                        "illuminahiseq_mirnaseq"              "genome_wide_snp_6"
#'                                        "humanmethylation450"                 "mda_rppa_core"
#'                                        "illuminahiseq_rnaseqv2"              "illuminaga_dnaseq_curated"
#'                                        "illuminaga_dnaseq_automated"         "illuminaga_dnaseq_cont_automated"
#'                                        "mixed_dnaseq_curated"                "illuminaga_mirnaseq"
#'                                        "illuminahiseq_dnaseqc"               "illuminahiseq_wgbs"
#'                                        "illuminahiseq_rnaseq"                "illuminahiseq_totalrnaseqv2"
#'                                        "illuminaga_dnaseq"                   "humanmethylation27"
#'                                        "agilentg4502a_07_3"                  "illuminahiseq_dnaseq_automated"
#'                                        "illuminahiseq_dnaseq_cont_automated" "miRNASeq"
#'                                        "microsat_i"                          "illuminaga_rnaseq"
#'                                        "illuminaga_rnaseqv2"                 "RNASeq"
#'                                        "solid_dnaseq"                        "mixed_dnaseq_automated"
#'                                        "minbio"                              "abi"
#'                                        "ht_hg-u133a"                         "hg-cgh-244a"
#'                                        "hg-cgh-415k_g4124a"                  "humanhap550"
#'                                        "illuminadnamethylation_oma002_cpi"   "illuminadnamethylation_oma003_cpi"
#'                                        "huex-1_0-st-v2"                      "agilentg4502a_07_1"
#'                                        "agilentg4502a_07_2"                  "h-mirna_8x15k"
#'                                        "h-mirna_8x15kv2"                     "mixed_dnaseq_cont"
#'                                        "mixed_dnaseq"                        "mixed_dnaseq_cont_curated"
#'                                        "hg-u133_plus_2"                      " "
#'                                        "human1mduo"                          "cgh-1x1m_g4447a"
#'                                        "illuminaga_mrna_dge"                 "solid_dnaseq_curated"
#'
#' @param level level 1 2 3
#'
#'@param barcode Get for each version all barcodes?
#'               Default FALSE
#'               Takes a lot of time to get barcode
#' @param version version code -TO DO-
#'
#' @param file - link reference data matrix
#'
#' @param qOutput place where the query is saved to be downloaded automatically.
#'        The folder can be specified in both TCGAQuery and TCGADownload
#'
#' @author Tiago
#' @import stringr
#' @export
TCGAVersion <- function(tumor = "all",
                        centerType = "all",
                        center = "all",
                        platform = "all",
                        level = "all",
                        version = "all",
                        barcode=F,
                        file = system.file("extdata/dataFolders.rda",
                                           package="TCGAbiolinks"),
                        qOutput = "data/version/"){

  if(!exists("dataFolders")) load(file)
  dataFolders <- get("dataFolders", envir=environment())

  info.tcga <- get.data.folder(tumor,centerType,center,platform)
  magetab <- info.tcga[grep("mage-tab", info.tcga[,"Folder"]),]
  data    <- info.tcga[grep("Level_",   info.tcga[,"Folder"]),]

  # Search resulted in only one result
  if(is.null(nrow(data))){
    platform.url <- dirname(data["Manifest"])
  } else {
    # Search resulted in more than one result
    platform.url <- unique(dirname(data[,"Manifest"]))
  }
  if(is.null(nrow(magetab))){
    magetab.url <- dirname(magetab["Manifest"])
  } else {
    # Search resulted in more than one result
    magetab.url <- unique(dirname(magetab[,"Manifest"]))
  }

  if(length(platform.url) == 0) {
    message("No results found")
    return (NULL)
  } else{
    message(paste("Found", length(platform.url), "Version of", platform, sep = " "))
    message("Looking for metadata...")
  }

  dir.create(qOutput, showWarnings = F)
  version  <- as.data.frame(matrix(0,length(platform.url),9))
  colnames(version) <- c("Version","Disease","Platform","Level",
                         "Batch","Date","Samples","Total Size","Files")


  for(j in 1:length(platform.url)){

    message(paste("Version", j , "of", length(platform.url),
                  basename(platform.url[j]),
                  sep=" "))

    content <- DownloadHTML(platform.url[j])
    content <- as.matrix(unlist(strsplit(content, "  ")))
    content <- str_trim(content[content != ""])

    time  <- unique(content[grep(":",content)])
    # handle difference in minutes
    if(length(time)>1) time <- time[1]

    # TODO: size should be filtered - not all files are relevant
    # DONE: filtering happens in the getTotalSize function

    regex <- "^[0-9]+\\.?[0-9]*([K]{1}|[M]{1}|[G]{1})"
    sizes <- content[grep(regex,content)]
    info  <- rev(unlist(strsplit(basename(platform.url[j]), "\\.")))

    version$Disease[j]  <- unlist(strsplit(info[6], "_"))[2]
    version$Platform[j] <- info[5]
    version$Level[j]    <- info[4]
    version$Batch[j]    <- info[3]
    version$Date[j]     <- time
    version$Version[j]  <- basename(platform.url[j])
    version$SizeMB[j]   <- getTotalSize(sizes)
    version$Files[j]  <-   get.manifest.files(platform.url[j],qOutput)
    version$Samples[j]  <- length(sizes)

    if(barcode){
      version$Barcodes[j] <- get.barcodes(platform=version$Platform[j],
                                          disease=version$Disease[j],
                                          center=unlist(strsplit(info[6], "_"))[1]
                                          )
    }
    message('...done')
  }
  return(version)
}

get.manifest.files <- function(platform.url,qOutput){
  manifest <- paste0(platform.url,"/MANIFEST.txt")
  if(RCurl::url.exists(manifest)){
    download(manifest,
             destfile = paste0(qOutput,"filenames.txt"),
             mode="w",
             quiet = 1)
    manifest.content <- read.table(file = paste0(qOutput,"filenames.txt"),sep="")
    manifest.files <- as.character(manifest.content$V2)
    unlink(paste0(qOutput,"filenames.txt"))
    return (list(manifest.files))
  }
}

#' @title Get bar code info
#'
#' @description  Get bar code info using sdrf files
#'               Obs1: it takes a lot of time, is there a better way
#'               than downloading sdrf?
#'               Obs2: Not all folders has mage-tab file
#'               Solution 1: search for barcode pattern in the filenames
#'               Solution 2: maf files has it inside. Should use downloaded files.
#'
#' @param magetab.url path to mage-tab folder
#' @keywords internal
#'
get.barcodes <- function(platform,disease,center){

  idx <- which(grepl(platform,tcga.barccodes$Platform) &
                 tcga.barccodes$Disease == disease &
                 grepl(center,tcga.barccodes$Receiving.Center))

  return (list(unique(as.character(tcga.barccodes[idx,]$Barcode))))
}

get.data.folder <- function(tumor,centerType,center,platform){

  ifelse(tumor != "all",x <- subset(dataFolders, dataFolders[,"Tumor"] == tolower(tumor)),x<-dataFolders)

  if(centerType != "all" && is.null(nrow(x))) x <- subset(x, x["CenterType"] == tolower(centerType))
  if(centerType != "all" && !is.null(nrow(x))) x <- subset(x, x[,"CenterType"] == tolower(centerType))

  if(center != "all" && is.null(nrow(x)))  x <- subset(x, x["Center"] == tolower(center))
  if(center != "all" && !is.null(nrow(x)))  x <- subset(x, x[,"Center"] == tolower(center))

  if(platform != "all" && is.null(nrow(x)))  x <- subset(x, x["Platform"] == tolower(platform))
  if(platform != "all" && !is.null(nrow(x))) x <- subset(x, x[,"Platform"] == tolower(platform))

  return (x)
}
# param sizeList - list of sizes in KB, MB
# output total size in MB
getTotalSize <- function(sizeList){
  sizeK  <- sizeList[grep("K",sizeList)]
  sizeK  <- sizeK[as.numeric(gsub("K","",sizeK))>5] #choose a proper limit
  sizeM  <- sizeList[grep("M",sizeList)]
  sizeG  <- sizeList[grep("G",sizeList)]
  totalK <- round(sum(as.numeric(gsub("K","",sizeK))) / 1000)
  totalM <- sum(as.numeric(gsub("M","",sizeM)))
  totalG <- sum(as.numeric(gsub("G","",sizeG))) * 1000

  return(paste0((totalM + totalK + totalG), "Mb"))
}
