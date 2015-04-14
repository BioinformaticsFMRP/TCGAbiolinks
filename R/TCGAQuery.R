#' @title TCGA query
#'
#' @description  Crossfinding of file locations for downloading (TCGADownload)
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
#' @param version version code -TO DO-
#'
#' @param i - interactive -TO DO-
#'
#' @param file - link reference data matrix
#'
#' @param qOutput place where the query is saved to be downloaded automatically.
#'        The folder can be specified in both TCGAQuery and TCGADownload
#'
#' @examples
#' \dontrun{
#'   TCGAQuery(tumor = "all",centerType = "all",center = "all",
#'             platform = "all",level = "all",version = "all",
#'             i = F,file = "data/dataFolders.rda",
#'             qOutput = "data/query/")
#' }
#'
#' @author Davide Garolini
#' @seealso TCGADownload
#' @export
#' @import downloader
#' @import stringr
TCGAQuery <- function(tumor = "all",
                      centerType = "all",
                      center = "all",
                      platform = "all",
                      level = "all",
                      version = "all",
                      i = F,
                      metadata = F,
                      file = system.file("extdata/dataFolders.rda",
                                         package="TCGAbiolinks"),
                      qOutput = "data/query/"){

  # laod tcga info if not done
  if(!exists("dataFolders")) load(file)
  dataFolders <- get("dataFolders", envir=environment())

  dir.create(path = qOutput, showWarnings = F, recursive = T)

  if(!i){
    ifelse(tumor != "all",x <- subset(dataFolders, dataFolders[,"Tumor"] == tolower(tumor)),x<-dataFolders)

    if(centerType != "all" && is.null(nrow(x))) x <- subset(x, x["CenterType"] == tolower(centerType))
    if(centerType != "all" && !is.null(nrow(x))) x <- subset(x, x[,"CenterType"] == tolower(centerType))

    if(center != "all" && is.null(nrow(x)))  x <- subset(x, x["Center"] == tolower(center))
    if(center != "all" && !is.null(nrow(x)))  x <- subset(x, x[,"Center"] == tolower(center))


    if(platform != "all" && is.null(nrow(x)))  x <- subset(x, x["Platform"] == tolower(platform))
    if(platform != "all" && !is.null(nrow(x))) x <- subset(x, x[,"Platform"] == tolower(platform))

    if(metadata){
      if(!is.null(nrow(x))){
        magetab <- x[grep("mage-tab", x[,"Folder"]),]
      }else{
        magetab <- x[grep("mage-tab", x["Folder"]),]
      }
    }

    if(!is.null(nrow(x))){
      x <- x[grep("Level_",   x[,"Folder"]),]
    }else{
      x <- x[grep("Level_",   x["Folder"]),]
    }


    if(level != "all"){
      if(is.null(nrow(x))){
        l <- grep(paste("Level_", as.character(level),sep=""), x["Folder"])
        if(length(l)==0) x <- NULL
      }else{
        l <- grep(paste("Level_", as.character(level),sep=""), x[,"Folder"])
        x <- x[l,]
      }
    }
    #     if(version != "all"){ #something smarter could be done (maybe "lastUp" as keyword)
    #       if(is.null(nrow(x)))l <- grep(paste("Level_",as.character(version),sep=""),x["Folder"])
    #       if(!is.null(nrow(x)))l <- grep(paste("Level_",as.character(version),sep=""),x[,"Folder"])
    #       if(length(l)>1) x<-x[l,]
    #       if(length(l)==0)x<-NULL
    #     }

    if(length(x)==0){
      stop("Nothing found. Check the proper spelling in the documentation.")
    }else if(is.null(nrow(x))) {
      print("Found: 1 folder. Start downloading filenames:")
    }else{
      print(paste("Found:", length(x[,1]), "folders. Preparing filenames:",sep=" "))
    }



    if(metadata && length(magetab)!=0 ) {
      print("Found magetab metadata. Downloading..")
      for(k in 1:length(magetab[,"Manifest"])){
        if(RCurl::url.exists(dirname(magetab[,"Manifest"][k]))){
          download(dirname(magetab[,"Manifest"][k]),
                   "temp.html",
                   mode="wb",
                   quiet = 1)
          tmp <- htmlTreeParse("temp.html")
          mt <- capture.output(tmp)
          unlink("temp.html")
          rm(tmp)
        }else{
          stop("Can't find URL. Please check the web site or the internet connection.")
        }
        mt <- as.matrix(unlist(strsplit(mt, "  ")))
        mt <- str_trim(mt[mt != ""])
        mt <- mt[grep(".sdrf.txt",mt)]
        mt <- sapply(strsplit(mt, ">"), function(y) y[1])
        mt <- sapply(strsplit(mt, "="), function(y) y[2])
        mt <- sub('\"','',mt)
        mt <- sub('\"','',mt)
        mt <- paste(dirname(magetab[,"Manifest"][k]),mt, sep = "/")
        download(mt,
                 destfile = paste0(qOutput,
                                   paste0(unlist(strsplit(dirname(magetab[,"Manifest"][k]), split = '/', fixed = TRUE))[14],
                                          ".sdrf.txt")),
                 mode = "w",
                 quiet = 1)
        print(paste("Downloaded",k," metadata out of",length(magetab[,"Manifest"]), sep = " "))
      }
    }else if(metadata){
      #.maf handling

    }




    queryURI = NULL
    if(is.null(nrow(x))){
      fails <- 0
      if(RCurl::url.exists(x["Manifest"])){
        download(x["Manifest"],
                 destfile = paste0(qOutput,"/filenames.txt"),
                 mode="w",
                 quiet = 1)
        queryURI <- paste(unlist(strsplit(x["Manifest"], split='MANIFEST.txt', fixed=TRUE)),
                          as.character(read.table(file = paste0(qOutput,"/filenames.txt"))[2]$V2),sep="")
        print("Donwloaded.")
      }else{
        message(paste0("ERROR: ", x["Manifest"]," not found"))
        fails <- 1
      }
    }else{
      fails <- 0
      for(j in 1:length(x[,"Tumor"])){

        if(RCurl::url.exists(x[,"Manifest"][j])){
          download(x[,"Manifest"][j],
                 destfile = paste0(qOutput,"/filenames.txt"),
                 mode="w",
                 quiet = 1) #character. The mode with which to write the file.
                            #Useful values are "w", "wb" (binary), "a" (append) and "ab".
                            #Only used for the "internal" method.

        print(paste("Downloaded:",j,"out of",length(x[,"Tumor"]),sep=" "))

        queryURI<-c(queryURI,paste(unlist(strsplit(x[,"Manifest"][j], split='MANIFEST.txt', fixed=TRUE)),
                                   as.character(read.table(file = paste(qOutput,"/filenames.txt",sep=""))[2]$V2),sep=""))
        unlink(paste0(qOutput,"/filenames.txt"))
        }
        else{
          message(paste0("ERROR: ", x[,"Manifest"][j] ," not found"))
          fails <- fails + 1
        }
      }
    }
  }
  print(paste("We found",length(queryURI) - fails ,"files",sep=" "))

  if(metadata && length(magetab)!=0){
    print("Returning metadata from .sdrf files")
    n<-list.files(path = qOutput)
    n<-n[grep(".sdrf",n)]
    u <- NULL
    for(i in 1:length(n)){
      if(!is.null(u)) u<-rbind(u,read.delim(paste(qOutput,n[i],sep= "/"),sep = "\t"))
      else u<-read.delim(paste(qOutput,n[i],sep= "/"),sep = "\t")
      #unlink(paste(qOutput,n[i],sep= "/"))
    }
  }else if(metadata){
    print("Returning metadata from .maf files")
  }
  save(queryURI, file = paste0(qOutput,"/fileURLs.rda"))
  #todo - add the showing of the result in a human readable way
  #if(metadata) return(u)
}
