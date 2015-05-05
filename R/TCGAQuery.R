#' @title TCGA query
#'
#' @description  Crossfinding of file locations for downloading (TCGADownload)
#' @param tumor tumor code between "acc"  "blca" "brca" "cesc" "chol" "cntl" "coad" "dlbc" "esca" "fppp" "gbm"
#'                                 "hnsc" "kich" "kirc" "kirp" "laml" "lcml" "lgg"  "lihc" "lnnh" "luad" "lusc"
#'                                 "meso" "misc" "ov"   "paad" "pcpg" "prad" "read" "sarc" "skcm" "stad" "tgct"
#'                                 "thca" "thym" "ucec" "ucs"  "uvm"
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
#' @param added.since Date to search for files (since)
#' @param added.up.to Date to search for files (up.to)
#' @param level level 1 2 3
#' @param listSample list of barcodes to be considered in the search
#' @example inst/examples/queryExamples.R
#' @seealso TCGADownload
#' @export
#' @import downloader
TCGAQuery <- function(tumor=NULL,
                      platform=NULL,
                      added.since=NULL,
                      added.up.to=NULL,
                      listSample=NULL,
                      level=NULL
){
  if(!(is.element(tumor,disease.table$abbreviation) | is.null(tumor))){
    message("Disease not found. Chosse between:")
    message(paste(disease.table$abbreviation, collapse = " "))
    stop("Invalid tumor")
  }
  if(!(is.element(platform,platform.table$alias) | is.null(platform))){
    message("Platform not found. Chosse between:")
    message(paste(platform.table$alias, collapse = " "))
    stop("Invalid platform")
  }

  if(!(is.element(level,c("1","2","3")) | is.null(level))){
    message("Levelnot found. Chosse between:'1', '2' or '3'")
    stop("Invalid platform")
  }

  if(!is.null(added.since)){
    d <- try( as.Date( added.since, format= "%m/%d/%Y" ) )
    if( class( d ) == "try-error" || is.na( d ) ) print( "Date format should be mm/dd/YYYY" )
  }
  if(!is.null(added.up.to)){
    d <- try( as.Date( added.up.to, format= "%m/%d/%Y" ) )
    if( class( d ) == "try-error" || is.na( d ) ) print( "Date format should be mm/dd/YYYY" )
  }

  # to be improved
  if(!is.null(listSample)){
    #archives <- c()
    files <- c()
    for(i in seq_along(listSample)){
      # table with barcode id
      # example: query=BiospecimenBarcode[@barcode=TCGA-28-2499*]
      message("Searching for barcode files...")
      message(paste("Barcode:",listSample[i]))
      db <- get.barcode.table(listSample[i])

      # Improvement: using portion analyte in order to select platform
      if(!is.null(platform)){
        analyte <- c()
        pat <- "TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-[[:alnum:]]{3}-[[:alnum:]]{2}"
        if(is.element(tolower(platform),tolower(dna.plat))){
          analyte = c(analyte,"D")
        }
        if(is.element(tolower(platform),tolower(rna.plat))){
          analyte = c(analyte,"R")
        }
        if(is.element(tolower(platform),tolower(total.rna.plat))){
          analyte = c(analyte,"T")
        }
        if(is.element(tolower(platform),tolower(wgarubcon.plat))){
          analyte = c(analyte,"G")
        }
        if(is.element(tolower(platform),tolower(mirna.plat))){
          analyte = c(analyte,"H")
        }
        if(is.element(tolower(platform),tolower(wgaqiagen2.plat))){
          analyte = c(analyte,"X")
        }
        if(is.element(tolower(platform),tolower(wgaqiagen1.plat))){
          analyte = c(analyte,"W")
        }
        idx <- c()
        for(i in seq_along(analyte)){
          aux <- grep(paste0(pat,"[",analyte[i],"]"),as.character(db$barcode))
          idx <- union(idx,aux)
        }
        print(length(idx))
        db <- db[idx,]
      }
      if(!is.null(platform)){
        pat <- "TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-"
        idx <- grep(platform,plat.center$platform)
        centers <- as.character(plat.center[idx,"center"])
        idx <- c()
        for(i in seq_along(centers)){
          aux <- grep(paste0(pat,centers[i]),as.character(db$barcode))
          idx <- union(idx,aux)
        }
        print(length(idx))
        db <- db[idx,]
      }

      # get getBcrArchiveCollection table
      for(i in seq_along(db$id)){
        #aux <- getBcrArchiveCollection(db[i,"id"])
        #archives <- rbind(archives,aux)
        aux <- get.samples.files(db[i,"id"])
        aux$barcode <- db[i,"barcode"]
        if(!is.null(aux)){
          files <- rbind(files,aux)
        }
      }
    }
    x <- subset(files, files$isLatest == 1)
    x <- x[!duplicated(x), ]
    x <- x[,order(names(x))]
    if(!is.null(platform)){
      x <- subset(x, grepl(tolower(platform),tolower(x$Platform)))
    }
    if(!is.null(platform)){
      x <- subset(x, tolower(Disease) == tolower(tumor))
    }
    if(!is.null(level)){
      x <- subset(x, grepl(paste0("Level_",level),name))
    }

  }
  else{
    message("CREATING TABLE")
    x <- create.tcga.table (platform=platform,type=level,disease=tumor)
  }
  if(!is.null(added.since)){
    x <- subset(x, as.Date(addedDate) > as.Date(added.since,"%m/%d/%Y"))
  }
  if(!is.null(added.up.to)){
    x <- subset(x, as.Date(addedDate) < as.Date(added.up.to,"%m/%d/%Y"))
  }

  return(x)
}
