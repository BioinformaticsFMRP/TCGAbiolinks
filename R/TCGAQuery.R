getBarcode <- function(table){
  tcga.db <- get("tcga.db")
  root <- "https://tcga-data.nci.nih.gov"
  allFiles <- c()
  mages <-  tcga.db[grep("mage-tab",tcga.db$name),]
  max <- length(table[,1])
  pb <- txtProgressBar(min = 0, max = max , style = 3)

  for (i in seq_along(table[,1])) {
    message(i," - ",table[i,]$deployLocation)
    if (table[i,]$deployStatus == "Available") {

      # get the mage file for the entry
      mage <- subset(mages, mages$Disease == table[i,]$Disease &
                       mages$Platform == table[i,]$Platform &
                       mages$Center == table[i,]$Center)

      # if no mage file: not found (case PLatform:ABI)
      if (dim(mage)[1] == 0) {
        # For IlluminaGA and ABI platforms
        # Barcodes are in maf files
        if ((grepl("DNASeq",table[i,]$Platform) &
               grepl("Level",table[i,]$name))
        ) {
          barcodes <- tcga.get.barcode(table[i,])

          if (is.null(barcodes)){
            table[i,]$deployStatus <- "Not found"
            next
          }
          table[i,]$deployStatus <- barcodes
          setTxtProgressBar(pb, i)
        } else if (table[i,]$Platform == "diagnostic_images" |
                     table[i,]$Platform == "pathology_reports" |
                     table[i,]$Platform == "tissue_images"
        ) {
          print("image or report")

          folder <- gsub(".tar.gz","",table[i,]$deployLocation)
          files <- getFileNames(paste0(root,folder))
          if (is.null(files))
            next
          files <- files[grepl("TCGA",files)]
          if (table[i,]$Platform == "pathology_reports") {
            barcode <- substr(files,1,12)
          } else {
            barcode <- substr(files,1,23)
          }
          table[i,]$deployStatus <- paste0(unique(barcode), collapse = ",")
        } else if (table[i,]$Platform == "ABI") {
          print("ABI")
          folder <- gsub(".tar.gz","",table[i,]$deployLocation)
          files <- getFileNames(paste0(root,folder))
          if (is.null(files))
            next

          maf <- files[grep("somatic.maf",files)]
          print(maf)
          if(length(maf) == 0){
            table[i,]$deployStatus <- "Not found"
            next
          }
          if ( !file.exists(maf)) {
            download(paste0(root,folder,"/",maf), maf, quiet = TRUE)
          }
          df <- read.delim(file = maf,
                           sep = "\t", comment.char = "#",
                           stringsAsFactors = FALSE,
                           fileEncoding = "latin1",
                           row.names = NULL)
          barcode <- union(df$TUMOR_SAMPLE_ID,
                           df$MATCH_NORM_SAMPLE_ID)
          if (is.null(barcode)) {
            barcode <- union(df$Tumor_Sample_Barcode,
                             df$Matched_Norm_Sample_Barcode)
            table[i,]$deployStatus <- paste0(barcode, collapse = ",")
            unlink(maf)
          } else {
            table[i,]$deployStatus <- "Not found"
          }
        }
        next
      }

      # In case we have two files
      # This should not happen after filtering by center
      # probably this code can be remove until next comment
      file <- basename(mage$deployLocation)
      x <- stringr::str_replace_all(file, "[^[:alnum:]]", "")
      y <- stringr::str_replace_all(table[i,]$baseName, "[^[:alnum:]]", "")
      idx <- grep(y,x)

      if (length(idx) > 0) {
        file <- file[idx]
        mage <- mage[idx,]
      }
      # get name of files to be deleted
      allFiles <- c(allFiles, file)
      print(file)
      if ( !file.exists(file)) {
        download(paste0(root,mage$deployLocation), file, quiet = TRUE)
      }
      folder <- gsub(".tar.gz","",file)
      if ( !file.exists(folder)) {
        untar(file)
      }
      files <- list.files(folder)
      print("has mage-tab")

      # For this platform the barcodes are in array_design files
      if (table[i,]$Platform == "MDA_RPPA_Core") {
        sdrf <- files[grep("array_design",files)]
        # case with 2 array_design BLCA
        if (length(sdrf) > 1) {
          sdrf <- sdrf[1]
        }
        df <- read.delim(file = file.path(folder,sdrf),
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         fileEncoding = "latin1")

        if (is.element("Sample.description",colnames(df))) {
          barcode <- unique(as.character(df$Sample.description))
        } else {
          barcode <- unique(as.character(df$Biospecimen.Barcode))
        }
        table[i,]$deployStatus <- paste0(barcode, collapse = ",")

      } else {
        # Platform is not MDA_RPPA_Core
        sdrf <- files[grep("sdrf",files)]
        df <- read.delim(file = file.path(folder,sdrf),
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         fileEncoding = "latin1")

        # The mage file could is common to the same
        # disease, center and platform.
        # So we're going to fill them so we don't need
        # To re-read the barcodes
        for (j in seq_along(table[,1])) {
          if (table[j,]$Disease == table[i,]$Disease &&
                table[j,]$Platform == table[i,]$Platform &&
                table[j,]$Center == table[i,]$Center) {
            aux <- grep("Comment..TCGA.Archive.Name",colnames(df))

            barcode <- data.frame()
            x <- table[j,]$name

            for (z in seq_along(aux)){
              barcode <- rbind(barcode,subset(df, x == df[,aux[z]]))
            }

            barcode <- unique(as.character(barcode$Comment..TCGA.Barcode.))
            barcode <- barcode[barcode != "->"]
            table[j,]$deployStatus <- paste0(barcode,collapse = ",")
          }
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  setTxtProgressBar(pb, max)
  close(pb)
  colnames(table)[4] <- "barcode"
  # Removing the files mess
  # porbably this can be used up in the code
  unlink(allFiles)
  folders <- gsub(".tar.gz","",allFiles)
  unlink(folders, recursive = TRUE)

  return(table)
}

# ------------------------ Encode search
#' @title TCGAQuery
#' @description
#'    Searches in the tcga database
#' @param tumor Disease Examples:
#' \tabular{lllll}{
#'OV   \tab BRCA \tab CESC \tab ESCA \tab PCPG\cr
#'LUSC \tab LGG  \tab SKCM \tab KICH \tab CHOL\cr
#'GBM  \tab UCEC \tab PRAD \tab PAAD \tab THYM\cr
#'KIRC \tab THCA \tab SARC \tab LAML \tab TGCT\cr
#'COAD \tab KIRP \tab HNSC \tab ACC  \tab UVM \cr
#'READ \tab BLCA \tab DLBC \tab UCS  \tab FPPP\cr
#'LUAD \tab LIHC \tab STAD \tab MESO \tab CNTL
#'}
#' @param platform Example:
#' \tabular{ll}{
#'CGH- 1x1M_G4447A                   \tab IlluminaGA_RNASeqV2   \cr
#'AgilentG4502A_07                  \tab IlluminaGA_mRNA_DGE    \cr
#'Human1MDuo                        \tab HumanMethylation450    \cr
#'HG-CGH-415K_G4124A                \tab IlluminaGA_miRNASeq    \cr
#'HumanHap550                       \tab IlluminaHiSeq_miRNASeq \cr
#'ABI                               \tab H-miRNA_8x15K  \cr
#'HG-CGH-244A                       \tab SOLiD_DNASeq                       \cr
#'IlluminaDNAMethylation_OMA003_CPI \tab IlluminaGA_DNASeq_automated        \cr
#'IlluminaDNAMethylation_OMA002_CPI \tab HG-U133_Plus_2                 \cr
#'HuEx- 1_0-st-v2 \tab Mixed_DNASeq                  \cr
#'H-miRNA_8x15Kv2 \tab IlluminaGA_DNASeq_curated      \cr
#'MDA_RPPA_Core   \tab IlluminaHiSeq_TotalRNASeqV2    \cr
#'HT_HG-U133A     \tab IlluminaHiSeq_DNASeq_automated \cr
#'diagnostic_images                 \tab microsat_i                     \cr
#'IlluminaHiSeq_RNASeq              \tab SOLiD_DNASeq_curated           \cr
#'IlluminaHiSeq_DNASeqC             \tab Mixed_DNASeq_curated           \cr
#'IlluminaGA_RNASeq                 \tab IlluminaGA_DNASeq_Cont_automated   \cr
#'IlluminaGA_DNASeq                 \tab IlluminaHiSeq_WGBS             \cr
#'pathology_reports                 \tab IlluminaHiSeq_DNASeq_Cont_automated\cr
#'Genome_Wide_SNP_6                 \tab bio                            \cr
#'tissue_images                     \tab Mixed_DNASeq_automated         \cr
#'HumanMethylation27                \tab Mixed_DNASeq_Cont_curated      \cr
#'IlluminaHiSeq_RNASeqV2            \tab Mixed_DNASeq_Cont
#'}
#' @param level '1' '2' '3'
#' @param added.since 04- 14-2010
#' @param added.up.to 04- 14-2010
#' @param center center name
#' @param samples List of samples. Ex:c('TCGA-04-06-*','TCGA-04-08-*')
#' @example inst/examples/tcgaSearch.R
#' @export
#' @importFrom downloader download
#' @importFrom knitr kable
#' @return A dataframe with the results of the query
#'        (lastest version of the files)
tcgaQuery <- function(tumor = NULL, platform = NULL, added.since = NULL,
                      added.up.to = NULL, samples = NULL, center = NULL,
                      level = NULL) {

  disease.table   <- get("disease.table")
  platform.table  <- get("platform.table")
  center.table  <- get("center.table")
  db <-  get("tcga.db")

  if (!is.null(tumor)) {
    sapply(tumor, function(x){
      if (!(is.element(tolower(x),
                       tolower(disease.table$abbreviation)))) {
        suppressWarnings(
          df <- as.data.frame(matrix(sort(unique(disease.table$abbreviation)),
                                     ncol = 8))
        )
        print(kable(df, col.names = NULL, format = "pandoc",
                    caption = "TCGA tumors"))
        cat("=======================================================\n")
        cat("ERROR: Disease not found. Select from the table above.\n")
        cat("=======================================================\n")
      }
      return(NULL)
    })
  }
  if (!is.null(platform)) {
    sapply(platform, function(x){
      if (!(is.element(tolower(x), tolower(platform.table$name)))) {
        suppressWarnings(
          df <- as.data.frame(matrix(sort(unique(platform.table$name)),
                                     ncol = 3))
        )
        print(kable(df, col.names = NULL, format = "pandoc",
                    caption = "TCGA Platforms"))
        cat("=======================================================\n")
        cat("ERROR: Platform not found. Select from the table above.\n")
        cat("=======================================================\n")
        return(NULL)
      }
    })}

  if (!is.null(center)) {
    if (!(is.element(tolower(center), tolower(center.table$name)))) {
      suppressWarnings(
        df <- as.data.frame(matrix(sort(unique(center.table$name)),
                                   ncol = 3))
      )
      print(kable(df, col.names = NULL, format = "pandoc",
                  caption = "TCGA Centers"))
      cat("=======================================================\n")
      cat("ERROR: Center not found. Select from the table above.\n")
      cat("=======================================================\n")
      return(NULL)
    }
  }
  if (!is.null(level)) {
    if (!(is.element(level, c("1", "2", "3","mage-tab")))) {
      message("Level not found. Chosse between:'1', '2','3','mage-tab'")
      stop("Invalid platform")
    }
  }

  if (!is.null(added.since)) {
    d <- try(as.Date(added.since, format = "%Y/%m/%d"))
    if (class(d) == "try-error" || is.na(d)) {
      print("Date format should be YYYY-mm-dd")
    }
  }
  if (!is.null(added.up.to)) {
    d <- try(as.Date(added.up.to, format = "%Y/%m/%d"))
    if (class(d) == "try-error" || is.na(d)) {
      print("Date format should be YYYY-mm-dd")
    }
  }


  if(!is.null(tumor)){
    id <- sapply(tumor, function(x){
      grepl(x, db$Disease, ignore.case = TRUE)
    })
    id <- apply(id, 1,any)
    db <-  db[id,]
  }
  if(!is.null(platform)){
    id <- sapply(platform, function(x){
      grepl(x, db$Platform, ignore.case = TRUE)
    })
    id <- apply(id, 1,any)
    db <-  db[id,]
  }
  if(!is.null(center)){
    id <- sapply(center, function(x){
      grepl(x, db$Center, ignore.case = TRUE)
    })
    id <- apply(id, 1,any)
    db <-  db[id,]
  }
  if(!is.null(level)){
    id <- grep(paste0("Level_", level), db$name)
    if(length(id) > 0){
      db <-  db[id,]
    }
  }

  if (!is.null(added.since)) {
    db <- subset(db, as.Date(db$addedDate) > as.Date(added.since,
                                                     "%m/%d/%Y"))
  }
  if (!is.null(added.up.to)) {
    db <- subset(db, as.Date(db$addedDate) < as.Date(added.up.to,
                                                     "%m/%d/%Y"))
  }

  # to be improved
  idx <- c()
  if(!is.null(samples)){
    for(i in seq_along(samples)){
      aux <- grep(samples[i],db$barcode)
      idx <- union(idx, aux)
    }
    db <- db[idx,]
  }
  return(db)
}

