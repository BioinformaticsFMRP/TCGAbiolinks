
tcgaQueryApi <- function(tumor = NULL, platform = NULL, added.since = NULL,
                         added.up.to = NULL, samples = NULL, level = NULL) {
  disease.table   <- get("disease.table")
  platform.table  <- get("platform.table")
  dna.plat        <- get("dna.plat")
  rna.plat        <- get("rna.plat")
  total.rna.plat  <- get("total.rna.plat")
  mirna.plat      <- get("mirna.plat")
  wgarubcon.plat  <- get("wgarubcon.plat")
  wgaqiagen1.plat <- get("wgaqiagen1.plat")
  wgaqiagen2.plat <- get("wgaqiagen2.plat")
  plat.center     <- get("plat.center")

  if (!is.null(tumor)) {
    if (!(is.element(tolower(tumor),
                     tolower(disease.table$abbreviation)))) {
      df <- as.data.frame(matrix(sort(unique(disease.table$abbreviation)),
                                 ncol = 8))
      print(kable(df, col.names = NULL, format = "pandoc",
                  caption = "TCGA tumors"))
      cat("=======================================================\n")
      cat("ERROR: Disease not found. Select from the table above.\n")
      cat("=======================================================\n")
      return(NULL)
    }
  }
  if (!is.null(platform)) {
    if (!(is.element(tolower(platform), tolower(platform.table$name)))) {
      df <- as.data.frame(matrix(sort(unique(platform.table$name)),
                                 ncol = 3))
      print(kable(df, col.names = NULL, format = "pandoc",
                  caption = "TCGA Platforms"))
      cat("=======================================================\n")
      cat("ERROR: Platform not found. Select from the table above.\n")
      cat("=======================================================\n")
      return(NULL)
    }
  }

  if (!is.null(level)) {
    if (!(is.element(level, c("1", "2", "3")))) {
      message("Level not found. Chosse between:'1', '2' or '3'")
      stop("Invalid platform")
    }
  }

  if (!is.null(added.since)) {
    d <- try(as.Date(added.since, format = "%m/%d/%Y"))
    if (class(d) == "try-error" || is.na(d)) {
      print("Date format should be mm/dd/YYYY")
    }
  }
  if (!is.null(added.up.to)) {
    d <- try(as.Date(added.up.to, format = "%m/%d/%Y"))
    if (class(d) == "try-error" || is.na(d)) {
      print("Date format should be mm/dd/YYYY")
    }
  }

  # to be improved
  if (!is.null(samples)) {
    # archives <- c()
    files <- c()
    for (i in seq_along(samples)) {
      # table with barcode id example:
      # query=BiospecimenBarcode[@barcode=TCGA-28-2499*]
      message("Searching for barcode files...")
      message(paste("Barcode:", samples[i]))
      db <- getBarcodeTable(samples[i])

      # Improvement: using portion analyte in order to select
      # platform
      if (!is.null(platform)) {
        analyte <- c()
        pat <- paste0("TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-",
                      "[[:alnum:]]{3}-[[:alnum:]]{2}")
        if (is.element(tolower(platform), tolower(dna.plat))) {
          analyte <- c(analyte, "D")
        }
        if (is.element(tolower(platform), tolower(rna.plat))) {
          analyte <- c(analyte, "R")
        }
        if (is.element(tolower(platform), tolower(total.rna.plat))) {
          analyte <- c(analyte, "T")
        }
        if (is.element(tolower(platform), tolower(wgarubcon.plat))) {
          analyte <- c(analyte, "G")
        }
        if (is.element(tolower(platform), tolower(mirna.plat))) {
          analyte <- c(analyte, "H")
        }
        if (is.element(tolower(platform), tolower(wgaqiagen2.plat))) {
          analyte <- c(analyte, "X")
        }
        if (is.element(tolower(platform), tolower(wgaqiagen1.plat))) {
          analyte <- c(analyte, "W")
        }
        idx <- c()
        for (i in seq_along(analyte)) {
          aux <- grep(paste0(pat, "[", analyte[i], "]"),
                      as.character(db$barcode))
          idx <- union(idx, aux)
        }
        print(length(idx))
        db <- db[idx, ]
      }
      if (!is.null(platform)) {
        pat <- paste0("TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-",
                      "[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-")
        with(plat.center,{
          idx <- grep(platform, plat.center$platform)
          centers <- as.character(plat.center[idx, "center"])
        })
        idx <- c()
        for (i in seq_along(centers)) {
          aux <- grep(paste0(pat, centers[i]),
                      as.character(db$barcode))
          idx <- union(idx, aux)
        }
        print(length(idx))
        db <- db[idx, ]
      }

      # get getBcrArchiveCollection table
      for (i in seq_along(db$id)) {
        # aux <- getBcrArchiveCollection(db[i,'id']) archives <-
        # rbind(archives,aux)
        aux <- getSamplesFiles(db[i, "id"])
        aux$barcode <- db[i, "barcode"]
        if (!is.null(aux)) {
          files <- rbind(files, aux)
        }
      }
    }
    x <- subset(files, files$isLatest == 1)
    x <- x[!duplicated(x), ]
    x <- x[, order(names(x))]
    if (!is.null(platform)) {
      x <- subset(x, grepl(tolower(x$platform), tolower(x$Platform)))
    }
    if (!is.null(platform)) {
      x <- subset(x, tolower(x$Disease) == tolower(tumor))
    }
    if (!is.null(level)) {
      x <- subset(x, grepl(paste0("Level_", level), x$name))
    }

  } else {
    message("CREATING TABLE")
    x <- createTcgaTable(platform = platform, type = level,
                         disease = tumor)
  }
  if (!is.null(added.since)) {
    x <- subset(x, as.Date(x$addedDate) > as.Date(added.since,
                                                  "%m/%d/%Y"))
  }
  if (!is.null(added.up.to)) {
    x <- subset(x, as.Date(x$addedDate) < as.Date(added.up.to,
                                                  "%m/%d/%Y"))
  }

  return(x)
}


getBarcode <- function(table){
  tcga.db <- get("tcga.db")
  root <- "https://tcga-data.nci.nih.gov"
  allFiles <- c()
  mages <-  tcga.db[grep("mage-tab",tcga.db$name),]
  #message(dim(mages))
  for (i in seq_along(table[,1])){
    #message(i)
    if (table[i,]$deployStatus == "Available") {
      #message(i,table[i,]$Platform)
      #mage <- createTcgaTable(disease = table[i,]$Disease,
      #                        platform = table[i,]$Platform,
      #                        center =  table[i,]$Center,
      #                        type = "mage-tab")

      mage <- subset(mages, mages$Disease == table[i,]$Disease &
                       mages$Platform == table[i,]$Platform &
                       mages$Center == table[i,]$Center)
      if (dim(mage)[1] == 0){
        table[i,]$deployStatus <- "Not found"
        next
      }
      file <- basename(mage$deployLocation)
      #print(file)
      #print(table[i,]$baseName)
      x <- stringr::str_replace_all(file, "[^[:alnum:]]", "")
      y <- stringr::str_replace_all(table[i,]$baseName, "[^[:alnum:]]", "")
      idx <- grep(y,x)

      if (length(idx) > 0) {
        file <- file[idx]
        mage <- mage[idx,]
      }
      print(file)
      allFiles <- c(allFiles, file)
      if ( !file.exists(file)) {
        download(paste0(root,mage$deployLocation), file, quiet = TRUE)
        untar(file)
      }
      folder <- gsub(".tar.gz","",file)
      files <- list.files(folder)
      print(table[i,]$Platform)
      if (table[i,]$Platform == "MDA_RPPA_Core") {
        sdrf <- files[grep("array_design",files)]
        # case with 2 array_design BLCA
        if(length(sdrf) > 1){
          sdrf <- sdrf[1]
        }
        df <- read.delim(file = file.path(folder,sdrf), sep = "\t",
                          stringsAsFactors = FALSE, fileEncoding="latin1")
        if(is.element("Sample.description",colnames(df))){
        barcode <- unique(as.character(df$Sample.description))
        } else {
          barcode <- unique(as.character(df$Biospecimen.Barcode))
        }
        table[i,]$deployStatus <- paste0(barcode, collapse = ",")
      } else {
        sdrf <- files[grep("sdrf",files)]
        df <- read.delim(file = file.path(folder,sdrf) ,sep = "\t",
                         stringsAsFactors = FALSE, fileEncoding="latin1")
        for (j in seq_along(table[,1])) {
          if (table[j,]$Disease == table[i,]$Disease &&
                table[j,]$Platform == table[i,]$Platform &&
                table[j,]$Center == table[i,]$Center) {
            aux <- grep("Comment..TCGA.Archive.Name",colnames(df))
            #len <- length(grep("Comment..TCGA.Archive.Name",colnames(df)))
            #print(len)
            #print(file.path(folder,sdrf))

            barcode <- data.frame()
            x <- table[j,]$name

            for (z in seq_along(aux)){
              barcode <- rbind(barcode,subset(df, x == df[,aux[z]]))
            }

            #           if(len == 4){
            #             x <- table[j,]$name
            #             barcode <- subset(df,
            #                               x == df$Comment..TCGA.Archive.Name..3 |
            #                                 x == df$Comment..TCGA.Archive.Name..2 |
            #                                 x == df$Comment..TCGA.Archive.Name. |
            #                                 x == df$Comment..TCGA.Archive.Name..1
            #             )
            #           }
            #
            #           if(len == 3){
            #             x <- table[j,]$name
            #             barcode <- subset(df,
            #                               x == df$Comment..TCGA.Archive.Name..2 |
            #                                 x == df$Comment..TCGA.Archive.Name. |
            #                                 x == df$Comment..TCGA.Archive.Name..1
            #             )
            #           }
            #           if(len == 2){
            #             x <- table[j,]$name
            #             barcode <- subset(df, x == df$Comment..TCGA.Archive.Name. |
            #                                 x == df$Comment..TCGA.Archive.Name..1
            #             )
            #           }
            #           if(len == 1){
            #             x <- table[j,]$name
            #             barcode <- subset(df, x == df$Comment..TCGA.Archive.Name.)
            #           }
            barcode <- unique(as.character(barcode$Comment..TCGA.Barcode.))
            barcode <- barcode[barcode != "->"]
            table[j,]$deployStatus <- paste0(barcode,collapse = ",")
          }
        }
      }
    }
  }
  # removing the mess
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

  db <- tcga.db
  if(!is.null(tumor)){
    id <- sapply(tumor, function(x){db$Disease == x} )
    id <- apply(id, 1,any)
    db <-  db[id,]
  }
  if(!is.null(platform)){
    id <- sapply(platform, function(x){db$Platform == x})
    id <- apply(id, 1,any)
    db <-  db[id,]
  }
  if(!is.null(center)){
    id <- sapply(center, function(x){db$Center == x})
    id <- apply(id, 1,any)
    db <-  db[id,]
  }
  if(!is.null(level)){
    id <- grep(paste0("Level_", level), db$name)
    if(length(idx) > 0){
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

