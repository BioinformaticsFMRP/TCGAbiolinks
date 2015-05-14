# @description get the arhive info from tcga api
getArchive <- function(id) {
  root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  query <- paste0("query=Archive&FileInfo[@id=", id,
                  "]&roleName=archiveCollection")
  url <- paste0(root, query)
  db <- tcgaGetTable(url)
  return(db)
}

getSamplesFiles <- function(id) {
  archives <- NULL
  root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  query <- paste0("query=FileInfo&BiospecimenBarcode[@id=",
                  id, "]&roleName=fileCollection")
  url <- paste0(root, query)
  db <- tcgaGetTable(url)
  if (!is.null(db)) {
    obsolete <- grep("-", db$md5sum)

    if (length(obsolete) > 0) {
      db <- db[-obsolete, ]
    }
    files <- unique(db$name)
    # for each file, get the highest ID
    for (i in seq_along(files)) {
      same.files <- grep(files[i], db$name)
      latest <- max(as.integer(db[same.files, ]$id))
      if (!exists("archives")) {
        archives <- getArchive(latest)
        archives$file <- files[i]

      } else {
        aux <- getArchive(latest)
        aux$file <- files[i]
        archives <- rbind(archives, aux)
      }
    }
    archives <- archives[, -c(10:15)]
    archives$addedDate <- as.Date(archives$addedDate, "%m-%d-%Y")
    archives <- tcgaDbAddCol(archives)
  }
  return(archives)
}

getBcrArchiveCollection <- function(id) {
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=Archive&BiospecimenBarcode[@id=",
                       id, "]&roleName=bcrArchiveCollection")
  url <- paste0(tcga.root, tcga.query)
  db <- tcgaGetTable(url)
  db$addedDate <- as.Date(db$addedDate, "%m-%d-%Y")
  db <- tcgaDbAddCol(db)
}
getBarcodeTable <- function(barcode) {
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=BiospecimenBarcode[@barcode=",
                       barcode, "][@isValid=true]")
  url <- paste0(tcga.root, tcga.query)
  db <- tcgaGetTable(url)
  return(db)
}
# This function created tcga.db warning to create the db
# takes almost 3 days to update it should take some minutes
# (function still needed to be done)
#' @importFrom downloader download
createTcgaTable <- function(disease = NULL, platform = NULL,
                            type = NULL, center = NULL) {
  # get all compressed archives
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  message("Downloading TCGA database")
  tcga.query <- paste0("query=Archive[isLatest=1]")

  extra <- ""
  pages <- 32
  if (!is.null(platform)) {
    extra <- paste0(extra, "[Platform[@name=", platform,
                    "]]")
    pages <- pages / 4
  }
  if (!is.null(disease)) {
    extra <- paste0(extra, "[Disease[@abbreviation=", disease,
                    "]]")
    pages <- pages / 4
  }
  if (!is.null(type)) {
    if (is.numeric(type)) {
      extra <- paste0(extra, "[ArchiveType[@type=", type, "]]")
    } else {
      extra <- paste0(extra, "[ArchiveType[@type=mage-tab]]")
    }
    if (!is.null(center)) {
      extra <- paste0(extra, "[Center[@name=", center, "]]")
    }

    pages <- pages / 4
  }

  url <- paste0(tcga.root, tcga.query, extra)
  print(url)
  db <- tcgaGetTable(url, pages)
  # remove useless cols
  db <- db[, 1:9]

  # remove protected data from the update
  idx <- grep("tcga4yeo", db$deployLocation)
  if (length(idx) > 0) {
    db <- db[-idx, ]
  }
  if (!is.null(db)) {
    # transoform date into same format
    db$addedDate <- as.Date(db$addedDate, "%m-%d-%Y")
    db <- tcgaDbAddCol(db)
  }
  return(db)
}

# Using the name create two collumns Platform and Disease
tcgaDbAddCol <- function(data) {
  disease.table  <- get("disease.table")
  platform.table  <- get("platform.table")
  center.table  <- get("center.table")
  data$Platform <- ""
  data$Disease <- ""
  data$Center <- ""
  diseases <- disease.table$abbreviation
  for (i in seq_along(diseases)) {
    idx <- grep(diseases[i], data$baseName)
    if (length(idx > 0)) {
      data[idx, ]$Disease <- diseases[i]
    }
  }

  for (i in seq_along(platform.table$name)) {
    idx <- grepl(platform.table[i, ]$name, data$baseName)
    idx <- idx & (data$Platform == "")
    if (any(idx)) {
      data[idx, ]$Platform <- platform.table[i, ]$name
    }
  }

  for (i in seq_along(center.table$name)) {
    idx <- grepl(center.table[i, ]$name, data$baseName)
    idx <- idx & (data$Center == "")
    if (any(idx)) {
      data[idx, ]$Center <- center.table[i, ]$name
    }
  }

  return(data)
}

# A function to get a tcga table from api input: url (max for
# progress bar) return: table
#' @importFrom stringr str_match
#' @importFrom downloader download
tcgaGetTable <- function(url, max = 0) {
  db <- NULL
  next.url <- url
  regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
  next.regex <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"

  if (max > 0) {
    pb <- txtProgressBar(min = 0, max = max, style = 3)
    i <- 0
  }
  # print(url)
  while (!is.na(next.url)) {
    download(next.url, "tcga.html", quiet = TRUE)
    html <- readLines("tcga.html")
    match <- str_match(html, regex)
    idx <- which(!is.na(match))
    if (length(idx) == 0) {
      next.url <- NA
      break
    }
    table <- readHTMLTable(toString(match[idx, ]), header = TRUE,
                           stringsAsFactors = FALSE)$"NULL"
    colnames(table) <- table[1, ]
    table <- table[-1, ]

    if (exists("db")) {
      db <- rbind(db, table)
    } else {
      db <- table
    }
    # get next table
    next.url <- str_match(html, next.regex)[6, ]
    next.url <- gsub("amp;", "", gsub("\">Next", "", next.url))

    if (max > 0) {
      # update progress bar
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
  }
  if (max > 0) {
    setTxtProgressBar(pb, max)
    close(pb)
  }
  return(db)
}


tcgaUpdate <- function(){

  tcga.db <-  get("tcga.db")

  # get new version of files
  new.db <-  createTcgaTable()
  # copy not modified ones
  for (i in seq_along(new.db[,1])){
    db <- subset(tcga.db,new.db[i,"name"] == tcga.db$name)
    new.db[i,"deployStatus"] <- db$barcode
  }

  new.db <- getBarcode(new.db)
  return(new.db)
  #tcga.db <- new.db
  #save(platform.table, disease.table, tcga.db, center.table,
  #     file = paste0(system.file("extdata", package = "TCGAbiolinks"),
  #                   "/dataFolders.rda")
  #     )
}
