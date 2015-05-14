# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach <- function (libname, pkgname){

  file = system.file("extdata/GRCh.rda",package = "TCGAbiolinks")
  load(file,envir = as.environment("package:TCGAbiolinks"))

  file = system.file("extdata/dataFolders.rda",package = "TCGAbiolinks")
  time <- file.info(file)$ctime
  if (file.exists(file)) {
    load(file, envir = as.environment("package:TCGAbiolinks"))
  } else {
    env <- as.environment("package:TCGAbiolinks")
    load.tcga(env)
  }

  if (!interactive() || stats::runif(1) > 0.1) return()
  welcome.message <- paste0(
    " =============================================================\n",
    " ______  ___  ____   ___                                        \n",
    "   ||   |    |      |   | |    o  __  |   o  _         __         \n",
    "   ||   |    | ___  |___| |__  | |  | |   | | | | |_/ |__         \n",
    "   ||   |___ |____| |   | |__| | |__| |__ | | |_| | \\  __|       \n",
    " ------------------------------------------------------------\n",
    " Search, download & analyse - TCGA                  \n",
    " Version:",utils::packageVersion("TCGAbiolinks"),"\n",
    " Last TCGAUpdate(): ",time,"\n",
    " ==============================================================\n"
  )
  packageStartupMessage(welcome.message)

}

# Updates tcga platform and diseases
# @param env package environment
#' @importFrom stringr str_match
#' @importFrom XML readHTMLTable
#' @importFrom downloader download
#' @keywords internal
load.tcga <- function(env) {
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"

  # Get platform table
  tcga.query <- "query=Platform"
  next.url <- paste0(tcga.root, tcga.query)
  download(next.url, "tcga.html", quiet = TRUE)
  regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
  html <- readLines("tcga.html")
  platform.table <- readHTMLTable(toString(str_match(html,regex)[6, ]),
                                  header = TRUE,
                                  stringsAsFactors = FALSE)$"NULL"
  colnames(platform.table) <- platform.table[1, ]
  platform.table <- platform.table[-1, 1:4]
  platform.table <- platform.table[order(platform.table$name,
                                         decreasing = TRUE),]
  env$platform.table <- platform.table

  # Get disease table
  tcga.query <- "query=Disease"
  next.url <- paste0(tcga.root, tcga.query)
  download(next.url, "tcga.html", quiet = TRUE)
  regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
  html <- readLines("tcga.html")
  match <- str_match(html, regex)
  idx <- which(!is.na(match))
  if (length(idx) > 0) {
    disease.table <- readHTMLTable(toString(match[idx, ]),
                                   header = TRUE,
                                   stringsAsFactors = FALSE)$"NULL"
    colnames(disease.table) <- disease.table[1, ]
    disease.table <- disease.table[-1, 1:3]
    env$disease.table <- disease.table
  }

  # Get center table
  tcga.query <- "query=Center"
  next.url <- paste0(tcga.root, tcga.query)
  download(next.url, "tcga.html", quiet = TRUE)
  html <- readLines("tcga.html")
  match <- str_match(html, regex)
  idx <- which(!is.na(match))
  if (length(idx) > 0) {
    center.table <- readHTMLTable(toString(match[idx, ]),
                                  header = TRUE,
                                  stringsAsFactors = FALSE)$"NULL"
    colnames(center.table) <- center.table[1, ]
    center.table <- center.table[-1, 1:3]
    env$center.table <- center.table
  }

  if (file.exists("tcga.html")) {
    file.remove("tcga.html")
  }
  # Get tcga folders with barcodes without private folders
  tcga.db <-  createTcgaTable()
  env$tcga.db <- tcga.db
  tcga.db <- getBarcode(tcga.db)

  env$tcga.db <- tcga.db
  save(platform.table, disease.table, tcga.db, center.table,
       file = paste0(system.file("extdata", package = "TCGAbiolinks"),
                     "/dataFolders.rda")
  )
}

# Get all files in the ftp directory @keywords internal
getFileNames <- function(url) {

  if (RCurl::url.exists(url)) {
    download(url,
             "temp.html",
             mode = "wb",
             quiet = 1)
    x <- capture.output(XML::htmlTreeParse("temp.html"))
    unlink("temp.html")
  } else {
    stop("Can't find URL. Please check the web site or the internet connection.")
  }

  x <- x[grep("href", x)]
  x = sapply(strsplit(x, ">"), function(y) y[2])
  x = sapply(strsplit(x, "<"), function(y) y[1])
  return(x)
}
