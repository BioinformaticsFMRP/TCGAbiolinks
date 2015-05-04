#' @title .onAttach
#' @description  Load required data into gloval enviroment
#' @keywords internal
.onAttach <- function (libname, pkgname){

  file = system.file("extdata/dataFolders.rda",package="TCGAbiolinks")
  time <- file.info(file)$ctime
  if(file.exists(file)){
    load(file,envir = as.environment("package:TCGAbiolinks"))
  } else {
    env <- as.environment("package:TCGAbiolinks")
    load.tcga(env)
  }
  cat("\014")
  welcome.message <- paste0(
    " =============================================================\n",
    " ______  ___  ____   ___                                        \n",
    "   ||   |    |      |   | |    o  __  |   o  _         __         \n",
    "   ||   |    | ___  |___| |__  | |  | |   | | | | |_/ |__         \n",
    "   ||   |___ |____| |   | |__| | |__| |__ | | |_| | \\  __|       \n",
    " ------------------------------------------------------------\n",
    " Search, download & analyse - TCGA                  \n",
    " Version:0.01 \n",
    " Last TCGAUpdate(): ",time,"\n",
    " ==============================================================\n"
  )
  packageStartupMessage(welcome.message)

}

#' @import XML stringr
load.tcga <- function(env){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- "query=Platform"
  next.url <- paste0(tcga.root,tcga.query)
  downloader::download(next.url,"tcga.html",quiet =T)
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  html <- readLines("tcga.html")
  platform.table <- readHTMLTable(toString(str_match(html,regex)[6,]),
                                  header = T,
                                  stringsAsFactors = FALSE)$'NULL'
  colnames(platform.table) <- platform.table[1,]
  env$platform.table <- platform.table[-1,1:4]

  tcga.query <- "query=Disease"
  next.url <- paste0(tcga.root,tcga.query)
  downloader::download(next.url,"tcga.html",quiet =T)
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  html <- readLines("tcga.html")
  disease.table <- readHTMLTable(toString(str_match(html,regex)[6,]),
                                  header = T,
                                  stringsAsFactors = FALSE)$'NULL'
  colnames(disease.table) <- disease.table[1,]
  env$disease.table <- disease.table[-1,1:4]

  if (file.exists('tcga.html')) {file.remove('tcga.html')}
  save(platform.table,disease.table,
       file = paste0(system.file("extdata", package="TCGAbiolinks"),"/dataFolders.rda")
  )
}
