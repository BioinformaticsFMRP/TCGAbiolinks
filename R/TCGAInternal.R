# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach <- function (libname, pkgname){

    file = system.file("extdata/GRCh.rda",package = "TCGAbiolinks")
    load(file,envir = as.environment("package:TCGAbiolinks"))

    file = system.file("extdata/wine.rda",package = "TCGAbiolinks")
    load(file,envir = as.environment("package:TCGAbiolinks"))


    file = system.file("extdata/dataGeneExpression.rda",
                       package = "TCGAbiolinks")
    load(file,envir = as.environment("package:TCGAbiolinks"))

    file = system.file("extdata/dataEnrichmentAnalysis.rda",
                       package = "TCGAbiolinks")
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
        " Query, download & analyze - TCGA                  \n",
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
    platform.table <- tcgaGetTable(next.url)
    platform.table <- platform.table[, 1:4]
    platform.table <- platform.table[order(platform.table$name,
                                           decreasing = TRUE),]

    # Get disease table
    tcga.query <- "query=Disease"
    next.url <- paste0(tcga.root, tcga.query)
    disease.table <- tcgaGetTable(next.url)
    disease.table <- disease.table[, 1:3]

    # Get center table
    tcga.query <- "query=Center"
    next.url <- paste0(tcga.root, tcga.query)
    center.table  <- tcgaGetTable(next.url)
    center.table <- center.table[, 1:3]

    env$platform.table <- platform.table
    env$center.table <- center.table
    env$disease.table <- disease.table

    # Get tcga folders with barcodes without private folders
    tcga.db <- createTcgaTable()
    env$tcga.db <- tcga.db
    save(tcga.db, paste0(system.file("extdata", package = "TCGAbiolinks"),
                         "/tcgadb.rda"))
    step <- 200
    for(i in seq(1, nrow(tcga.db) , by = step)){
        print(i)
        j <- i + step
        if(j > nrow(tcga.db)){
            j <-nrow(tcga.db)
        }
        tcga.db[i:j,]$deployStatus <- getBarcode(tcga.db[i:j,])$barcode
        save(tcga.db, paste0(system.file("extdata", package = "TCGAbiolinks"),
                             "/tcgadb.rda"))
    }

    env$tcga.db <- tcga.db
    save(platform.table, disease.table, tcga.db, center.table,
         file = paste0(system.file("extdata", package = "TCGAbiolinks"),
                       "/dataFolders.rda")
    )
}

# Get all files in the ftp directory @keywords internal
getFileNames <- function(url) {

    download(url,
             "temp.html",
             mode = "wb",
             quiet = 1)
    x <- capture.output(XML::htmlTreeParse("temp.html"))
    unlink("temp.html")
    x <- x[grep("href", x)]
    if (is.null(x)){
        return(NULL)
    }
    x = sapply(strsplit(x, ">"), function(y) y[2])
    if (is.null(x)){
        return(NULL)
    }
    x = sapply(strsplit(x, "<"), function(y) y[1])
    return(x)
}

tcga.get.barcode <- function(data){
    # todo considere next page =/
    # comparar com sdrf
    # for each tcga.db id get barcodes
    message("Downloading TCGA barcodes")
    all.barcode <- c()

    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query="
    tcga.query <- paste0("FileInfo&Archive[@id=",data$id,
                         "]&roleName=fileCollection")

    next.url <- paste0(tcga.root,tcga.query)
    files <- tcgaGetTable(next.url)
    #print(files)
    if (nrow(files) == 0) {
        return(NULL)
    }
    files <- files[,1:4]
    idx <- grep("somatic.maf",files$name)
    if (length(idx > 0)) {
        files <- files[idx,]
    }
    else {
        # no maf file
        # try in the name
        pat <- paste0("((((TCGA-[A-Z0-9]{2})-[A-Z0-9]{4})-[A-Z0-9]{3}-",
                      "[A-Z0-9]{3})-[A-Z0-9]{4}-[A-Z0-9]{2})")
        barcode <- str_match(files$name,pat)[,1]
        #message("Found in name")
        if (!all(is.na(barcode))) {
            message("Found in name")
            barcode <- barcode[!is.na(barcode)]
            barcode <- paste0(unique(barcode), collapse = ",")
            return(barcode)
        }

        return(NULL)
    }

    # maybe two maf files
    for (i in  seq_along(files$name)) {
        tcga.query <- paste0("BiospecimenBarcode&FileInfo[@id=",files[i,"id"],
                             "]&roleName=biospecimenBarcodeCollection")
        next.url <- paste0(tcga.root,tcga.query)
        print(next.url)
        barcode.table <- tcgaGetTable(next.url)
        barcode.table <- barcode.table[,1:8]
        all.barcode <- union(all.barcode, unique(barcode.table$barcode))
    }

    all.barcode <- paste0(unique(all.barcode), collapse = ",")

    return(all.barcode)
}

#' @import utils
#' @importFrom RCurl getURL
.DownloadURL <-
    function(Site){
        # setInternet2(use = TRUE)
        Site <- URLencode(Site)
        x=  getURL(Site, ssl.verifypeer = FALSE)
        x <- unlist(strsplit(x,"\n"))
        return(x)
    }


is.windows <- function() {
    Sys.info()["sysname"] == "Windows"
}

is.mac <- function() {
    Sys.info()["sysname"] == "Darwin"
}

is.linux <- function() {
    Sys.info()["sysname"] == "Linux"
}

