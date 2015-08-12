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
        extra <- paste0(extra, "[Platform[@name=", platform, "]]")
        pages <- pages / 4
    }
    if (!is.null(disease)) {
        extra <- paste0(extra, "[Disease[@abbreviation=", disease, "]]")
        pages <- pages / 4
    }
    if (!is.null(type)) {
        if (is.numeric(type)) {
            extra <- paste0(extra, "[ArchiveType[@type=", type, "]]")
        } else {
            extra <- paste0(extra, "[ArchiveType[@type=mage-tab]]")
        }
    }
    if (!is.null(center)) {
        extra <- paste0(extra, "[Center[@name=", center, "]]")
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
        db <- db[order(db$Disease,db$Platform,db$Center),]
    }
    return(db)
}

# Using the name create two collumns Platform and Disease
tcgaDbAddCol <- function(data) {
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
    db <- data.frame()

    next.url <- url
    regex <- "<table summary=\"Data Summary\".*</a></td></tr></table>"
    next.regex <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"

    if (max > 0) {
        pb <- txtProgressBar(min = 0, max = max, style = 3)
        i <- 0
    }
    while (!is.na(next.url)) {
        # As we have a limitation of 1000 connections
        # every 3 minutes, we need to verify it
        download(next.url,
                 "tcga.html",
                 quiet = TRUE)

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

        db <- rbind(db, table)
        # get next table
        next.url <- str_match(html, next.regex)[idx, ]
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
    unlink("tcga.html")
    return(db)
}

#' @title TCGAUpdate
#' @description Updates the preprocessed TCGA database used by the package.
#'
#' This function can update the the disease, platform, center, data table.
#' Updating disease, platform, center will take some seconds. And Updating the
#' data.table shoul take no more than 10 minutes.
#'
#' The package will be updated with lastest version of the table every week.
#' @return platform/center/disease/data tables will be updated in the package
#' @importFrom devtools use_data load_all
#' @keywords internal
TCGAUpdate <- function(){
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

    # get new version of files
    new.db <-  createTcgaTable()
    # copy not modified ones
    for (i in seq_along(new.db[,1])){
        db <- subset(tcga.db,new.db[i,"name"] == tcga.db$name)
        if (nrow(db) == 1){
            new.db[i,"deployStatus"] <- db$barcode
        }
    }

    idx <- ((new.db$deployStatus == "" |  new.db$deployStatus == "Not found" |
                 (new.db$deployStatus == "Available")) &
                !grepl("aux|mage-tab", new.db$name)
    )
    new.db[idx,]$deployStatus <- "Available"
    new.db[idx,]$deployStatus <- getBarcode(new.db[idx,])$barcode
    colnames(new.db)[4] <- "barcode"
    tcga.db <- new.db

    gene.location <- get.GRCh.bioMart()

    use_data(platform.table, disease.table, tcga.db, center.table,
             DAVID_BP_matrix,DAVID_CC_matrix,DAVID_MF_matrix,
             EAGenes,gene.location,listEA_pathways,
             lgg.subtype, gbm.subtype, luad.subtype,
             stad.subtype, brca.subtype, coad.subtype,
             internal = TRUE,overwrite = TRUE)
}


# Get latest Genome Reference Consortium Human Build And save
# it as Genomic Ranges
#' @importFrom biomaRt useMart getBM
#' @keywords internal
get.GRCh.bioMart <- function(genome="hg19") {

    if (genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "grch37.ensembl.org",
                           path = "/biomart/martservice" ,
                           dataset = "hsapiens_gene_ensembl")
    } else {
        # for hg38
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    }

    chrom <- c(1:22, "X", "Y")
    gene.location <- getBM(attributes = c("chromosome_name",
                                          "start_position",
                                          "end_position", "strand",
                                          "external_gene_name",
                                          "entrezgene"),
                           filters = c("chromosome_name"),
                           values = list(chrom), mart = ensembl)

    return(gene.location)
}
