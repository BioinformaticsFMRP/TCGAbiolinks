#' @title TCGA TCGAPrepare
#' @description
#'  Prepare data matrices for downstream analysis,
#'  ready to use also with other R packages
#'  DNA Methylation
#'   Row: matrix with genes/loci
#'   Cols: samples in columns
#' @return Data frame with data read
#' @examples
#'  TCGAPrepare2(query,dir="path")
#' @export
TCGAPrepare <- function(query, dir = NULL, type = NULL){

    if (is.null(dir)) {
        message("Argument dir is NULL. Plese provide the directory
            with the folders to be prepared. ")
        return(NULL)
    }
    if (length(unique(query$Platform)) > 1 |
        length(unique(query$Center)) > 2) {
        message("Sorry! But, for the moment, we can only prepare on type of
            platform per call")
        return(NULL)
    } else {
        platform <- unique(query$Platform)
    }

    files <- NULL
    dirs <- gsub(".tar.gz","",basename(query$deployLocation))
    for (i in seq_along(dirs)) {
        aux <- list.files(file.path(dir,dirs[i]), full.names = TRUE,
                          recursive = TRUE)
        files <- c(files, aux )
    }
    print(files)
    idx <- grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE",files)
    if(length(idx) > 0){
        files <- files[-idx]
    }
    print(files)
    if(!is.null(type)){
        files <- files[grep(type,files)]
    }

    df <- NULL
    if (grepl("humanmethylation",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE)
            sample <- gsub("\\.", "-", colnames(data)[2])
            colnames(data) <- data[1,]
            data <- data[-1,] # removing Composite Element REF
            colnames(data)[2] <- sample
            if (i == 1) {
                df <- data[, c(1, 3:5, 2)]
            } else {
                df <- merge(df, data[, c(1, 2)],
                            by = "Composite Element REF")
            }
        }
        rownames(df) <- df$Composite.Element.REF
        message("Removing X Y chromossomes")
        df <- df[df$Chromosome != "X" & df$Chromosome != "Y", ]
        # methylation$Chromosome <- NULL

        # remove NA lines
        message("Removing NA Lines")
        df <- na.omit(df)

    }

    if (grepl("mda_rppa_core",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
            sample <- gsub("\\.", "-", colnames(data)[2])
            colnames(data) <- data[1,]
            data <- data[-1,] # removing Composite Element REF
            colnames(data)[2] <- sample
            print(files[i])

            print(sample)
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = "Composite Element REF")
            }
        }

        # get array_design.txt from mage folder
        # and change uuid by Barcode
        uuid <- colnames(df)
        idx <- grep("Sample|Control",uuid)
        uuid <- uuid[-idx]
        map <- mapuuidbarcode(uuid)
        idx <- which(colnames(df) %in% map$uuid)
        colnames(df)[idx] <- map$barcode

    }
    # case: header has barcode
    # Line 2 is useless
    if (grepl("agilent",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
            data <- data[-1,] # removing Composite Element REF
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = "Hybridization REF")
            }
        }
    }

    if (grepl("illuminaga",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE,
                               comment.char = "#",fill = TRUE)
            data <- data[-1,] # removing Composite Element REF
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = "Hybridization REF")
            }
        }
    }

    if(grepl("illuminahiseq_rnaseqv2|illuminahiseq_totalrnaseqv2",
             tolower(platform))){
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE,
                               comment.char = "#",fill = TRUE)
            data <- data[-1,] # removing Composite Element REF
            uuid <- unlist(strsplit(files[i],"\\."))[9]
            map <- mapuuidbarcode(uuid)
            colnames(data)[2] <- map$barcode
            if (i == 1) {
                df <- data[,c(1,2)]
            } else {
                df <- merge(df, data[,c(1,2)],by = "gene_id")
            }

        }
    }
    return(df)
}


#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapuuidbarcode <- function(uuids){
    # Using tcga api: https://goo.gl/M1uQLR
    serv <- "https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch"
    header <- c('Content-Type' = 'text/plain')
    uuids <- paste0(uuids, collapse = ",")
    ans <- fromJSON(postForm(serv,
                             .opts = list(postfields = uuids,
                                          httpheader = header,
                                          ssl.verifypeer = FALSE)))

    # transform to dataframe
    x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame))[[1]])
    barcodes <- seq(1,length(x),2)
    uuid <- seq(2,length(x),2)
    x <- data.frame(x[uuid],as.character(x[barcodes]),
                    stringsAsFactors = FALSE)
    colnames(x) <- c("uuid","barcode")
    return(x)
}

#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom RJSONIO fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapbarcodeuuid <- function(barcode){
    # Using tcga api: https://goo.gl/M1uQLR
    barcode <- paste0(barcode, collapse = ",")
    serv <- paste0("https://tcga-data.nci.nih.gov/",
                   "uuid/uuidws/mapping/json/barcode/batch")
    header <- c('Content-Type' = 'text/plain')
    ans <- fromJSON(postForm(serv,
                             .opts = list(postfields = barcode,
                                          httpheader = header,
                                          ssl.verifypeer = FALSE)))

    # transform to dataframe
    x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame))[[1]])
    barcodes <- seq(1,length(x),2)
    uuid <- seq(2,length(x),2)
    x <- data.frame(x[uuid],as.character(x[barcodes]),
                    stringsAsFactors = FALSE)
    colnames(x) <- c("uuid","barcode")
    return(x)
}

# Get sdrf file/array_design of a line
# example
# query <- TCGAQuery(tumor = "BRCA")
# getMagecontent(query[1,])
# Obs: delete the file after reading
#      is it better to save it?
getMage <- function(line){
    tcga.db <- get("tcga.db")
    root <- "https://tcga-data.nci.nih.gov"
    path <- "mages"
    dir.create(path,showWarnings = FALSE)
    mages <-  tcga.db[grep("mage-tab",tcga.db$name),]
    # get the mage file for the entry
    mage <- subset(mages, mages$Disease == line$Disease &
                       mages$Platform == line$Platform &
                       mages$Center == line$Center)

    if (dim(mage)[1] != 0) {
        file <- file.path(path,basename(mage$deployLocation))
        if ( !file.exists(file)) {
            download(paste0(root,mage$deployLocation), file, quiet = TRUE)
        }
        folder <- gsub(".tar.gz","",file)
        if ( !file.exists(folder)) {
            untar(file,exdir = "mages")
        }
        files <- list.files(folder)
        if (line$Platform == "MDA_RPPA_Core") {
            sdrf <- files[grep("array_design",files)]
        } else {
            # Platform is not MDA_RPPA_Core
            sdrf <- files[grep("sdrf",files)]
        }
        if (length(sdrf) > 1) {
            sdrf <- sdrf[1]
        }

        df <- read.delim(file = file.path(folder,sdrf),
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         fileEncoding = "latin1")
    }
    return(df)
}

