#' @title TCGA TCGAPrepare
#' @description
#'  Prepare data matrices for downstream analysis,
#'  ready to use also with other R packages
#'  DNA Methylation
#'   Row: matrix with genes/loci
#'   Cols: samples in columns
#' @return Data frame with data read
#' @param query Data frame as the one returned from TCGAQuery
#' @param dir Directory with the files
#' @param type File to prepare.
#' @examples
#' sample <- "TCGA-06-0939-01A-01D-1228-05"
#' query <- TCGAQuery(tumor = "GBM",samples = sample, level = 3)
#' TCGADownload(query,path = "exampleData",samples = sample, quiet = TRUE)
#' prepared <- TCGAPrepare(query, dir="exampleData")
#' @export
#' @importFrom stringr str_match str_trim
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
    idx <- grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE",files)
    if (length(idx) > 0) {
        files <- files[-idx]
    }
    if (!is.null(type)) {
        files <- files[grep(type,files)]
        if(length(files) == 0){
            message("No files of that type found")
            return(NULL)
        }
    }


    df <- NULL
    if (grepl("humanmethylation",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE)
            sample <- gsub("\\.", "-", colnames(data)[2])
            colnames(data) <- gsub(" ", "\\.", data[1,])
            data <- data[-1,] # removing Composite Element REF
            colnames(data)[2] <- sample
            if (i == 1) {
                df <- data[, c(1, 3:5, 2)]
            } else {
                df <- merge(df, data[, c(1, 2)],
                            by = "Composite.Element.REF")
            }
        }
        rownames(df) <- df$Composite.Element.REF
        # remove NA lines
        df[,3:ncol(df)] <- sapply(df[,3:ncol(df)], as.numeric)
    }

    if (grepl("mda_rppa_core",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
            sample <- gsub("\\.", "-", colnames(data)[2])
            colnames(data) <- data[1,]
            data <- data[-1,] # removing Composite Element REF
            colnames(data)[2] <- sample

            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = "Composite Element REF")
            }
        }
        rownames(df) <- df[,1]
        df[,1] <- NULL
        # get array_design.txt from mage folder
        # and change uuid by Barcode
        uuid <- colnames(df)
        print(uuid)
        idx <- grep("Sample|Control",uuid)
        if(length(idx) > 0){
            uuid <- uuid[-idx]
        }
        map <- mapuuidbarcode(uuid)
        idx <- which(colnames(df) %in% map$uuid)
        colnames(df)[idx] <- as.character(map$barcode)
    }
    # case: header has barcode
    # Line 2 is useless
    if (grepl("agilent|H-miRNA_8x15K",platform, ignore.case = TRUE)) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], skip=2,
                               stringsAsFactors = FALSE)
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = colnames(df)[1])
            }
        }
        rownames(df) <- df[,1]
        df[,1] <- NULL
        colnames(df) <- sapply(files,
                               function(x) {
                                   read.table(x,
                                              stringsAsFactors = FALSE,
                                              nrows=1)[3]
                               })
    }

    if (grepl("illuminaga|illuminadnamethylation_oma",
              platform, ignore.case = TRUE)) {

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
        rownames(df) <- df[,1]
        df[,1] <- NULL
    }

    if (grepl("illuminahiseq_rnaseqv2|illuminahiseq_totalrnaseqv2",
              tolower(platform))) {
        regex <- paste0("[:alnum:]{8}-[:alnum:]{4}",
                        "-[:alnum:]{4}-[:alnum:]{4}-[:alnum:]{12}")
        uuid <- str_match(files,regex)
        map <- mapuuidbarcode(uuid)
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE,
                               comment.char = "#",fill = TRUE)
            data <- data[-1,] # removing Composite Element REF
            x <- subset(map, uuid ==uuid[i])
            colnames(data)[2] <- as.character(x$barcode)
            if (i == 1) {
                df <- data[,c(1,2)]
            } else {
                df <- merge(df, data[,c(1,2)],by = colnames(df)[1])
            }
        }
        rownames(df) <- df[,1]
        df[,1] <- NULL
    }

    if (grepl("illuminahiseq_mirnaseq",platform, ignore.case = TRUE)) {
        files <- files[grep("mirna",files)]
        if(length(files) == 0){
            message("No mirna files of that type found")
            return(NULL)
        }

        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- str_match(files,regex)

        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE,
                               comment.char = "#",fill = TRUE)
            if (i == 1) {
                df <- data[,"read_count"]
            } else {
                df <- cbind(df, data[,"read_count"])
            }
        }
        colnames(df) <- as.character(barcode)
        rownames(df) <- data[,1]
    }

    if (grepl("bio",platform,ignore.case = TRUE)) {
        if (!is.null(type)) {
            files <- files[grep(type,files)]
        }
        if (length(files) == 1) {
            df <- read.table(files, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE,
                             comment.char = "#",fill = TRUE)
            message("Please, remove useless info")
            # removing CDE line and duplicated header
        } else {
            message("We're preaparing for the moment only one clinical file")
            return(NULL)
        }
    }

    if (grepl("genome_wide_snp_6",tolower(platform))) {
        files <- files[-grep("nocnv", files)]

        #load("hg19genes.RData")
        genes <- sort(unique(hg19genes$gene_name))

        df <- matrix(0, nrow = length(genes), ncol = length(files))
        rownames(df) <- genes

        colNames <- rep("", ncol(df)) #check barcode
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE)
            ####Check barcodes
            colNames[i] <- data$Sample[1]
            infofiles <- paste("n. ",i," of ", length(files), " done..",sep="" )
            print(infofiles)
            for(j in 1:nrow(data)){
                gg <- sort(unique(hg19genes[hg19genes$start >= data$Start[j] & hg19genes$end <= data$End[j], "gene_name"]))
                df[gg, i] <- df[gg, i] + data$Segment_Mean[j]

                #print(j)
            }
        }
        colnames(df) <- colNames
    }

    return(df)
}


#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapuuidbarcode <- function(uuid){
    # Using tcga api: https://goo.gl/M1uQLR
    serv <- "https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch"
    header <- c('Content-Type' = 'text/plain')
    uuids <- paste0(uuid, collapse = ",")
    ans <- fromJSON(postForm(serv,
                             .opts = list(postfields = uuids,
                                          httpheader = header,
                                          ssl.verifypeer = FALSE)))
    if(length(uuid) == 1){
        x <- data.frame(ans$uuidMapping$barcode,ans$uuidMapping$uuid,
                        stringsAsFactors = FALSE)
        colnames(x) <- c("barcode","uuid")
    } else {
        # Extract patient barcode from sample barcode.
        x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame)))
    }
    return(x)
}

#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapbarcodeuuid <- function(barcode){
    # Using tcga api: https://goo.gl/M1uQLR
    barcodes <- paste0(barcode, collapse = ",")
    serv <- paste0("https://tcga-data.nci.nih.gov/",
                   "uuid/uuidws/mapping/json/barcode/batch")
    header <- c('Content-Type' = 'text/plain')
    ans <- fromJSON(postForm(serv,
                             .opts = list(postfields = barcodes,
                                          httpheader = header,
                                          ssl.verifypeer = FALSE)))

    if(length(barcode) == 1){
        x <- data.frame(ans$uuidMapping$barcode,ans$uuidMapping$uuid,
                        stringsAsFactors = FALSE)
        colnames(x) <- c("barcode","uuid")
    } else {
        # Extract patient barcode from sample barcode.
        x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame)))
    }
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

