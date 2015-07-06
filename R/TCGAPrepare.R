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
#' @param save Save a rda object with the prepared object? Default FALSE
#' @param filename Name of the saved file
#' @param toPackage For whihc package are you preparing the data?
#' Default: NULL. Options: "ELMER"
#' @examples
#' sample <- "TCGA-06-0939-01A-01D-1228-05"
#' query <- TCGAQuery(tumor = "GBM",samples = sample, level = 3)
#' TCGADownload(query,path = "exampleData",samples = sample, quiet = TRUE)
#' prepared <- TCGAPrepare(query, dir="exampleData")
#' @export
#' @importFrom stringr str_match str_trim
#' @import data.table
#' @import utils
TCGAPrepare <- function(query,
                        dir = NULL,
                        type = NULL,
                        save = FALSE,
                        filename = NULL,
                        toPackage = NULL){

    gene.location   <- get("gene.location",
                           envir =  as.environment("package:TCGAbiolinks"))

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

    # Get all files from directory except MANIFEST, README, etc
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

    # Filter by type
    if (!is.null(type)) {
        files <- files[grep(type,files)]
        if(length(files) == 0){
            message("No files of that type found")
            return(NULL)
        }
    }

    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    df <- NULL


    if (grepl("humanmethylation",tolower(platform))) {

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            sample <- colnames(data)[2]
            setnames(data,gsub(" ", "\\.", data[1,]))
            data <- data[-1,] # removing Composite Element REF
            setnames(data,2,sample)

            if (i == 1) {
                setcolorder(data,c(1, 3:5, 2))
                df <- data
            } else {
                data <- subset(data,select = c(1,2))
                df <- merge(df, data, by = "Composite.Element.REF")
            }

            setTxtProgressBar(pb, i)
        }

        setDF(df)
        rownames(df) <- df$Composite.Element.REF
        df$Composite.Element.REF <- NULL
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
            setTxtProgressBar(pb, i)
        }
        rownames(df) <- df[,1]
        df[,1] <- NULL
        # get array_design.txt from mage folder
        # and change uuid by Barcode
        uuid <- colnames(df)
        idx <- grep("Sample|Control",uuid)
        if(length(idx) > 0){
            uuid <- uuid[-idx]
        }
        map <- mapuuidbarcode(uuid)
        idx <- which(colnames(df) %in% map$uuid)
        colnames(df)[idx] <- as.character(map$barcode)
    }


    if (grepl("illuminadnamethylation_oma",
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
            setTxtProgressBar(pb, i)
        }

        rownames(df) <- df[,1]
        df[,1] <- NULL
    }

    if (tolower(platform) == "illuminaga_rnaseq" ||
        tolower(platform) == "illuminahiseq_rnaseq") {
        # Barcode in the name
        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- str_match(files,regex)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)

            x <- colnames(data)
            setnames(data,colnames(data)[2:4],paste0(colnames(data)[2:4],"_",barcode[i]))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            setTxtProgressBar(pb, i)
        }
    }

    if (tolower(platform) == tolower("HT_HG-U133A")) {
        # Barcode in the mage file
        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t", skip= 1,
                          stringsAsFactors = FALSE)

            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            setTxtProgressBar(pb, i)
        }
        names <- sapply(files,
                        function(x) {
                            y <- fread(x, header = FALSE,
                                       stringsAsFactors = FALSE,
                                       nrows=1)$V2
                            idx <- grep(y,mage$Hybridization.Name)
                            mage[idx,]$Comment..TCGA.Barcode.
                        })
        setnames(df,2:ncol(df),names)
    }

    if (tolower(platform) == tolower("HG-U133_Plus_2") ||
        grepl("H-miRNA_8x15K|agilent",platform, ignore.case = TRUE)) {
        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)

            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            setTxtProgressBar(pb, i)
        }
        df <- df[-1,]

    }
    if (grepl("rnaseqv2",platform, ignore.case = TRUE)) {

        if(is.null(type) || (type != "rsem.genes.results" &&
                             type != "rsem.isoforms.results" &&
                             type != "rsem.genes.normalized_results" &&
                             type != "rsem.isoforms.normalized_results")
           ){
            msg <- paste0("Plase select a type. \n Possibilities:\n",
                          " = rsem.genes.results\n",
                          " = rsem.isoforms.results\n",
                          " = rsem.genes.normalized_results\n",
                          " = rsem.isoforms.normalized_results\n"
                          )
            message(msg)
            return()
        }

        if(type == "rsem.genes.results")               pat <- "rsem.genes.results"
        if(type == "rsem.isoforms.results")            pat <- "rsem.isoforms.results"
        if(type == "rsem.genes.normalized_results")    pat <- "rsem.genes.normalized_results"
        if(type == "rsem.isoforms.normalized_results") pat <- "rsem.isoforms.normalized_results"

        files <- files[grep(pat,files, perl = TRUE)]

        regex <- paste0("[:alnum:]{8}-[:alnum:]{4}",
                        "-[:alnum:]{4}-[:alnum:]{4}-[:alnum:]{12}")
        uuid <- str_match(files,regex)
        map <- mapuuidbarcode(uuid)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            x <- subset(map, uuid == uuid[i])
            setnames(data,2, as.character(x$barcode))
            data <- subset(data,select=c(1,2))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = colnames(df)[1])
            }
            setTxtProgressBar(pb, i)
        }
        setDF(df)
        rownames(df) <- df[,1]
        df[,1] <- NULL
    }

    if (grepl("illuminahiseq_mirnaseq",platform, ignore.case = TRUE)) {

        if(is.null(type) || (type != "hg19.mirna" && type != "mirna")){
            msg <- paste0("Plase select a type. \n Possibilities:\n",
                          " = hg19.mirna\n = mirna")
            message(msg)
            return()
        }

        if(type == "hg19.mirna")   pat <- "(hg19.)mirna"
        if(type == "mirna")        pat <- "(?<!hg19\\.)mirna"

        files <- files[grep(pat,files, perl = TRUE)]

        if(length(files) == 0){
            message("No mirna files of that type found")
            return(NULL)
        }

        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- str_match(files,regex)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            data <- subset(data,select=c(1:3))
            setnames(data,2:ncol(data),
                     paste0(as.character(barcode[i]),".",colnames(data)[2:ncol(data)]))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, by=colnames(data)[1])
            }
        }
        setDF(df)
        rownames(df) <- df[,1]
        df <- df[,-1]
    }

    if (grepl("bio",platform,ignore.case = TRUE)) {
        if (!is.null(type)) {
            files <- files[grep(type,files)]
        }
        if (length(files) == 1) {
            df <- read.table(files, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE,
                             comment.char = "#",fill = TRUE, skip =1)
            regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                            "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
            idx <- grepl(regex,df$bcr_patient_uuid)
            df <- df[idx,]
        } else {
            message("We're preaparing for the moment only one clinical file")
            return(NULL)
        }
    }
    if (grepl("genome_wide_snp_6",tolower(platform))) {

        close(pb)
        message("Preparing h19 files...")
        idx <- grep("nocnv|hg18", files)
        if(length(idx)>0){
            files <- files[-idx]
        }

        pb <- txtProgressBar(min = 0, max = length(files), style = 3)

        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        #load("hg19genes.RData")
        genes <- sort(unique(gene.location$external_gene_name))

        df <- matrix(0, nrow = length(genes), ncol = length(files))
        rownames(df) <- genes

        colNames <- rep("", ncol(df)) #check barcode
        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            ####Check barcodes
            colNames[i] <- data$Sample[1]

            for(j in 1:nrow(data)){
                gg <- sort(unique(gene.location[gene.location$start_position >= data$Start[j]
                                                & gene.location$end_position <= data$End[j],
                                                "external_gene_name"]))
                df[gg, i] <- df[gg, i] + data$Segment_Mean[j]
            }
            setTxtProgressBar(pb, i)
        }
        id <- data.frame(id=colNames)
        #idx <- unlist(lapply(colNames,function(x){grep(x, mage$Derived.Array.Data.Matrix.File.2)}))
        names <- merge(id,mage,by.x="id",by.y="Hybridization.Name")
        colnames(df) <- names$Comment..TCGA.Barcode.
    }
    close(pb)

    if(save){
        if(is.null(filename)){
            filename <- paste0(platform,"_",gsub(" ","_",Sys.time()),".rda")
        }
        save(df,file = filename)
    }


    if(!is.null(toPackage)){
        df <- prepareToPackage(df, platform,toPackage)
    }

    return(df)
}

# This function will help the user to prepare the data to an specific package
prepareToPackage <- function(df, platform, toPackage){

    if(grepl("elmer", toPackage, ignore.case = TRUE)){

        if (grepl("illuminahiseq_rnaseqv2|illuminahiseq_totalrnaseqv2",
                  platform, ignore.case = TRUE)){
            message("============ Pre-pocessing expression data =============")
            message(paste0("1 - expression = log2(expression + 1): ",
                           "To linearize \n    relation between ",
                           "methylation and expression"))
            df <- log2(df+1)
            message("2 - rownames  (gene|loci) => ('ID'loci) ")
            aux <- strsplit(rownames(df),"\\|")
            GeneID <- unlist(lapply(aux,function(x) x[2]))
            row.names(df) <- paste0("ID",GeneID)
            df <- as.matrix(df)
        }

        if (grepl("humanmethylation", platform, ignore.case = TRUE)) {
            message("============ Pre-pocessing methylation data =============")
            msg <- paste0("1 - Removing Columns: \n  * Gene_Symbol  \n",
                          "  * Chromosome  \n  * Genomic_Coordinate")
            message(msg)
            df <- subset(df,select = 4:ncol(df))
            msg <- paste0("2 - Removing probes with ",
                          "NA values in more than 0.80% samples")
            message(msg)
            df <- df[rowMeans(is.na(df))<0.2,]
        }
    }

    if(grepl("limma", toPackage, ignore.case = TRUE)){
        message("TBD")
    }

    if(grepl("methylmix", toPackage, ignore.case = TRUE)){
        message("TBD")
    }

    if(grepl("biomics", toPackage, ignore.case = TRUE)){
        message("TBD")
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
            if(!is.windows()){
                download(paste0(root,mage$deployLocation), file,
                         quiet = TRUE)
            } else {
                download(paste0(root,mage$deployLocation), file,
                         quiet = TRUE,method="auto")
            }
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

