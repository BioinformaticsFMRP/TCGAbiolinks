#' @title TCGA Download
#' @description Download data previously selected using the TCGASeach
#' @param data TCGAQuery output
#' @param path Directory to save the downloaded data
#' @param type Filter the files that will be downloaded by type
#' @param quiet Supress output messages?. Default: FALSE
#' @param samples List of samples to download
#' @param force Download files even if it was already downladed? Default: FALSE
#' @seealso TCGAQuery
#' @examples
#'    samples <- c("TCGA-26-1442-01A-01R-1850-01")
#'    query <- TCGAQuery(tumor = "gbm", platform = "IlluminaHiSeq_RNASeqV2",
#'    level = "3", samples = samples)
#'    TCGADownload(query,path = "dataDemo2",samples = samples,
#'                 type ="rsem.genes.results")
#' @export
#' @importFrom downloader download
#' @return Download tcga into path
TCGADownload <- function(data = NULL, path = ".", type = NULL, samples = NULL,
                         quiet = FALSE, force = FALSE) {

    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    root <- "https://tcga-data.nci.nih.gov"

    # Downloading the folder
    if (is.null(type) && is.null(samples) ) {
        for (i in 1:nrow(data)) {

            file <- paste0(path, "/", basename(data[i, "deployLocation"]))
            cat(paste0("Downloading:",
                       basename(data[i, "deployLocation"]),"\n"))
            if (force || !file.exists(file)) {
                if(is.windows()){
                    suppressWarnings(
                        download(paste0(root, data[i, "deployLocation"]),
                                 file, quiet, method = "auto")
                    )
                } else {
                    download(paste0(root, data[i, "deployLocation"]),
                             file, quiet)

                }
                untar(file, exdir = path)
            }
        }
    } else {
        # Downloading files
        for (i in 1:nrow(data)) {

            folder <- gsub(".tar.gz","",basename(data[i,]$deployLocation))
            dir.create(file.path(path,folder), showWarnings = FALSE)
            url <- gsub(".tar.gz","",data[i,]$deployLocation)
            files <- getFileNames(paste0(root,url))
            idx <- grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE|Name|Size|Parent|Last",files)
            files <- files[-idx]

            if(!is.null(type)){
                files <- files[grepl(type,files)]
                if (length(files) == 0) {
                    next
                }
            }

            if(!is.null(samples)){
                files <- filterFiles(data[i,],samples,files)
            }
            cat(paste0("Downloading:", length(files), " files","\n"))

            for (i in seq_along(files)) {
                if (force || !file.exists(file.path(path,folder,files[i]))) {
                    if(is.windows()){
                        suppressWarnings(
                            download(paste0(root,url,"/",files[i]),
                                     file.path(path,folder,files[i]),
                                     quiet,method = "auto"
                            )
                        )
                    } else {
                        download(paste0(root,url,"/",files[i]),
                                 file.path(path,folder,files[i]),quiet)
                    }
                }
            }
        }
    }
}

# Filter files by barcode
filterFiles <- function(data,samples,files){

    barcodeName <- paste("IlluminaHiSeq_RNASeq",
                         "humanmethylation",
                         "H-miRNA_8x15K",
                         "images",
                         "SOLiD_DNASeq",
                         "pathology_reports",
                         "IlluminaDNAMethylation",
                         "HG-CGH-244A",
                         "HG-CGH-415K_G4124A",
                         "HG-U133_Plus_2",
                         "IlluminaGA_DNASeq_automated",
                         "IlluminaGA_miRNASeq",
                         "IlluminaGA_mRNA_DGE",
                         "IlluminaGA_RNASeq",
                         "IlluminaHiSeq_DNASeqC",
                         "IlluminaHiSeq_miRNASeq",
                         "IlluminaHiSeq_RNASeq",
                         "IlluminaHiSeq_WGBS", sep = "|")

    uuidName <- paste("RNASeqV2",
                      "MDA_RPPA_Core",
                      sep = "|")

    mageName <-  paste("AgilentG4502A",
                       "CGH-1x1M_G4447A",
                       "Genome_Wide_SNP_6",
                       "HT_HG-U133A",
                       sep = "|")


    # case uuid in name
    if(grepl(uuidName,data$Platform, ignore.case = TRUE)){
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        uuid <- str_match(files,regex)[,1]
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        idx <- unique(unlist(lapply(samples, function(x) grep(x,map$barcode))))
        idx <- which(uuid %in% map[idx,]$uuid)
        files <- files[idx]
    } else if (grepl("IlluminaGA_DNASeq_curated",data$Platform) & data$Center == "broad.mit.edu"){
    # Exception - two uuids in the name
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        files <- files[grep(regex,files)]
        uuid <- unlist(str_match_all(files,regex))
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        idx <- unique(unlist(lapply(samples, function(x) grep(x,map$barcode))))
        idx <- which(uuid %in% map[idx,]$uuid)
        files <- files[ceiling(idx / 2)]
    } else if(grepl(barcodeName, data$Platform, ignore.case = TRUE)) {
        idx <- unique(unlist(lapply(samples, function(x) grep(x,files))))
        files <- files[idx]
    } else if(grepl(mageName, data$Platform, ignore.case = TRUE)) {
        mage <- getMage(data)
        idx <- unlist(lapply(samples,
                             function(x) grep(x,mage$Comment..TCGA.Barcode.)))
        idx <- unique(idx)
        mage <- mage[idx,]
        idx <- grep("Derived.Array.Data.Matrix.File|Array.Data.File",colnames(mage))
        names <- unique(unlist(mage[,idx]))
        idx <- unique(unlist(lapply(names, function(x) grep(x,files))))
        files <- files[idx]
    }
    return(files)
}
