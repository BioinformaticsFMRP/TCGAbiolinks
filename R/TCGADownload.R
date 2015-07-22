#' @title Download the data from TCGA using as reference the output from TCGAquery
#' @description
#'      The TCGAdownload function will download the data using as reference
#'      the the lines of the TCGAquery search result.
#'
#'      There is an option to download the entire tar.gz folder or download
#'      specific files using the \emph{type} parameter or the \emph{samples}
#'      parameter
#'
#'      The outpufiles will be saved into the path parameters. If this path does
#'      not exists the package will try to create the directories.
#'
#'      By default, if a sample was already downloaded the function will not
#'      download again, unless the force parameter is set to \code{TRUE}
#'
#' @param data The TCGAquery output
#' @param path Directory to save the downloaded data
#' @param type Filter the files that will be downloaded by
#'  type. Example:"rsem.genes.results"
#' @param samples List of samples to download data
#' @param force Download files even if it was already downladed?
#' Default: \code{FALSE}
#' @seealso \code{\link{TCGAquery}} for searching the data to download
#'
#' \code{\link{TCGAprepare}} for preparing the data for the user into
#' a Summarized experiment object, or a matrix.
#' @examples
#' samples <- c("TCGA-26-1442-01A-01R-1850-01")
#' query <- TCGAquery(tumor = "gbm",
#'                    platform = "IlluminaHiSeq_RNASeqV2",
#'                    level = "3",
#'                    samples = samples)
#' TCGAdownload(query,path = "RNA",
#'              samples = samples,
#'              type ="rsem.genes.results")
#' @export
#' @importFrom downloader download
#' @return Download TCGA data into the given path
#' @family data functions
TCGAdownload <- function(data = NULL, path = ".", type = NULL, samples = NULL,
                         force = FALSE) {

    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    root <- "https://tcga-data.nci.nih.gov"

    # Downloading the folder
    if (is.null(type) && is.null(samples) ) {
        message(rep("-=",(nchar(file.path(path))/2) + 4))
        message(paste0("| Downloading:", dim(data), " folders"))
        message(paste0("| Path:", file.path(path)))
        message(rep("-=",(nchar(file.path(path))/2) + 4))
        pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)
        for (i in 1:nrow(data)) {

            file <- paste0(path, "/", basename(data[i, "deployLocation"]))
            cat(paste0("\nDownloading:",
                       basename(data[i, "deployLocation"]),"\n"))
            if (force || !file.exists(file)) {
                if(is.windows()){
                    suppressWarnings(
                        download(paste0(root, data[i, "deployLocation"]),
                                 file, quiet = TRUE, method = "auto")
                    )
                } else {
                    suppressWarnings(
                        download(paste0(root, data[i, "deployLocation"]),
                                 file, quiet = TRUE)
                    )

                }
                untar(file, exdir = path)
            }
            setTxtProgressBar(pb, i)
        }
        close(pb)
    } else {
        # Downloading files

        for (i in 1:nrow(data)) {

            folder <- gsub(".tar.gz","",basename(data[i,]$deployLocation))

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

            if(length(files) > 0){
                dir.create(file.path(path,folder), showWarnings = FALSE)

                message(rep("-=",(nchar(file.path(path,folder))/2)+4))
                message(paste0("| Downloading:", length(files), " files"))
                message(paste0("| Path:", file.path(path,folder)))
                message(rep("-=",(nchar(file.path(path,folder))/2)+4))
                pb <- txtProgressBar(min = 0, max = length(files), style = 3)
            }

            for (i in seq_along(files)) {
                if (force || !file.exists(file.path(path,folder,files[i]))) {
                    message(paste0("[",i,"] ", files[i],"\n"))
                    if(is.windows()){
                        suppressWarnings(
                            download(paste0(root,url,"/",files[i]),
                                     file.path(path,folder,files[i]),
                                     quiet = TRUE,method = "auto"
                            )
                        )
                    } else {
                        suppressWarnings(
                            download(paste0(root,url,"/",files[i]),
                                     file.path(path,folder,files[i]),
                                     quiet = TRUE)
                        )
                    }
                }
                setTxtProgressBar(pb, i)
            }
            if(length(files) > 0){
                close(pb)
            }
        }
    }

}

# Filter files by barcode
filterFiles <- function(data,samples,files){

    # If it if maf it is better to let download all files
    #maf <- paste("SOLiD_DNASeq_curated",
    #             "SOLiD_DNASeq", # partial barcode 2
    #             "illuminaga_dnaseq",
    #             sep = "|"
    #)

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
                         "IlluminaHiSeq_RNASeq", sep = "|")

    uuidName <- paste("RNASeqV2",
                      "MDA_RPPA_Core",
                      sep = "|")

    mageName <-  paste("AgilentG4502A",
                       "CGH-1x1M_G4447A",
                       "Genome_Wide_SNP_6",
                       "HT_HG-U133A",
                       "IlluminaHiSeq_WGBS",
                       sep = "|")


    if (grepl(uuidName,data$Platform, ignore.case = TRUE)) {
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        uuid <- str_match(files,regex)[,1]
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        idx <- unique(unlist(lapply(samples, function(x) grep(x,map$barcode))))
        idx <- which(tolower(uuid) %in% map[idx,]$uuid)
        files <- files[idx]
    } else if (grepl("IlluminaGA_DNASeq_curated|illuminaga_dnaseq_automated",data$Platform) & data$Center == "broad.mit.edu") {
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
        idx <- grep("Derived.Array.Data.Matrix.File|Array.Data.File|Derived.Data.File",
                    colnames(mage))
        names <- unique(unlist(mage[,idx]))
        idx <- unique(unlist(lapply(names, function(x) grep(x,files))))
        files <- files[idx]
    }
    return(files)
}
