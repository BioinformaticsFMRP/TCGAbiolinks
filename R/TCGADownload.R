#' @title TCGA Download
#' @description Download data previously selected using the TCGASeach
#' @param data TCGASearch output
#' @param path location of the final data saving
#' @param type Get files with type pattern instead of downloading all the files
#' @param quiet Should download supress output messages. Default: false
#' @param samples List of samples to download
#' @seealso TCGASearch
#' @examples
#'    samples <- c("TCGA-06-0125-01A-01D-A45W-05")
#'    query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",
#'    level = "3", samples = samples)
#'    TCGADownload(query,path = "dataDemo2")
#' @export
#' @importFrom downloader download
#' @return Download tcga into path
TCGADownload <- function(data = NULL, path = ".", type = NULL, samples = NULL,
                         quiet = FALSE) {

    dir.create(path, showWarnings = FALSE)
    root <- "https://tcga-data.nci.nih.gov"

    # Downloading the folder
    if (is.null(type) && is.null(samples) ) {
        for (i in 1:nrow(data)) {

            file <- paste0(path, "/", basename(data[i, "deployLocation"]))
            message(paste0("Downloading:",
                           basename(data[i, "deployLocation"])))
            if (!file.exists(file)) {
                download(paste0(root, data[i, "deployLocation"]),
                         file, quiet,method="auto")
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

            if(!is.null(type)){
                files <- files[grepl(type,files)]
                if (length(files) == 0) {
                    next
                }
            }

            if(!is.null(samples)){
                files <- filterFiles(data[i,],samples,files)
            }
            message(paste0("Downloading:", length(files), " files"))

            for (i in seq_along(files)) {
                if (!file.exists(file.path(path,folder,files[i]))) {
                    download(paste0(root,url,"/",files[i]),
                             file.path(path,folder,files[i]),quiet,method="auto")
                }
            }
        }
    }
}

# Filter files by barcode
filterFiles <- function(data,samples,files){

    # case uuid in name
    if(grepl("IlluminaHiSeq_RNASeqV2",data$Platform)){
        mage <- getMage(data)
        idx <- unlist(lapply(samples,
                             function(x) grep(x,mage$Comment..TCGA.Barcode.)))
        idx <- unique(idx)
        names <- mage[idx,]$Extract.Name
        idx <- unique(unlist(lapply(names, function(x) grep(x,files))))
        files <- files[idx]
        return(files)
    }
    # case barcode in name
    if(grepl("IlluminaHiSeq_RNASeq|humanmethylation",data$Platform,
             ignore.case = TRUE)){
        idx <- unique(unlist(lapply(samples, function(x) grep(x,files))))
        files <- files[idx]
        return(files)
    }

}
