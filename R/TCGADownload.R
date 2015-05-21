#' @title TCGA Download
#' @description Download data previously selected using the TCGASeach
#' @param data TCGASearch output
#' @param path location of the final data saving
#' @param type Get files with type pattern instead of downloading all the folder
#' @seealso TCGASearch
#' @examples
#'    samples <- c("TCGA-06-0125-01A-01D-A45W-05")
#'    query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",
#'    level = "3", samples = samples)
#'    TCGADownload(query,path = "dataDemo2")
#' @export
#' @importFrom downloader download
#' @return Download tcga into path
TCGADownload <- function(data = NULL, path = ".", type = NULL, samples = NULL) {

  dir.create(path, showWarnings = FALSE)
  root <- "https://tcga-data.nci.nih.gov"

  # Downloading the folder
  if (is.null(type)) {
    for (i in 1:nrow(data)) {

      file <- paste0(path, "/", basename(data[i, "deployLocation"]))
      message(paste0("Downloading:",
                     basename(data[i, "deployLocation"])))
      if (!file.exists(file)) {
        download(paste0(root, data[i, "deployLocation"]),
                 file)
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
      files <- files[grepl(type,files)]
      if (length(files) == 0) {
        next
      }

      if(!is.null(samples)){
          files <- filterFiles(data[i,],samples,files)
      }
      message(paste0("Downloading:", length(files), " files"))

      for (i in seq_along(files)) {
        if (!file.exists(files[i])) {
          download(paste0(root,url,"/",files[i]),
                   file.path(path,folder,files[i]))
        }
        #message(paste0("Downloading:", files[i], " of ", length(files), " files"))

      }
    }
  }
}

# Filter files by barcode
filterFiles <- function(data,samples,files){

    if(grep("IlluminaHiSeq",data$Platform)){
        mage <- getMage(data)
        idx <- unique(unlist(lapply(samples,
                    function(x) grep(x,mage$Comment..TCGA.Barcode.))))
        names <- mage[idx,]$Extract.Name
        idx <- unique(unlist(lapply(names,function(x) grep(x,files))))
        files <- files[idx]
        return(files)
    }
}
