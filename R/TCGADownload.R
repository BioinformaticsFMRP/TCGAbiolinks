#' @title TCGA Download
#' @description Download data previously selected using the TCGAQuery function
#' @param data TCGAQUery output
#' @param path location of the final data saving
#' @seealso TCGAQuery
#' @examples
#' \dontrun{
#'          TCGADownload(data,"folder")
#' }
#' @export
#' @import downloader
TCGADownload <- function(data=NULL,path=".")
{
  dir.create(path,showWarnings = F)
  root <- "https://tcga-data.nci.nih.gov"
  if(!("file" %in% colnames(data))){
    message("Downloading folders")
    for(i in 1:nrow(data)){
      file <- paste0(path,"/",basename(data[i,"deployLocation"]))
      message(paste0("Downloading:", basename(data[i,"deployLocation"])))
      if(!file.exists(file)){
        downloader::download(paste0(root,data[i,"deployLocation"]),file)
        untar(file, exdir = path)
      }
    }
  }
  else{
    message("Downloading files")
    for(i in 1:nrow(data)){
      file <- paste0(path,"/",basename(data[i,"file"]))
      message(paste0("Downloading:", basename(data[i,"file"])))
      if(!file.exists(file)){
        downloader::download(paste0(root,gsub(".tar.gz","",data[i,"deployLocation"]),"/",data[i,"file"]),file)
      }
    }
  }
}
