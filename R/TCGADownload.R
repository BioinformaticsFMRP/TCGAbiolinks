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
  for(i in 1:nrow(data)){
    file <- paste0(path,"/",basename(data[i,"deployLocation"]))
    message(paste0("Downloading:", basename(data[i,"deployLocation"])))
    downloader::download(paste0(root,x[i,"deployLocation"]),file)
    untar(file, exdir = path)
  }
}
