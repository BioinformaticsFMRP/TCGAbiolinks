#' @title TCGA Download
#' @description Download data previously selected using the TCGAQuery function
#' @param dirURL the location of the query saved URL links.
#' @param finalDir location of the final data saving
#' @param earlyStop If the number of file is too large
#'                 for you need you can choose a subset
#' @author Davide
#' @seealso TCGAQuery
#' @examples
#' \dontrun{
#'          TCGADownload(dirURL="data/query/",
#'                      finalDir = "data/final",
#'                      earlyStop = 0)
#' }
#' @export
#' @import downloader
TCGADownload <- function(dirURL="data/query/", finalDir="data/final", earlyStop = 0){

  # load query files if not yet done
  if(!exists("queryURI")) load(paste0(dirURL,"fileURLs.rda"))
  queryURI <- get("queryURI", envir=environment())

  finalDir <- createDir(finalDir)
  t1 = Sys.time()
  if(earlyStop==0 || earlyStop > length(queryURI)) earlyStop <- length(queryURI)
  for(k in 1:earlyStop){
    dir <- paste(finalDir, basename(dirname(queryURI[k])),sep="/")
    dir.create(dir, showWarnings = F)
    download(queryURI[k],
             destfile = paste(dir,
                              basename(queryURI[k]),
                              sep="/"),
             mode="w",
             quiet = 1)
    print(paste("Downloaded",k,"out of",length(queryURI),sep=" "))
    if(k%%10==0) print(paste("Estimated time for the end of the process:",
                             abs((as.numeric(difftime(Sys.time(), t1, units="min"))/k)*earlyStop -
                               as.numeric(difftime(Sys.time(), t1, units="min"))), "minutes", sep=" "))
  }
  rm(t1)
  message(paste0('Donwload completed. Files saved in "', finalDir,'"'))
}
