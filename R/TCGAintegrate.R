#' @title TCGAintegrate
#' @description Common samples among platforms with 2 by 2 comparison
#' @param query output of TCGAQuery
#' @export
#' @return table with common samples among platforms from TCGAQuery
#' @examples
#' query <- TCGAQuery(tumor = "brca",level = 3)
#' matSamples <- TCGAintegrate(query)
TCGAintegrate <- function(query){

    querySamples <- samplesfilter(query)

    matSamples <- matrix(0,length(querySamples),length(querySamples))
    colnames(matSamples) <- names(querySamples)
    rownames(matSamples) <- names(querySamples)

    for (i in 1:nrow(matSamples)) {
        for (j in 1:nrow(matSamples)) {
            if (i != j) {
                x <- as.vector(substr(querySamples[[i]],1,12))
                y <- as.vector(substr(querySamples[[j]],1,12))
                matSamples[i,j] <- length(intersect(x,y))
            } else {
                matSamples[i,j] <- length(querySamples[[i]])
            }
        }
    }
    return(matSamples)
}

#' @title Filtering sample output from TCGAQuery
#' @description
#'    Filtering sample output from TCGAQuery
#' @param query metaData output from TCGAQuery
# @examples
# query <- TCGAQuery(tumor = "brca",level = 3)
# querySamples <- samplesfilter(query)
# @export
# @return list of samples for a tumor
samplesfilter <- function(query) {

    # Find unique platforms
    plat <- sort(unique(query$Platform))
    TumorDataList <- vector("list", length(plat))
    names(TumorDataList) <- plat

    for (idx in 1:length(TumorDataList)) {
        currPlatform <- query[query$Platform == names(TumorDataList)[idx],]
        TumorDataList[[idx]] <- as.matrix(unlist(strsplit(currPlatform$barcode,",")))
    }
    return(TumorDataList)
}

