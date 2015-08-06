#' @title Filtering common samples among platforms from TCGAquery for the same tumor
#' @description In order to help the user to have an overview of the number of
#' samples in commun we created the function `TCGAquery_integrate` that will receive the
#' data frame returned from `TCGAquery` and produce a matrix n platforms x n platforms
#' with the values of samples in commum.
#' @param query is the output of TCGAquery
#' @export
#' @return table with common samples among platforms from TCGAquery
#' @examples
#' query <- TCGAquery(tumor = 'brca',level = 3)
#' matSamples <- TCGAquery_integrate(query)
TCGAquery_integrate <- function(query) {

    querySamples <- TCGAquery_samplesfilter(query)

    matSamples <- matrix(0, length(querySamples), length(querySamples))
    colnames(matSamples) <- names(querySamples)
    rownames(matSamples) <- names(querySamples)

    for (i in 1:nrow(matSamples)) {
        for (j in 1:nrow(matSamples)) {
            if (i != j) {
                x <- as.vector(substr(querySamples[[i]], 1, 12))
                y <- as.vector(substr(querySamples[[j]], 1, 12))
                matSamples[i, j] <- length(intersect(x, y))
            } else {
                matSamples[i, j] <- length(querySamples[[i]])
            }
        }
    }
    return(matSamples)
}

#' @title Filtering sample output from TCGAquery
#' @description
#'    Filtering sample output from TCGAquery
#' @param query metaData output from TCGAquery
#' @examples query <- TCGAquery(tumor = 'brca',level = 3)
#' querySamples <- TCGAquery_samplesfilter(query)
#' @export
#' @return list of samples for a tumor
TCGAquery_samplesfilter <- function(query) {

    # Find unique platforms
    plat <- sort(unique(query$Platform))
    TumorDataList <- vector("list", length(plat))
    names(TumorDataList) <- plat

    for (idx in 1:length(TumorDataList)) {
        currPlatform <- query[query$Platform == names(TumorDataList)[idx], ]
        TumorDataList[[idx]] <- as.matrix(
            unlist(strsplit(currPlatform$barcode, ","))
        )
    }
    return(TumorDataList)
}
