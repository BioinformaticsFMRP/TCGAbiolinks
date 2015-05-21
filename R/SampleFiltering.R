#' @title Filtering sample output from TCGAQuery
#' @description
#'    Filtering sample output from TCGAQuery
#' @param query metaData output from TCGAQuery
#' @examples
#' query <- TCGAQuery(tumor = "brca",level = 3)
#' querySamples <- samplesfilter(query)
#' @export
#' @return list of samples for a tumor
samplesfilter <- function(query) {

    # Find unique platforms
    Platforms <- sort(unique(query$Platform))
    TumorDataList <- vector("list", length(sort(unique(query$Platform))))
    names(TumorDataList) <- sort(unique(query$Platform))

    for (idx in 1:length(TumorDataList)) {
        currPlatform <- query[query$Platform == names(TumorDataList)[idx],]
        currPlatform <- currPlatform[order(currPlatform$addedDate,
                                           decreasing = TRUE),]
        # taking last version updated
        sampleUpdated <- currPlatform[1,]
        TumorDataList[[idx]] <- as.matrix(unlist(
            strsplit(sampleUpdated$barcode,",")))
    }
    return(TumorDataList)
}
