#' @title TCGAintegrate
#' @description Common samples among platforms with 2 by 2 comparison
#' @param querySamples output of TCGAQuery filtered
#' @export
#' @return table with common samples among platforms from TCGAQuery
#' @examples
#' TCGAintegrate(c("TCGA-A7-A0D9-01A-31D-A060-02",
#'                 "TCGA-B6-A0RG-01A-11D-A060-02",
#'                  "TCGA-B6-A0RG-01A-11D-A060-01"))
TCGAintegrate<- function(querySamples){

    matSamples <- matrix(0,length(querySamples),length(querySamples))
    colnames(matSamples) <- names(querySamples)
    rownames(matSamples) <- names(querySamples)

    for ( i in 1: nrow(matSamples)){
        for( j in 1: nrow(matSamples)){
            if(i!=j){ matSamples[i,j] <- length(
                intersect(substr(querySamples[[i]],1,12),
                          substr(querySamples[[j]],1,12)))}
        else{ matSamples[i,j] <- length(intersect(
            querySamples[[i]], querySamples[[j]]))}
        }

    }
    return(matSamples)
}
