#' @title Prepare CEL files into an AffyBatch.
#' @description Prepare CEL files into an AffyBatch.
#' @param ClinData
#' @param PathFolder
#' @param TabCel
#' @importFrom affy ReadAffy
#' @importFrom affy rma
#' @importFrom Biobase exprs
#' @examples
#' \dontrun{
#' to add example
#' }
#' @export
#' @return Normalizd Expression data from Affy eSets
TCGAprepare_Affy <- function(ClinData, PathFolder, TabCel){

    affy_batch <- ReadAffy(filenames=as.character(paste(TabCel$samples, ".CEL", sep="")))

    eset <- rma(affy_batch)

    mat <- exprs(eset)

    return(mat)

    }
