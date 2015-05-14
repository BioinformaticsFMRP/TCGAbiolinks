#' Download data of samples from TCGA
#'
#' TCGAbiolinks allows you to Download data of samples from TCGA
#'
#' The functions you're likely to need from \pkg{TCGAbiolinks} is
#' \code{\link{tcgaDownload}}, \code{\link{tcgaQuery}} and
#' \code{\link{tcgaUpdate}}. Otherwise refer to the vignettes to see
#' how to format the documentation.
#'
#' @docType package
#' @name TCGAbiolinks
NULL

#' TCGA disease table
#' @docType data
#' @keywords internal
#' @name disease.table
#' @format A data frame with 37 rows and 4 variables
NULL

#' TCGA platforms table
#' @docType data
#' @keywords internal
#' @name platform.table
#' @format A data frame with 79 rows and 4 variables
NULL


#' TCGA center table
#' @docType data
#' @keywords internal
#' @name center.table
#' @format A data frame with 36 rows and 3 variables
NULL

#' TCGA metadata database.
#' The data set contains the following fields:
#' \itemize{
#' \item addedDate Date when sample was added to database
#' \item baseName name of the sample folder
#' \item deployLocation Path of the sample folder
#' \item barcode list of barcode
#' \item id TCGA id
#' \item isLatest Is the lasted version of the sample?
#' \item name Sample name
#' \item revision Sample revision
#' \item serialIndex Sample serial Index
#' \item Center Center name
#' \item Platform Sample platform
#' \item Disease Sample disease
#'}
#' @docType data
#' @keywords internal
#' @name tcga.db
#' @format A data frame with 4738 rows and 12 variables
NULL
