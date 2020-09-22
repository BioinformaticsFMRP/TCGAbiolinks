#' TCGA samples with their Tumor Purity measures
#'
#' A dataset containing the Sample Ids from TCGA tumor purity measured according to 4 estimates attributes of 9364 tumor patients
#'
#' @format A data frame with 9364 rows and 7 variables:
#' \describe{
#'   \item{Sample.ID}{Sample ID from TCGA barcodes, character string}
#'   \item{Cancer.type}{Cancer type, character string}
#'   \item{ESTIMATE}{uses gene expression profiles of 141 immune genes and 141 stromal genes, 0-1 value}
#'   \item{ABSOLUTE}{uses somatic copy-number data (estimations were available for only 11 cancer types), 0-1 value}
#'   \item{LUMP}{(leukocytes unmethylation for purity), which averages 44 non-methylated immune-specific CpG sites, 
#' 0-1value}
#'   \item{IHC}{as estimated by image analysis of haematoxylin and eosin stain slides produced by the Nationwide 
#' Childrens Hospital Biospecimen Core Resource, 0-1 value}
#'   \item{CPE}{derived consensus measurement as the median purity level after normalizing levels from all methods to 
#' give them equal means and s.ds, 0-1 value}
#'   ...
#' }
#' @source \url{https://images.nature.com/original/nature-assets/ncomms/2015/151204/ncomms9971/extref/ncomms9971-s2.xlsx}

"Tumor.purity"
