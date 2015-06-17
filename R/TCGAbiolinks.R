#' Download data of samples from TCGA
#'
#' TCGAbiolinks allows you to Download data of samples from TCGA
#'
#' The functions you're likely to need from \pkg{TCGAbiolinks} is
#' \code{\link{TCGADownload}}, \code{\link{TCGAQuery}}.
#'Otherwise refer to the vignettes to see
#' how to format the documentation.
#'
#' @docType package
#' @name TCGAbiolinks
NULL

#' DAVID GO table (BP)
#' @docType data
#' @keywords internal
#' @name DAVID_BP_matrix
#' @format A data frame with 2423 rows and 2 variables
NULL

#' DAVID GO table (CC)
#' @docType data
#' @keywords internal
#' @name DAVID_CC_matrix
#' @format A data frame with 517 rows and 2 variables
NULL

#' DAVID GO table (MF)
#' @docType data
#' @keywords internal
#' @name DAVID_MF_matrix
#' @format A data frame with 1157 rows and 2 variables
NULL

#' TCGA PlatformAndAssociatedData data frame indicating the different
#' characteristics of the data e.g. tumour, type, species.
#' @docType data
#' @keywords internal
#' @name PlatformAndAssociatedData
#' @format A data frame with 615 rows and 12 variables
NULL

#' TCGA Table with version, number of samples and size (Mbyte) of
#' BRCA IlluminaHiSeq_RNASeqV2 Level 3
#' @docType data
#' @keywords internal
#' @name BRCA_RNASeqV2_version
#' @format A data frame with 12 rows and 4 variables
NULL

#' EAGenes
#' @docType data
#' @keywords internal
#' @name EAGenes
#' @format A data frame with 20038 rows (genes) and 6 variables
NULL

#' tabPackage1
#' @docType data
#' @keywords internal
#' @name tabPackage1
#' @format A data frame with 40 rows and 3 variables
NULL

#' tabPackage2
#' @docType data
#' @keywords internal
#' @name tabPackage2
#' @format A data frame with 40 rows and 3 variables
NULL

#' wine
#' @docType data
#' @keywords internal
#' @name wine
#' @format A data frame with 178 rows and 13 variables
NULL

#' wine.class
#' @docType data
#' @keywords internal
#' @name wine.class
#' @format vector factor with with 178 variables and 3 levels
NULL

#' TCGA data matrix BRCA
#' @docType data
#' @keywords internal
#' @name dataBRCA
#' @format A data frame with 20531 rows (genes) and 50 variables (samples)
NULL

#' TCGA data matrix BRCA DEGs
#' @docType data
#' @keywords internal
#' @name dataDEGsFiltLevel
#' @format A data frame with 3649 rows and 6 variables
NULL

#' geneInfo for normalization
#' @docType data
#' @keywords internal
#' @name geneInfo
#' @format A data frame with 20531 rows and 2 variables
NULL

#' hg19genes data frame from Human (hg19) genome indicating different
#' characteristics of genes e.g. chr, region, start, end, gene_name
#' @docType data
#' @keywords internal
#' @name hg19genes
#' @format A data frame with 869204 rows and 5 variables
NULL

#' list EA pathways
#' @docType data
#' @keywords internal
#' @name listEA_pathways
#' @format A data frame with 589 rows and 4 variables
NULL

#' TCGA data matrix BRCA DEGs Pubmed
#' @docType data
#' @keywords internal
#' @name tabDEGsTFPubmed
#' @format A data frame with 10 rows and 8 variables
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

#' Genome Reference Consortium Human retrived from biomaRt library
#' The data set contains the following fields:
#' \itemize{
#' \item ensembl_gene_id
#' \item ensembl_transcript_id
#' \item chromosome_name
#' \item start_position
#' \item end_position
#' \item strand external_gene_name
#' \item external_transcript_name
#' \item external_gene_source
#' \item external_transcript_source_name
#' \item hgnc_id
#' \item entrezgene
#'}
#' @docType data
#' @keywords internal
#' @name gene.location
#' @format A data frame with 24474 rows and 12 variables
NULL
