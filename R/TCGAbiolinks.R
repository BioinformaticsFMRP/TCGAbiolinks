#' Clinical data TCGA BRCA
#' @docType data
#' @keywords internal
#' @name clinBRCA
#' @format A data frame with 1061 rows and 109 variables
NULL

#' tabSurvKMcompleteDEGs
#' @docType data
#' @keywords internal
#' @name tabSurvKMcompleteDEGs
#' @format A data frame with 200 rows and 7 variables
NULL



#' TCGA data matrix BRCA
#' @docType data
#' @keywords internal
#' @name dataBRCA
#' @format A data frame with 20531 rows (genes) and 50 variables (samples)
NULL

#' TCGA data SummarizedExperiment READ
#' @docType data
#' @keywords internal
#' @name dataREAD
#' @format A SummarizedExperiment of READ with 2 samples
NULL

#' TCGA data matrix READ
#' @docType data
#' @keywords internal
#' @name dataREAD_df
#' @format A data frame with 20531 rows (genes) and 2 variables (samples)
NULL

#' TCGA data matrix BRCA DEGs
#' @docType data
#' @keywords internal
#' @name dataDEGsFiltLevel
#' @format A data frame with 3649 rows and 6 variables
NULL

#' geneInfo for normalization of RNAseq data
#' @docType data
#' @keywords internal
#' @name geneInfo
#' @format A data frame with 20531 rows and 2 variables
NULL

#' geneInfoHT for normalization of HTseq data
#' @description
#' Code to generate the data
#' ```{R, eval = F}
#' library(EDASeq)
#' library(biomaRt)
#' #get ensembl gene IDs for hg38
#' ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#' biomart_getID <- getBM(attributes = c("ensembl_gene_id"), mart = ensembl)
#' #get gene length and GC content for all IDs
#'
#' step <- 500
#' geneInfoHT <- plyr::adply(seq(1,length(biomart_getID$ensembl_gene_id),step),.margins = 1,.fun = function(x){
#'     begin <- x
#'    end <- x + step
#'     if(end > length(biomart_getID$ensembl_gene_id)) end <- length(biomart_getID$ensembl_gene_id)
#'     file <- paste0("geneInfoHT_from_",begin,"_to_",end,".rda")
#'     if(!file.exists(file)){
#'         df <- getGeneLengthAndGCContent(biomart_getID$ensembl_gene_id[begin:end] , org="hsa", mode = c("biomart"))
#'         save(df,file = file)
#'     } else {
#'         df <- get(load(file))
#'     }
#'     df
#' },.progress = "time")
#' saveRDS(getdata, file = "getGLGC_download.RDS")a
#' save(getdata, file = "getGLGC_download.rda")
#' #Save output as data frame with correct header names
#' geneInfoHT <- data.frame(
#'     geneLength = getdata[,1] ,
#'     gcContent = getdata[,2]
#' )
#' #Save final table
#' save(geneInfoHT, file = "data/geneInfoHT.rda")
#' ```
#' @docType data
#' @keywords internal
#' @name geneInfoHT
#' @format A data frame with 23486 rows and 2 variables
NULL

#' TCGA batch information from Biospecimen Metadata Browser
#' @docType data
#' @keywords internal
#' @name batch.info
#' @format A data frame with 11382 rows and 3 variables
NULL

#' TCGA CHOL MAF transformed to maftools object
#' @docType data
#' @keywords internal
#' @name chol_maf
#' @format An object of class  MAF
NULL

#' TCGA CHOL MAF
#' @docType data
#' @keywords internal
#' @name bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.maf
#' @format A tibble: 3,555 x 34
NULL

#' MSI data for two samples
#' @docType data
#' @keywords internal
#' @name msi_results
#' @format A data frame: 2 rows, 4 columns
NULL

#' A RangedSummarizedExperiment two samples with gene expression data from vignette
#' aligned against hg38
#' @docType data
#' @keywords internal
#' @name gbm.exp.harmonized
#' @format A RangedSummarizedExperiment: 56963 genes, 2 samples
NULL

#' A RangedSummarizedExperiment two samples with gene expression data from vignette
#' aligned against hg19
#' @docType data
#' @keywords internal
#' @name gbm.exp.legacy
#' @format A RangedSummarizedExperiment: 21022 genes, 2 samples
NULL

#' A DNA methylation RangedSummarizedExperiment for 8 samples (only first 20 probes)
#' aligned against hg19
#' @docType data
#' @keywords internal
#' @name met.gbm.27k
#' @format A RangedSummarizedExperiment: 20 probes, 8 samples
NULL

#' A list of data frames with clinical data parsed from XML (code in vignettes)
#' @docType data
#' @keywords internal
#' @name clinical.biotab
#' @format A list with 7 elements
NULL

#' A data frame with all TCGA molecular subtypes
#' @docType data
#' @keywords internal
#' @name pancan2018
#' @format A data frame with 7,734 lines and 10 columns
NULL

#' A numeric vector with stem cell-like signature trained on PCBC's dataset
#' @docType data
#' @keywords internal
#' @name SC_PCBC_stemSig
#' @format A numeric vector with 12956 genes
NULL

#' A numeric vector with SC-derived mesoderm (MESO) signature trained on PCBC's dataset
#' @docType data
#' @keywords internal
#' @name MESO_PCBC_stemSig
#' @format A numeric vector with 12956 genes
NULL


#' A numeric vector with SC-derived ectoderm (ECTO) signature trained on PCBC's dataset
#' @docType data
#' @keywords internal
#' @name ECTO_PCBC_stemSig
#' @format A numeric vector with 12956 genes
NULL


#' A numeric vector with SC-derived definitive endoderm (DE) signature trained on PCBC's dataset
#' @docType data
#' @keywords internal
#' @name DE_PCBC_stemSig
#' @format A numeric vector with 12956 genes
NULL

#' A numeric vector with stem cell (SC)-derived embryoid bodies (EB) signature trained on PCBC's dataset
#' @docType data
#' @keywords internal
#' @name EB_PCBC_stemSig
#' @format A numeric vector with 12956 genes
NULL

#' Result of gliomaclassifier function
#' @docType data
#' @keywords internal
#' @name classification
#' @format A list of data frames
NULL
