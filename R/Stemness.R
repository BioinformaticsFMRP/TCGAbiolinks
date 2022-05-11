#' @title Generate Stemness Score based on RNASeq (mRNAsi stemness index) Malta et al., Cell, 2018
#' @description TCGAanalyze_Stemness generate the mRNAsi score
#' @param stemSig is a vector of the stemness Signature generated using gelnet package.
#' Please check the data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5902191/
#' \itemize{
#' \item \link[TCGAbiolinks]{SC_PCBC_stemSig} - Stemness Score
#' \item \link[TCGAbiolinks]{DE_PCBC_stemSig} - endoderm score
#' \item \link[TCGAbiolinks]{EB_PCBC_stemSig} -  embryoid bodies score
#' \item \link[TCGAbiolinks]{ECTO_PCBC_stemSig} - ectoderm score
#' \item \link[TCGAbiolinks]{MESO_PCBC_stemSig} - mesoderm score
#' }
#' @param dataGE is a matrix of Gene expression (genes in rows, samples in cols) from TCGAprepare
#' @param colname.score Column name of the output. Default "stemness_score"
#' @export
#' @return table with samples and selected score
#' @examples
#'  # Selecting TCGA breast cancer (10 samples) for example stored in dataBRCA
#'  dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)
#'
#'  # quantile filter of genes
#'  dataFilt <- TCGAanalyze_Filtering(
#'    tabDF = dataNorm,
#'    method = "quantile",
#'    qnt.cut =  0.25
#'  )
#'  Stemness_score <- TCGAanalyze_Stemness(
#'      stemSig = SC_PCBC_stemSig,
#'      dataGE = dataFilt,
#'      colname.score = "SC_PCBC_stem_score"
#'   )
#'  ECTO_score <- TCGAanalyze_Stemness(
#'      stemSig = ECTO_PCBC_stemSig,
#'      dataGE = dataFilt,
#'      colname.score = "ECTO_PCBC_stem_score"
#'   )
#'   MESO_score <- TCGAanalyze_Stemness(
#'      stemSig = MESO_PCBC_stemSig,
#'      dataGE = dataFilt,
#'      colname.score = "MESO_PCBC_stem_score"
#'   )
TCGAanalyze_Stemness <- function(
    stemSig,
    dataGE,
    colname.score = "stemness_score"
) {
    reads <- dataGE
    X <- reads
    w <- stemSig
    commonStemsigGenes <- intersect(names(w), rownames(X))

    X <- X[commonStemsigGenes, ]
    w <- w[rownames(X)]

    # Score the Matrix X using Spearman correlation.

    s <- apply(X, 2, function(z) {
            cor(z, w, method = "sp", use = "complete.obs")
        })

    ## Scale the scores to be between 0 and 1
    s <- s - min(s)
    s <- s / max(s)

    dataSce_stemness <- cbind(s)

    score <- data.frame("Sample" = colnames(reads))
    rownames(score) <- colnames(reads)

    score[[colname.score]] <- 0
    score[rownames(dataSce_stemness), colname.score] <- as.numeric(dataSce_stemness)

    return(score)
}
