#' @title Shows a summary (version, date, number of samples, size of the data) of
#' all versions of data for a given tumor and platform.
#' @description  As every new version of data has diferent samples and size
#' the TCGAquery_Version function shows a summary (version, date, number of samples,
#' size of the data) of all versions of data for a given tumor and platform.
#'
#' @param tumor  a character string indicating the cancer type for
#'        which to download data. Options include ACC, BLCA, BRCA,
#'        CESC, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML,
#'        LGG, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, SARC, SKCM, STAD,
#'        THCA, UCEC, UCS. Look at https://tcga-data.nci.nih.gov/tcga/
#'        for Available Cancer Types.
#' @param platform illuminahiseq_rnaseq, agilentg4502a_07_3,
#'        illuminahiseq_rnaseqv2, humanmethylation27, humanmethylation450,
#'         illuminaga_mirnaseq, genome_wide_snp_6
#' @importFrom rvest html html_text
#' @importFrom stringr str_split str_trim
#' @examples
#' TCGAquery_Version('GBM','illuminahiseq_rnaseqv2')
#' @export
#' @return Data frame with version, date, number of samples,size of
#'         the platform and tumor
TCGAquery_Version <- function(tumor = NULL, platform = NULL) {

    if (is.null(tumor) && is.null(platform)) {
        message("Please provide one tumor and platform")
    }

    # Get last version of files
    queryLastVersion <- TCGAquery(tumor = c(tumor),
                                  platform = c(platform),
                                  level = 3)
    last_version <- as.numeric(str_sub(queryLastVersion$name[1],-3,-3))

    BarcodeList <- vector("list",last_version)
    ret <-  queryLastVersion[1,c(1,2,7)]
    ret$name <- last_version
    listSamplesVer <- TCGAquery_samplesfilter(queryLastVersion)
    ret$Samples <- length(unlist(listSamplesVer))
    ret$BarcodeList <- listSamplesVer

    ret <- as.data.frame(ret)
    colnames(ret) <- c("Date","BaseName","Version","Samples","BarcodeList")

    last_version <- last_version - 1
    for (idx in last_version:0) {
        queryVers <- TCGAquery(tumor = c(tumor),
                               platform = c(platform),
                               level = 3,
                               version = list(c(platform,tumor,idx)))
        listSamplesVer <- TCGAquery_samplesfilter(queryVers)
        aux <- data.frame(as.character(queryVers[1,1]),
                          queryVers[1,2],
                 idx,
                 length(unlist(listSamplesVer)))
        aux$BarcodeList <- listSamplesVer
        list(listSamplesVer)
        colnames(aux) <- c("Date","BaseName","Version","Samples","BarcodeList")
        ret <- rbind(ret,aux)
    }

    colnames(ret) <- c("Date","BaseName","Version","Samples","BarcodeList")

    return(ret)
}

getTotalSize <- function(sizeList) {
    sizeK <- sizeList[grep("K", sizeList)]
    sizeM <- sizeList[grep("M", sizeList)]
    sizeG <- sizeList[grep("G", sizeList)]
    totalK <- round(sum(as.numeric(gsub("K", "", sizeK)))/1000)
    totalM <- sum(as.numeric(gsub("M", "", sizeM)))
    totalG <- sum(as.numeric(gsub("G", "", sizeG))) * 1000

    return(totalM + totalK + totalG)
}
