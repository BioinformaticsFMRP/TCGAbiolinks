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
#' TCGAquery_Version('LGG','agilentg4502a_07_3')
#' @export
#' @return Data frame with version, date, number of samples,size of
#'         the platform and tumor
TCGAquery_Version <- function(tumor = NULL, platform = NULL) {

    if (is.null(tumor) && is.null(platform)) {
        message("Please provide one tumor and platform")
    }
    query <- TCGAquery(tumor, platform)
    root <- "https://tcga-data.nci.nih.gov/"
    path <- paste0(root, unique(dirname(query$deployLocation)))
    html <- html(path)
    text <- html_text(html)
    lines <- unlist(str_split(text, "\n"))
    folders <- lines[grep("./", lines)]
    folders <- folders[-grep("Index of", folders)]
    folders <- folders[-grep("mage-tab", folders)]
    ret <- data.frame(Version = NULL, Date = NULL,
                      Samples = NULL, SizeMB = NULL)
    message(paste0("Found ", length(folders), " Version of ", platform))
    for (i in seq_along(folders)) {
        regex <- "^.*/"
        Version <- as.character(str_trim(str_match(folders[i], regex)))
        ret[i, "Version"] <- Version
        message(Version)
        regex <- "[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}"
        ret[i, "Date"] <- as.character(str_match(folders[i], regex))
        subdir <- html_text(html(file.path(path, Version)))
        x <- unlist(str_split(subdir, "\n"))

        if (platform == "illuminahiseq_rnaseq") {
            x <- x[grep("gene.quantification", x)]
        }
        if (platform == "agilentg4502a_07_3") {
            x <- x[grep("tcga_level3", x)]
        }
        if (platform == "illuminahiseq_rnaseqv2") {
            x <- x[grep("rsem.genes.results", x)]
        }
        if (platform == "humanmethylation27") {
            x <- x[grep("HumanMethylation27", x)]
        }
        if (platform == "humanmethylation450") {
            x <- x[grep("HumanMethylation450", x)]
        }
        if (platform == "illuminaga_mirnaseq") {
            x <- x[grep("mirna.quantification", x)]
        }
        if (platform == "genome_wide_snp_6") {
            x <- x[grep("hg19.seg", x)]
        }

        regex <- "[0-9]+\\.?[0-9]*([K]{1}|[M]{1}|[G]{1})"
        size <- as.character(str_match(x, regex)[, 1])
        ret[i, "Samples"] <- length(size)
        size <- getTotalSize(size)
        ret[i, "SizeMB"] <- size
    }
    message("==================  FOUND ==================")
    message("Platform: ", platform)
    message("Level 1 versions: ", length(grep("Level_1", ret$Version)))
    message("Level 2 versions: ", length(grep("Level_2", ret$Version)))
    message("Level 3 versions: ", length(grep("Level_3", ret$Version)))
    message("Mage versions: ", length(grep("mage-tab", ret$Version)))
    message("============================================")
    ret <- ret [ order( substr(gsub("-","",ret$Date),1,8),decreasing =F ), ]
    rownames(ret) <- NULL
    BarcodeList <- vector("list",nrow(ret))
    names(BarcodeList) <- ret$Version

    for( idx in 1: nrow(ret)){
        queryVers <- TCGAquery(tumor = c(tumor),
                               platform = c(platform),
                               level = 3,
                               version = list(c(platform,tumor,idx-1)))
        listSamplesVer <- TCGAquery_samplesfilter(queryVers)
        BarcodeList[[idx]] <- as.character(unlist(listSamplesVer))
        ret[idx,"Samples"] <- length(unlist(listSamplesVer))
    }

    ans <- list(TableVersion = ret, BarcodeList = BarcodeList)

    return(ans)
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
