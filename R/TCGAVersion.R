#' @title TCGA Version
#' @description  TCGA Version
#' @param Tumor  a character string indicating the cancer type for
#'        which to download data. Options include ACC, BLCA, BRCA,
#'        CESC, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML,
#'        LGG, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, SARC, SKCM, STAD,
#'        THCA, UCEC, UCS. Look at https://tcga-data.nci.nih.gov/tcga/
#'        for Available Cancer Types.
#' @param PlatformType illuminahiseq_rnaseq, agilentg4502a_07_3,
#'        illuminahiseq_rnaseqv2, humanmethylation27, humanmethylation450,
#'         illuminaga_mirnaseq, genome_wide_snp_6
# @param PlatformAndAssociatedData data frame 615 observations of 12 variables,
#        indicating the different characteristics of the data
#        e.g. tumour, type, species.
#' @importFrom rvest html html_text
#' @importFrom stringr str_split
#' @examples
#' TCGAVersion("LGG","illuminahiseq_rnaseqv2")
#' @export
#' @return Data frame with version, date, number of samples,size of
#'         the platform and tumor
TCGAVersion <- function(tumor = NULL, platform = NULL){

    if(is.null(tumor) && is.null(platform)){
        message("Please provide one tumor and platform")
    }
    query <- TCGAQuery(tumor,platform)
    root <- "https://tcga-data.nci.nih.gov/"
    path <- paste0(root,unique(dirname(query$deployLocation)))
    html <- html(path)
    text <- html_text(html)
    lines <- unlist(str_split(text,"\n"))
    folders <- lines[grep(".tar.gz ",lines)]
    ret <- data.frame(date = query$addedDate,
                         name=query$name,
                         samples = unlist(lapply(query$barcode,
                                function(x){length(unlist(strsplit(x,",")))}))
                      )
    for(i in seq_along(query$name)){
        idx <- grep(query[i,"name"],folders)
        ret[i,"hours"] <-  as.character(str_match(folders[idx],"[0-9]{2}:[0-9]{2}"))
        regex <- "[0-9]+\\.?[0-9]*([K]{1}|[M]{1}|[G]{1})"
        ret[i,"size"] <-  as.character(str_match(folders[idx],regex)[1,1])
    }
    message("==================  FOUND ==================")
    message("Platform: ", platform)
    message("Level 1 versions: ", length(grep("Level_1", query$name)))
    message("Level 2 versions: ", length(grep("Level_2", query$name)))
    message("Level 3 versions: ", length(grep("Level_3", query$name)))
    message("Mage versions: ", length(grep("mage-tab", query$name)))
    message("============================================")
    return(ret)
}
