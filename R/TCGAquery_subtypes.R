#' @title Retrieve multiple molecular subtypes according last publications
#' @description Retrieve multiple molecular subtypes according last publications
#' @param tumor is character such as cancer from TCGA. Examples:
#' \tabular{lllll}{
#'OV   \tab BRCA \tab CESC \tab ESCA \tab PCPG\cr
#'LUSC \tab LGG  \tab SKCM \tab KICH \tab CHOL\cr
#'GBM  \tab UCEC \tab PRAD \tab PAAD \tab THYM\cr
#'KIRC \tab THCA \tab SARC \tab LAML \tab TGCT\cr
#'COAD \tab KIRP \tab HNSC \tab ACC  \tab UVM \cr
#'READ \tab BLCA \tab DLBC \tab UCS  \tab FPPP\cr
#'LUAD \tab LIHC \tab STAD \tab MESO \tab CNTL
#'}
#' @param path Directory to save the downloaded data
#' @importFrom rvest html
#' @importFrom stringr str_match
#' @importFrom xlsx read.xlsx2
#' @examples
#' GBM_subtypes <- TCGAquery_subtypes(tumor = "gbm",path ="dataGBM")
#' @export
#' @return data.frame with information about molecular cancer subtypes
TCGAquery_subtypes <- function(tumor = NULL, path = ".") {

    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    root <- "https://tcga-data.nci.nih.gov"

    pg <- html("https://tcga-data.nci.nih.gov/docs/publications/")
    pg <- pg %>% html_nodes("a") %>% html_attr("href")
    pg <- pg[grep("https://tcga-data.nci.nih.gov/docs/publications/",pg)]
    year <- str_match(pg, "[0-9]{4}")
    cancer <- gsub("_","",str_match(pg, "[a-z]{3,4}\\_"))

    df <- data.frame(path=as.character(pg),
                     year=as.numeric(year),
                     cancer=as.character(cancer))
    df <- df[order(df$year,decreasing=T),]
    df <- df[!duplicated(df$cancer),]
    rownames(df) <- tolower(df$cancer)

    site2 <- as.character(df[tolower(tumor),"path"])

    tmp <- .DownloadURL(site2)

    tmp2 <- tmp[grep("xls",tolower(tmp))]

    if(length(grep("subtype", tmp)) > 1){
        tmp2 <- tmp[grep("subtype",tolower(tmp))]
    }

    tmp3 <- tmp2[grep("xls",tolower(tmp2))]

    if( length(tmp3)!=1){
        tmp3 <- tmp3[1]
    }

    tmp3 <- gsub("<li>","",tmp3)
    tmp4 <- unlist(strsplit(tmp3,">"))[1]
    tmp5 <- unlist(strsplit(tmp4,"<a href="))[2]
    tmp6 <- unlist(strsplit(tmp5,"xlsx"))[1]

    Filelocation <-  unlist(strsplit(tmp6,df[tolower(tumor),"year"]))[2]
    if( !is.na(Filelocation)){
        Filelocation <- substr(Filelocation,2,nchar(Filelocation))
    }

    if( is.na(Filelocation)){
        Filelocation <- substr(tmp6,2,nchar(tmp6))
    }


    FileSubtypes <- paste0(df[tolower(tumor),"path"],Filelocation,"xlsx")

    file <- paste0(path, "/", gsub("/","_",Filelocation) ,"xlsx")

    if(is.windows()){
        suppressWarnings(
            downloader::download(FileSubtypes,file, quiet = TRUE, method = "wininet")
        )
    } else {
        suppressWarnings(
            downloader::download(FileSubtypes,file, quiet = TRUE)
        )
    }

    table <- read.xlsx2(file,1,stringsAsFactors = NULL)

    return(table)
}


