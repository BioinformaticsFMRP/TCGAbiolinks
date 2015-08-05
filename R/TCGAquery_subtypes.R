#' @title Retrieve multiple molecular subtypes according last publications
#' @description Retrieve multiple molecular subtypes according last publications
#' @param tumor is character such as cancer from TCGA. Examples:
#' \tabular{lllll}{
#' blca \tab gbm  \tab hnsc \tab kich \tab kirc\cr
#' laml \tab lgg  \tab luad \tab skcm \tab stad\cr
#' thca \tab ucec
#'}
#' @param path Directory to save the downloaded data
#' @importFrom downloader download
#' @importFrom rvest html html_nodes html_attr
#' @importFrom stringr str_match
#' @importFrom xlsx read.xlsx2
#' @examples
#' \dontrun{
#' GBM_path_subtypes <- TCGAquery_subtypes(tumor = "gbm",path ="dataGBM")
#' # to prepare run:
#' # require(xlsx)
#' # GBM_subtypes <- read.xlsx2(GBM_path_subtypes,1,stringsAsFactors = NULL)
#' }
#' @export
#' @return data.frame with information about molecular cancer subtypes
TCGAquery_subtypes <- function(tumor = NULL, path = ".") {

    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    root <- "https://tcga-data.nci.nih.gov"
    pg <- html("https://tcga-data.nci.nih.gov/docs/publications/")

    pg_nodes <- html_nodes(x = pg, css = "a")
    pg <- as.character(html_attr(x = pg_nodes, name = "href"))

    #pg <- pg %>% html_nodes("a") %>% html_attr("href")
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
    pg2 <- html(site2)
    pg_nodes2 <- html_nodes(x = pg2, css = "a")
    pg2 <- as.character(html_attr(x = pg_nodes2, name = "href"))

#   pg2 <- pg2 magrittr::%>% html_nodes("a") magrittr::%>% html_attr("href")
    pg3 <- pg2[grep("xls",pg2)]

    pgyear <- as.character(df[tolower(tumor),"year"])

    if( length(pg3)!=1){ pg4 <- pg3[1] }
    if( length(pg3)==1){ pg4 <- pg3}

    FileSubtypes <- paste0(df[tolower(tumor),"path"],pg4)
    filetoDown <- paste0(path, "/", gsub("/","_",pg4))

    if( length(grep("subtype",pg3))==1){
        pg4 <- pg3[grep("subtype",pg3)]
        FileSubtypes <- pg4
        pg5 <- unlist(strsplit(pg4, pgyear))[2]
        filetoDown <- paste0(path, "/", gsub("/","_",substr(pg5,2,nchar(pg5))))
    }

    suppressWarnings(
        download(FileSubtypes,filetoDown, quiet = TRUE,  mode="wb")
    )
    #table <- xlsx::read.xlsx2(file,1,stringsAsFactors = NULL)
    return(filetoDown)
}
