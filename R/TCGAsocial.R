#' @title TCGAsocial
#' @description Find number of downloads about a package in CRAN or BioC
#' @param listPackage list of package to find
# @import RCurl
#' @export
#' @return table with number of downloads about a package
TCGAsocial<- function(listPackage){

    siteBioC <- "http://www.bioconductor.org/packages/stats/index.html"

    TablePackage <- matrix(0, length(listPackage), 2)
    rownames(TablePackage)<-listPackage
    colnames(TablePackage) <- c("Package","NumberDownload")
    TablePackage <- as.data.frame(TablePackage)
    TablePackage$Package <- listPackage

    tmp <- .DownloadURL(siteBioC)

    for ( i in 1:nrow(TablePackage)){
        packagetofind <- listPackage[i]
        pos <- grep(tolower(packagetofind), tolower(tmp ))
        pos1 <- pos[1]
        tmp3 <- tmp[pos1]
        tmp3a <- gsub("</A></TD></TR>","", as.matrix(unlist(strsplit(tmp3,"&nbsp;")))[2])
        tmp4 <- as.numeric(substr(as.character(tmp3a),2,nchar(tmp3a)-1))
        TablePackage[i,"NumberDownload"] <- tmp4
    }

    TablePackage <- TablePackage[order(TablePackage$NumberDownload,decreasing=TRUE),]

    return(TablePackage)

}

