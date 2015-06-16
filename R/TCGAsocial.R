#' @title TCGAsocial
#' @description Find number of downloads about a package in CRAN or BioC
#' @param siteToFind website related to TCGA or a package
#' @param listPackage list of package to find
#' @param KeyInfo is a key to find related to TGGA or package
# @import RCurl
#' @export
#' @return table with number of downloads about a package
TCGAsocial<- function(siteToFind, listPackage=NULL,KeyInfo=NULL){

    if( siteToFind == "bioconductor.org"){
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
    }


    if( siteToFind == "biostars.org"){
        siteBioStar <- "https://www.biostars.org/local/search/page/?page=1&sort=update&limit=all%50time&q="
        siteQuestions <- "https://www.biostars.org/p"
        siteBioStarKey <- paste(siteBioStar,KeyInfo,sep="")

        tmp <- .DownloadURL(siteBioStarKey)
        tmp2 <- tmp[ grep("<h4>",tmp)]

        TableQuestions <- matrix(0, length(tmp2), 3)
        rownames(TableQuestions)<-paste("q",c(1:length(tmp2)),sep="")
        colnames(TableQuestions) <- c("question","BiostarsSite","PackageSuggested")
        TableQuestions <- as.data.frame(TableQuestions)

        for ( i in 1:nrow(TableQuestions)){
            #print(i)
            questiontofind <- tmp2[i]
            questiontofind <- gsub("<h4>","", questiontofind)
            qst_find_site <-gsub("<a href=","", as.matrix(unlist(strsplit(questiontofind,">")))[1])
            qst_find_site2 <-gsub("<a href=","", as.matrix(unlist(strsplit(qst_find_site,"p")))[2])

            qst_find_site2_sub <- substr(qst_find_site2,1, 7)
            TableQuestions[i,"BiostarsSite"] <- qst_find_site2_sub
            newsite_tofind <- paste(siteQuestions,qst_find_site2_sub,sep="")

            tmpPack <- .DownloadURL(newsite_tofind)

            if( length(grep("package",tolower(tmpPack)))!=0){
                pos <- grep("package",tolower(tmpPack))
                if( length(pos)!=1){
                    pos <- pos[1]
                }
                tmpPackage <-  tmpPack[pos]
                PackageSuggested <- tmpPackage
                #  print(PackageSuggested)
                TableQuestions[i,"PackageSuggested"] <- substr(PackageSuggested,1, 64)
            }

            tmp3a <- gsub("</a","", as.matrix(unlist(strsplit(questiontofind,">")))[2])
            tmp3a <- gsub("&#39;", ".", tmp3a)
            TableQuestions[i,"question"] <- tmp3a
        }
        site3 <- "http://www.bioconductor.org/packages/3.1/bioc/"
        tmp <- .DownloadURL(site3)
        tmpPack <-  tmp[grep(".html",tmp)]
        tmpPackMatrix <- as.matrix(tmpPack)
        posA <- grep("bioconductor",tolower(tmpPackMatrix))[1]+1
        posB <- grep("bioconductor",tolower(tmpPackMatrix))[2]-1
        tmpPackMatrixNew <- tmpPackMatrix[posA:posB,]

        TableQuestions <- TableQuestions[order(TableQuestions$PackageSuggested,decreasing=TRUE),]
        TablePackage <- TableQuestions
    }



    return(TablePackage)

}

