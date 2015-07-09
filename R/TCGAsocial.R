#' @title TCGAsocial
#' @description Find number of downloads about a package in CRAN or BioC
#' @param siteToFind website ("bioconductor.org", "biostars.org", "stackoverflow)
#' related to TCGA or a package
#' @param listPackage list of package to find in bioconductor
#' @param KeyInfo is a key to find in biostars related to TGGA or package
#' @export
#' @return table with number of downloads about a package
#' @examples
#' TCGAsocial("bioconductor.org","BiocCheck")
TCGAsocial <- function(siteToFind=NULL, listPackage=NULL,KeyInfo=NULL){


    # Find all packages in bioconductor
    site3 <- "http://www.bioconductor.org/packages/3.1/bioc/"
    tmp <- .DownloadURL(site3)
    tmpPack <-  tmp[grep(".html",tmp)]
    tmpPackMatrix <- as.matrix(tmpPack)
    posA <- grep("bioconductor",tolower(tmpPackMatrix))[1]+1
    posB <- grep("bioconductor",tolower(tmpPackMatrix))[2]-1
    tmpPackMatrixNew <- tmpPackMatrix[posA:posB,]
    tmpPackMatrixNew <- gsub("<td><a","",tmpPackMatrixNew)
    tmpPck2 <- lapply(tmpPackMatrixNew, function(x) gsub("</a></td>","", as.matrix(unlist(strsplit(x,">")))[2]))
    tmpPck3 <- as.matrix(unlist(tmpPck2))
    tmpPck2a <- lapply(tmpPck2, function(x) gsub("</a","", x))
    BiocPackageList <- as.matrix(unlist(tmpPck2a))

    if( siteToFind == "bioconductor.org"){

        if(is.null(listPackage)) {
            msg  <- paste0(
                "\nPlease, provide a listofPackage argument\n",
                "Example:  TCGAsocial('bioconductor.org','BiocCheck')")
            stop(msg)
        }

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

        if(is.null(KeyInfo)) {
            msg  <- paste0(
                "\nPlease, provide a KeyInfo argument\n",
                "Example: TCGAsocial('biostars.org', KeyInfo='methylation')")
            stop(msg)
        }


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

                PackMat <- sapply(BiocPackageList, grepl, tmpPackage, ignore.case=TRUE)
                if(sum(PackMat)>=1){
                    print(which(PackMat == TRUE))
                    PackageSuggested <- paste(names(PackMat[which(PackMat == TRUE)]),collapse=";")
                    #TableQuestions[i,"PackageSuggested"] <- PackageSuggested
                    TableQuestions[i,"PackageSuggested"] <- substr(PackageSuggested,1, 64)
                }
                # print(PackageSuggested)
            }

            tmp3a <- gsub("</a","", as.matrix(unlist(strsplit(questiontofind,">")))[2])
            tmp3a <- gsub("&#39;", ".", tmp3a)
            TableQuestions[i,"question"] <- tmp3a
        }

        TableQuestions <- TableQuestions[order(TableQuestions$PackageSuggested,decreasing=TRUE),]
        TablePackage <- TableQuestions
    }

    if( siteToFind == "support.bioconductor.org"){

        if(is.null(KeyInfo)) {
            msg  <- paste0(
                "\nPlease, provide a KeyInfo argument\n",
                "Example: TCGAsocial('support.bioconductor.org', KeyInfo='TCGA')")
            stop(msg)
        }

        sitesupportBioc_part1 <- "https://support.bioconductor.org/local/search/page/?page="
        sitesupportBioc_part2 <-"&sort=New%20answers&limit=All%20time&q="
        TablePackage <-NULL


        for( pg in 1:2){
            message(paste("pag",pg),sep="")
            #sitesupportBioc <- "https://support.bioconductor.org/local/search/page/?q="
            sitesupportBioc<- paste(sitesupportBioc_part1,pg,sitesupportBioc_part2,sep="")
            siteQuestions <- "https://support.bioconductor.org/p"
            sitesupportBiocKey <- paste(sitesupportBioc,tolower(KeyInfo),sep="")

            tmp <- .DownloadURL(sitesupportBiocKey)
            tmp2 <- tmp[ grep("<h4>",tmp)]

            TableQuestions <- matrix(0, length(tmp2), 3)
            rownames(TableQuestions)<-paste("q",c(1:length(tmp2)),sep="")
            colnames(TableQuestions) <- c("question","supportBiocsSite","PackageSuggested")
            TableQuestions <- as.data.frame(TableQuestions)

            for ( tbi in 1:nrow(TableQuestions)){

                questiontofind <- tmp2[tbi]

                if(length(grep("MISSING",questiontofind))!=1){

                questiontofind <- gsub("<h4>","", questiontofind)
                qst_find_site <-gsub("<a href=","", as.matrix(unlist(strsplit(questiontofind,">")))[1])
                qst_find_site2 <-gsub("<a href=","", as.matrix(unlist(strsplit(qst_find_site,"p")))[2])

                qst_find_site2_sub <- substr(qst_find_site2,1, 7)
                TableQuestions[tbi,"supportBiocsSite"] <- qst_find_site2_sub
                newsite_tofind <- paste(siteQuestions,qst_find_site2_sub,sep="")

                tmpPack <- .DownloadURL(newsite_tofind)
                # starting from Question
                # and removing similar post inside webpage
                tmpPack<-tmpPack[grep("Question",tmpPack)[3]:grep("Similar",tmpPack)]


                if( length(grep("package",tolower(tmpPack)))!=0){
                    pos <- grep("package",tolower(tmpPack))

                    if( length(pos)==1){
                        pos <- pos[1]
                        tmpPackage <-  tmpPack[pos]
                        PackMat <- sapply(BiocPackageList, grepl, tmpPackage, ignore.case=TRUE)
                    }

                    if( length(pos)!=1){
                        PackMatNew <- NULL
                        for( ip in 1: length(pos)){
                            tmpPackage <-  tmpPack[pos][ip]
                            PackMat <- sapply(BiocPackageList, grepl, tmpPackage, ignore.case=TRUE)
                            PackMatNew <- c(PackMatNew,PackMat)
                        }
                        PackMat<-PackMatNew
                    }

                    if(sum(PackMat)>=1){
                        print(which(PackMat == TRUE))
                        PackageSuggested <- paste(names(PackMat[which(PackMat == TRUE)]),collapse=";")
                        #TableQuestions[tbi,"PackageSuggested"] <- PackageSuggested
                        # print(PackageSuggested)
                        TableQuestions[tbi,"PackageSuggested"] <- substr(PackageSuggested,1, 64)

                    }
                }

                tmp3a <- gsub("</a","", as.matrix(unlist(strsplit(questiontofind,">")))[2])
                tmp3a <- gsub("&#39;", ".", tmp3a)
                TableQuestions[tbi,"question"] <- tmp3a
            }

            TableQuestions <- TableQuestions[order(TableQuestions$PackageSuggested,decreasing=TRUE),]
            TablePackage <- rbind(TablePackage,TableQuestions)
        }
        }
        TablePackage <- TablePackage[order(TablePackage$PackageSuggested,decreasing=TRUE),]
        TablePackage <- TablePackage[!duplicated(TablePackage$supportBiocsSite),]
    }

    return(TablePackage)

}

