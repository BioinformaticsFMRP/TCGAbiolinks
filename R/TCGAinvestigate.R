#' @title TCGAinvestigate
#' @description Find most studied TF in pubmed related to a specific cancer,
#'  disease, or tissue
#' @param tumor tumor
#' @param dataDEGsFiltLevelTF dataDEGsFiltLevelTF
#' @param topgenes topgenes
#' @importFrom RCurl url.exists curlVersion
#' @export
#' @return table with number of pubmed related to tfs.
TCGAinvestigate<- function(tumor,dataDEGsFiltLevelTF,topgenes){
    site <- "http://www.ncbi.nlm.nih.gov/pubmed/?term="

    # GenesTofix <- c("JUN","HR","HOMEZ",
    #                 "ANKAR", "REST", "BATF", "MAX", "ECD", "FOS")
    # dataDEGsFiltLevelTF <- dataDEGsFiltLevelTF[
    #                        setdiff(dataDEGsFiltLevelTFs$mRNA,GenesTofix),]

    dataDEGsFiltLevelTF <- dataDEGsFiltLevelTF[1:topgenes,]
    Pubmed <- matrix(0, nrow(dataDEGsFiltLevelTF), 1)
    PMID <- matrix(0, nrow(dataDEGsFiltLevelTF), 1)


    dataDEGsFiltLevelTF <- cbind(dataDEGsFiltLevelTF,Pubmed,PMID)
    dataDEGsFiltLevelTF<-as.data.frame(dataDEGsFiltLevelTF)

    for (k in 1:nrow( dataDEGsFiltLevelTF)){
        CurrentGene <- dataDEGsFiltLevelTF$mRNA[k]
        site2 <- paste(site,CurrentGene, "+", tumor,sep="")

        if(interactive() && ("ssl" %in% names(curlVersion()$features)) &&
           url.exists(site2)) {
            x = tryCatch(RCurl::getURL(site2), error = function(e) {
                RCurl::getURL(site2, ssl.verifypeer = FALSE) })
        }



        if ( length(grep("No items found.",x))!=1){

            if (length(grep("Display Settings",x))==1){
                x6 <- 1
                dataDEGsFiltLevelTF[k,"PMID"] <- substr(gsub("</dt> <dd>","",
                                        unlist(strsplit(x,"PMID:"))[2]),1,8)

            }


            if (length(grep("result_count",x))==1){
                x2a <- unlist(strsplit(x,"result_count"))[2]

                tmpPMID2 <- unlist(strsplit(x2a,"UidCheckBox"))
                tmpPMID3 <- tmpPMID2[grep("<span>",tmpPMID2)]
                dataDEGsFiltLevelTF[k,"PMID"] <- as.character(paste(
                    substr(tmpPMID3,1,8),collapse="; "))


                x3a <-  unlist(strsplit(x2a,"</h2>"))[1]

                if( length(grep("of",x3a))!=1){
                    x6 <- as.numeric(unlist(strsplit(x3a,": "))[2])
                } else { x6 <- as.numeric(unlist(strsplit(x3a,"of "))[2]) }
            }


            if (length(grep("following term was not found",x)) == 1) { x6 <- 0 }

            if (length(grep("Search instead for",x)) == 1) { x6 <- 1 }

            if (CurrentGene == "JUN")   { x6 <- 0 }
            if (CurrentGene == "HR")    { x6 <- 0 }
            if (CurrentGene == "HOMEZ") { x6 <- 0 }
            if (CurrentGene == "ANKAR") { x6 <- 0 }
            if (CurrentGene == "REST")  { x6 <- 0 }
            if (CurrentGene == "BATF")  { x6 <- 0 }
            if (CurrentGene == "MAX")   { x6 <- 0 }
            # if (CurrentGene =="FOS"){       x6 <- 0     }
            if (CurrentGene == "ECD")   { x6 <- 0 }
            # HR, HOMEZ, ANKAR, REST

            dataDEGsFiltLevelTF[k,"Pubmed"] <- x6
            print(paste("Cancer ", tumor, "with TF n. ",k, "of " ,
                        nrow( dataDEGsFiltLevelTF)," : ", CurrentGene,
                        "found n. ", x6, "pubmed."))

        }

        else{
            print(paste("Cancer ", tumor, "with TF n. ",k, "of " ,
                        nrow( dataDEGsFiltLevelTF)," : ", CurrentGene,
                        "no item found in pubmed."))
            dataDEGsFiltLevelTF[k,"Pubmed"] <- 0
        }

    }

    dataDEGsFiltLevelTF <- dataDEGsFiltLevelTF[order(dataDEGsFiltLevelTF$Pubmed,
                                                     decreasing = TRUE),]

    if( sum(dataDEGsFiltLevelTF$Pubmed == 1) != 0) {
        dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 1,][ which( nchar(dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 1,]$PMID) > 8),"PMID"] <- substr(dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 1,][ which( nchar(dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 1,]$PMID) > 8),"PMID"],1,8)
    }
    if( sum(dataDEGsFiltLevelTF$Pubmed == 0) != 0){
        dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 0,"PMID"] <- 0
    }

    tabDEGsTFPubmed$Tumor <- round(tabDEGsTFPubmed$Tumor)
    tabDEGsTFPubmed$Normal <- round(tabDEGsTFPubmed$Normal)
    tabDEGsTFPubmed$Delta <- round(tabDEGsTFPubmed$Delta)


    return(dataDEGsFiltLevelTF)
}
