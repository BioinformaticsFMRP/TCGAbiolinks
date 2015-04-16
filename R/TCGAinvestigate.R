#' @title TCGA investigate
#'
#' @description  Input (gene, cancer) return number of pubmed related to a specific cancer, disease, or tissue where the gene was studied.
#' @param tumor tumor code between "acc"  "blca" "brca" "cesc" "chol" "cntl" "coad" "dlbc" "esca" "fppp" "gbm"
#'                                 "hnsc" "kich" "kirc" "kirp" "laml" "lcml" "lgg"  "lihc" "lnnh" "luad" "lusc"
#'                                 "meso" "misc" "ov"   "paad" "pcpg" "prad" "read" "sarc" "skcm" "stad" "tgct"
#'                                 "thca" "thym" "ucec" "ucs"  "uvm", or tissues such as "breast", "liver", "colon", etc.
#'
#' @param CANCER_diff_normUQ_level_tf table out from DEA
#' @param topgenes number of genes of which investigate on Pubmed
#'
#' @examples
#' \dontrun{
#' TCGAinvestigate(tumor=BRCA,CANCER_diff_normUQ_level_tf,topgenes=100)
#' }
#'
#' @author Colaprico et al
#' @seealso TCGADownload
#' @export
#' @import Rcurl


TCGAinvestigate<- function(tumor,CANCER_diff_normUQ_level_tf,topgenes){
  site <- "http://www.ncbi.nlm.nih.gov/pubmed/?term="



  CANCER_diff_normUQ_level_tf <- CANCER_diff_normUQ_level_tf[1:topgenes,]
  Pubmed <- matrix(0, nrow(CANCER_diff_normUQ_level_tf), 1)
  PMID <- matrix(0, nrow(CANCER_diff_normUQ_level_tf), 1)


  CANCER_diff_normUQ_level_tf <- cbind(CANCER_diff_normUQ_level_tf,Pubmed,PMID)
  CANCER_diff_normUQ_level_tf<-as.data.frame(CANCER_diff_normUQ_level_tf)

  for (k in 1:nrow( CANCER_diff_normUQ_level_tf)){
    CurrentGene <- CANCER_diff_normUQ_level_tf$Gene[k]
    site2 <- paste(site,CurrentGene, "+", tumor,sep="")

    if(interactive() && ("ssl" %in% names(curlVersion()$features)) && url.exists(site2)) {
      x = tryCatch(getURL(site2), error = function(e) {
        getURL(site2, ssl.verifypeer = FALSE) })
    }



    if ( length(grep("No items found.",x))!=1){

      if (length(grep("Display Settings",x))==1){
        x6 <- 1
        CANCER_diff_normUQ_level_tf[k,"PMID"] <- substr(gsub("</dt> <dd>","",unlist(strsplit(x,"PMID:"))[2]),1,8)

      }


      if (length(grep("result_count",x))==1){
        x2a <- unlist(strsplit(x,"result_count"))[2]

        tmpPMID2 <- unlist(strsplit(x2a,"UidCheckBox"))
        tmpPMID3 <- tmpPMID2[grep("<span>",tmpPMID2)]
        CANCER_diff_normUQ_level_tf[k,"PMID"]<- as.character(paste(substr(tmpPMID3,1,8),collapse="; "))


        x3a <-  unlist(strsplit(x2a,"</h2>"))[1]

        if( length(grep("of",x3a))!=1){
          x6 <- as.numeric(unlist(strsplit(x3a,": "))[2])
        } else { x6 <- as.numeric(unlist(strsplit(x3a,"of "))[2]) }
      }


      if (length(grep("following term was not found",x))==1){       x6 <- 0     }

      if (length(grep("Search instead for",x))==1){       x6 <- 1     }

      if (CurrentGene =="JUN"){       x6 <- 0     }
      if (CurrentGene =="HR"){       x6 <- 0     }
      if (CurrentGene =="HOMEZ"){       x6 <- 0     }
      if (CurrentGene =="ANKAR"){       x6 <- 0     }
      if (CurrentGene =="REST"){       x6 <- 0     }
      if (CurrentGene =="BATF"){       x6 <- 0     }
      if (CurrentGene =="MAX"){       x6 <- 0     }
      # if (CurrentGene =="FOS"){       x6 <- 0     }
      if (CurrentGene =="ECD"){       x6 <- 0     }
      # HR, HOMEZ, ANKAR, REST

      CANCER_diff_normUQ_level_tf[k,"Pubmed"]<-x6
      print(paste("Cancer ", tumor, "with TF n. ",k, "of " ,nrow( CANCER_diff_normUQ_level_tf)," : ", CurrentGene, "found n. ", x6, "pubmed."))

    }

    else{
      print(paste("Cancer ", tumor, "with TF n. ",k, "of " ,nrow( CANCER_diff_normUQ_level_tf)," : ", CurrentGene, "no item found in pubmed."))
      CANCER_diff_normUQ_level_tf[k,"Pubmed"]<- 0
    }

  }

  CANCER_diff_normUQ_level_tf <- CANCER_diff_normUQ_level_tf[order(CANCER_diff_normUQ_level_tf$Pubmed,decreasing=T),]
  CANCER_diff_normUQ_level_tf[CANCER_diff_normUQ_level_tf$Pubmed == 1,][ which( nchar(CANCER_diff_normUQ_level_tf[CANCER_diff_normUQ_level_tf$Pubmed == 1,]$PMID) > 8),"PMID"] <- substr(CANCER_diff_normUQ_level_tf[CANCER_diff_normUQ_level_tf$Pubmed == 1,][ which( nchar(CANCER_diff_normUQ_level_tf[CANCER_diff_normUQ_level_tf$Pubmed == 1,]$PMID) > 8),"PMID"],1,8)
  CANCER_diff_normUQ_level_tf[CANCER_diff_normUQ_level_tf$Pubmed == 0,"PMID"]<-0



  return(CANCER_diff_normUQ_level_tf)
}
