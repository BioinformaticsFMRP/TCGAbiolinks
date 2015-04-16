#' @title TCGA DownstreamAnalysisGE
#'
#' @description TCGA DownstreamAnalysis using gene expression data as rnaseqv2 platform with complete pipeline (Query, Download, Prepare, Normalization, Filtering, Differential expression Analysis, Pubmed investigate)
#' @param tumor tumor code between "acc"  "blca" "brca" "cesc" "chol" "cntl" "coad" "dlbc" "esca" "fppp" "gbm"
#'                                 "hnsc" "kich" "kirc" "kirp" "laml" "lcml" "lgg"  "lihc" "lnnh" "luad" "lusc"
#'                                 "meso" "misc" "ov"   "paad" "pcpg" "prad" "read" "sarc" "skcm" "stad" "tgct"
#'                                 "thca" "thym" "ucec" "ucs"  "uvm", or tissues such as "breast", "liver", "colon", etc.
#'
#' @param PathFileTumor PathFileTumor
#' @param CancerName CancerName
#' @param Platform_type Platform_type
#' @param FileManifest FileManifest
#' @param tumorPubmed tumorPubmed
#' @param IPA_transcription_factors IPA_transcription_factors
#'
#'
#' @examples
#' \dontrun{
#' TCGADownstreamAnalysisGE(PathFileTumor,CancerName,Platform_type,FileManifest,tumorPubmed,IPA_transcription_factors)
#' }
#'
#' @author Colaprico et al.
#' @seealso TCGADownload
#' @export
#' @import Rcurl

TCGADownstreamAnalysisGE <- function(PathFileTumor,CancerName,Platform_type,FileManifest,tumorPubmed,IPA_transcription_factors){
  library(RCurl)

  Data_CANCER <- CreateMatrixRnaseqV2(PathFileTumor,CancerName,Platform_type,FileManifest)
  print(paste("STEP1..Matrix RnaseqV2 for ", CancerName, " with n. ", nrow(Data_CANCER), " genes and ", ncol(Data_CANCER),"samples...OK" ))

  Data_CANCER_normUQ <- RnaSeqNormalization(Data_CANCER,geneInfo)
  print(paste("STEP2..Normalization for ", CancerName, " with n. ", nrow(Data_CANCER_normUQ), " genes and ", ncol(Data_CANCER_normUQ),"samples...OK" ))

  Data_CANCER_normUQ_filt <- RnaSeqFilt(Data_CANCER_normUQ,0.25)
  print(paste("STEP3..Filter quantile for ", CancerName, " with n. ", nrow(Data_CANCER_normUQ_filt), " genes and ", ncol(Data_CANCER_normUQ_filt),"samples...OK" ))

  Data_CANCER_normUQ_filt_tumor <- SelectedSample(Data_CANCER_normUQ_filt,"tumor")
  Data_CANCER_normUQ_filt_normal <- SelectedSample(Data_CANCER_normUQ_filt,"normal")

  print(paste("STEP4..Sample found for ", CancerName, " with n. ", ncol(Data_CANCER_normUQ_filt_tumor), " tumor and ", ncol(Data_CANCER_normUQ_filt_normal),"normal samples...OK" ))

  if ( length(Data_CANCER_normUQ_filt_normal)!=0){
    perc <- ncol(Data_CANCER_normUQ_filt_normal)/ncol(Data_CANCER_normUQ_filt_tumor)
  } else { perc <- 0}


  if( perc >= 0.05 ){
    print(paste("Finding TR candidate with DEA with perc ", perc, sep =""))
    CANCER_diff_normUQ <- DEA_edge5(Data_CANCER_normUQ_filt_normal,Data_CANCER_normUQ_filt_tumor, "Normal", "Tumor")
    print(paste("STEP5.. DEA for ", CancerName, " completed."))

    CANCER_diff_normUQ_level <- CreateTabLevel( as.matrix(rownames(CANCER_diff_normUQ)),CANCER_diff_normUQ,"Tumor","Normal",Data_CANCER_normUQ_filt_tumor,Data_CANCER_normUQ_filt_normal,typeOrder=T)
    CANCER_diff_normUQ_level_tf <- CANCER_diff_normUQ_level[rownames(CANCER_diff_normUQ_level) %in% IPA_transcription_factors$Gene,]
    CANCER_diff_normUQ_level_tf <- CANCER_diff_normUQ_level_tf[order(CANCER_diff_normUQ_level_tf$Delta,decreasing=T),]
    CANCER_diff_normUQ_level_tf <- CANCER_diff_normUQ_level_tf[abs(CANCER_diff_normUQ_level_tf$logFC) >= 1,]

    CANCER_diff_normUQ_level_tf_pubmed <- FindPubmedTFgene(tumorPubmed,CANCER_diff_normUQ_level_tf,topgenes=nrow(CANCER_diff_normUQ_level_tf))

    print(paste("STEP6.. Pubmed for ", CancerName, " completed."))

    ans <- list(Data_CANCER = Data_CANCER, Data_CANCER_normUQ = Data_CANCER_normUQ,  Data_CANCER_normUQ_filt = Data_CANCER_normUQ_filt, CANCER_diff_normUQ = CANCER_diff_normUQ, CANCER_diff_normUQ_level = CANCER_diff_normUQ_level, CANCER_diff_normUQ_level_tf = CANCER_diff_normUQ_level_tf, CANCER_diff_normUQ_level_tf_pubmed = CANCER_diff_normUQ_level_tf_pubmed)
  }

  else {
    print(paste("Finding TR candidate with Top Expressed with perc ", perc, sep =""))
    CANCER_diff_normUQ_level_tf <- Data_CANCER_normUQ_filt[rownames(Data_CANCER_normUQ_filt) %in% IPA_transcription_factors$Gene,]
    tmp <- sort(rowMeans(CANCER_diff_normUQ_level_tf),decreasing=T)
    tmp <- as.matrix(tmp)
    colnames(tmp)<-"Tumor"
    tmp <- cbind(Gene = rownames(tmp),round(tmp))
    tmp <- as.data.frame(tmp)
    tmp$Tumor <-as.numeric( as.character(tmp$Tumor))
    tmp <-tmp[ tmp$Tumor > quantile(as.numeric(tmp$Tumor),0.50),]

    tmp2 <- FindPubmedTFgene2(tumorPubmed,tmp,topgenes=nrow(tmp))
    #tmp <- cbind(tmp, Pubmed = tmp2)
    CANCER_diff_normUQ_level_tf_pubmed <- tmp2

    ans <- list(Data_CANCER = Data_CANCER, Data_CANCER_normUQ = Data_CANCER_normUQ,  Data_CANCER_normUQ_filt = Data_CANCER_normUQ_filt,CANCER_diff_normUQ_level_tf = CANCER_diff_normUQ_level_tf, CANCER_diff_normUQ_level_tf_pubmed = CANCER_diff_normUQ_level_tf_pubmed )

  }

  return(ans)
}
