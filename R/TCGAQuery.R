#' @title TCGA query
#'
#' @description  Crossfinding of file locations for downloading (TCGADownload)
#' @param tumor tumor code between "acc"  "blca" "brca" "cesc" "chol" "cntl" "coad" "dlbc" "esca" "fppp" "gbm"
#'                                 "hnsc" "kich" "kirc" "kirp" "laml" "lcml" "lgg"  "lihc" "lnnh" "luad" "lusc"
#'                                 "meso" "misc" "ov"   "paad" "pcpg" "prad" "read" "sarc" "skcm" "stad" "tgct"
#'                                 "thca" "thym" "ucec" "ucs"  "uvm"
#'
#' @param platform platform code between  "clin"                                "bio"
#'                                        "biotab"                              "diagnostic_images"
#'                                        "pathology_reports"                   "tissue_images"
#'                                        "illuminahiseq_mirnaseq"              "genome_wide_snp_6"
#'                                        "humanmethylation450"                 "mda_rppa_core"
#'                                        "illuminahiseq_rnaseqv2"              "illuminaga_dnaseq_curated"
#'                                        "illuminaga_dnaseq_automated"         "illuminaga_dnaseq_cont_automated"
#'                                        "mixed_dnaseq_curated"                "illuminaga_mirnaseq"
#'                                        "illuminahiseq_dnaseqc"               "illuminahiseq_wgbs"
#'                                        "illuminahiseq_rnaseq"                "illuminahiseq_totalrnaseqv2"
#'                                        "illuminaga_dnaseq"                   "humanmethylation27"
#'                                        "agilentg4502a_07_3"                  "illuminahiseq_dnaseq_automated"
#'                                        "illuminahiseq_dnaseq_cont_automated" "miRNASeq"
#'                                        "microsat_i"                          "illuminaga_rnaseq"
#'                                        "illuminaga_rnaseqv2"                 "RNASeq"
#'                                        "solid_dnaseq"                        "mixed_dnaseq_automated"
#'                                        "minbio"                              "abi"
#'                                        "ht_hg-u133a"                         "hg-cgh-244a"
#'                                        "hg-cgh-415k_g4124a"                  "humanhap550"
#'                                        "illuminadnamethylation_oma002_cpi"   "illuminadnamethylation_oma003_cpi"
#'                                        "huex-1_0-st-v2"                      "agilentg4502a_07_1"
#'                                        "agilentg4502a_07_2"                  "h-mirna_8x15k"
#'                                        "h-mirna_8x15kv2"                     "mixed_dnaseq_cont"
#'                                        "mixed_dnaseq"                        "mixed_dnaseq_cont_curated"
#'                                        "hg-u133_plus_2"                      " "
#'                                        "human1mduo"                          "cgh-1x1m_g4447a"
#'                                        "illuminaga_mrna_dge"                 "solid_dnaseq_curated"
#'
#' @param added.since Date to search for files (since)
#' @param added.up.to Date to search for files (up.to)
#' @param level level 1 2 3
#' @param listSample list of barcodes to be considered in the search
#'
#' @param qOutput place where the query is saved to be downloaded automatically.
#'        The folder can be specified in both TCGAQuery and TCGADownload
#'
#' @examples
#' \dontrun{
#'   TCGAQuery(tumor = "all",centerType = "all",center = "all",
#'             platform = "all",level = "all",version = "all",
#'             i = F,file = "data/dataFolders.rda",
#'             qOutput = "data/query/")
#' }
#'
#' @author Davide
#' @seealso TCGADownload
#' @export
#' @import downloader
TCGAQuery <- function(tumor="all",
                       platform="all",
                       added.since=NULL,
                       added.up.to=NULL,
                       listSample=NULL,
                       level = "all"
){
  x <- tcga.db
  if(tumor!="all"){
    x <- subset(x, tolower(Disease) == tolower(tumor))
  }
  if(platform!="all"){
    x <- subset(x, tolower(Platform) == tolower(platform))
  }
  if(level!="all"){
    x <- subset(x, grepl(paste0("Level_",level),name) )
  }
  if(!is.null(added.since)){
    x <- subset(x, as.Date(addedDate,"%m/%d/%Y") > as.Date(added.since,"%m/%d/%Y"))
  }
  if(!is.null(added.up.to)){
    x <- subset(x, as.Date(addedDate,"%m/%d/%Y") < as.Date(added.up.to,"%m/%d/%Y"))
  }

  # to be improved
  if(!is.null(listSample)){
    for(i in seq_along(listSample)){
      aux <- grep(listSample[i],tcga.db$deployStatus)
      if(!exists("idx")){
        idx <- aux
      } else {
        idx <- union(idx, aux)
      }
    }
    x <- x[idx,]
  }
  return(x)
}
