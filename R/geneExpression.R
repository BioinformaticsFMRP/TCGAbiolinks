#' @title GenesCutID
#' @description
#'   GenesCutID
#' @param GeneList GeneList
#' @export
#' @return list of gene symbol without IDs
GenesCutID <- function(GeneList){
  GeneListCutID<- as.matrix(matrix(unlist(strsplit(as.character(GeneList),"|",fixed=T)),nrow(GeneList),2,byrow=T))[,1]
  return(as.matrix(GeneListCutID))
}

#' @title TimeUse
#' @description
#'   TimeUse
#' @param func func
#' @export
#' @return time execution of a function
TimeUse <-function(func){
  time1 <- proc.time() # mod1 to determine time calculation
  result <-func
  time2 <- proc.time()
  show(timeelapsed <- time2-time1)
}

#' @title Filtering rnaseq by quantile
#' @description
#'   Filtering rnaseq by quantile
#' @param TableRnaseq TableRnaseq
#' @param QuantileThresh QuantileThresh
#' @export
#' @return table filtered
RnaSeqFilt<- function(TableRnaseq,QuantileThresh ){
  GeneThresh <- as.numeric(quantile(rowMeans(TableRnaseq), QuantileThresh))
  geneFiltered <- names(which(rowMeans(TableRnaseq) > GeneThresh))
  Table_Rnaseq_Rawcount_Filt <- TableRnaseq[geneFiltered, ]
  return( Table_Rnaseq_Rawcount_Filt)
}

#' @title RnaSeqNormalization
#' @description
#'   RnaSeqNormalization
#' @param TCGA_RnaseqTable TCGA_RnaseqTable
#' @param geneInfo geneInfo
#' @importFrom EDASeq newSeqExpressionSet withinLaneNormalization betweenLaneNormalization counts
#' @export
#' @return table normalized
RnaSeqNormalization <- function(TCGA_RnaseqTable,geneInfo){

  rownames(TCGA_RnaseqTable) <- GenesCutID(as.matrix(rownames(TCGA_RnaseqTable)))
  TCGA_RnaseqTable <- TCGA_RnaseqTable[rownames(TCGA_RnaseqTable) != "?", ]
  TCGA_RnaseqTable<-TCGA_RnaseqTable[!duplicated(rownames(TCGA_RnaseqTable)), !duplicated(colnames(TCGA_RnaseqTable))]
  TCGA_RnaseqTable <- TCGA_RnaseqTable[, which(substr(colnames(TCGA_RnaseqTable), 14, 15) != "02")]
  geneInfo <- geneInfo[rownames(geneInfo) %in% rownames(TCGA_RnaseqTable), ]
  geneInfo <- geneInfo[!duplicated(rownames(geneInfo)), ]
  toKeep <- which(geneInfo[, "geneLength"] != 0)
  geneInfo <- geneInfo[toKeep, ]
  TCGA_RnaseqTable <- TCGA_RnaseqTable[toKeep, ]
  geneInfo <- as.data.frame(geneInfo)
  TCGA_RnaseqTable<-round(TCGA_RnaseqTable)

  timeEstimated<-format(ncol(TCGA_RnaseqTable)*nrow(TCGA_RnaseqTable)/80000,digits=2)
  print(messageEstimation<-paste("I Need about ", timeEstimated, "seconds for this Complete Normalization Upper Quantile [Processing 80k elements /s]  "))

  print("Step 1 of 4: newSeqExpressionSet ...")
  TimeUse(TCGA_RnaseqTable_norm <- newSeqExpressionSet(TCGA_RnaseqTable, featureData = geneInfo))
  print("Step 2 of 4: withinLaneNormalization ...")
  TimeUse(TCGA_RnaseqTable_norm <- withinLaneNormalization(TCGA_RnaseqTable_norm, "geneLength", which = "upper", offset = FALSE))
  print("Step 3 of 4: betweenLaneNormalization ...")
  TimeUse(TCGA_RnaseqTable_norm <- betweenLaneNormalization(TCGA_RnaseqTable_norm, which = "upper", offset = FALSE))
  print("Step 4 of 4: exprs ...")
  #TimeUse(TCGA_RnaseqTable_norm <- exprs(TCGA_RnaseqTable_norm))
  TimeUse(TCGA_RnaseqTable_norm <- counts(TCGA_RnaseqTable_norm))

  return(TCGA_RnaseqTable_norm)
}

#' @title Differentially expression analysis (DEA)
#' @description
#'    Differentially expression analysis (DEA) with gene expression
#' @param mat1 datamatrix normal
#' @param mat2 datamatrix tumor
#' @param Cond1type normal
#' @param Cond2type tumor
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags
#' @export
#' @return table with DEGs (diff.expr. genes)
DEArnaSEQ <- function(mat1,mat2,Cond1type,Cond2type) {

  TOC <- cbind(mat1,mat2)
  Cond1num <- ncol(mat1)
  Cond2num <- ncol(mat2)

  print(message1<-paste( "there are Cond1 type", Cond1type ,"in ", Cond1num, "samples"))
  print(message2<-paste( "there are Cond2 type", Cond2type ,"in ", Cond2num, "samples"))
  print(message3<-paste( "there are ", nrow(TOC) ,"species miRNA or genes "))

  timeEstimated<-format(ncol(TOC)*nrow(TOC)/30000,digits=2)
  print(messageEstimation<-paste("I Need about ", timeEstimated, "seconds for this DEA. [Processing 30k elements /s]  "))

  # Reading in the data and creating a DGEList object
  colnames(TOC) <- paste('s',1:ncol(TOC),sep="")
  #DGE <- DGEList(TOC,group=rep(c("Normal","Tumor"),c(NormalSample,TumorSample)))
  DGE <- DGEList(TOC,group=rep(c(Cond1type,Cond2type),c(Cond1num,Cond2num)))

  # Analysis using common dispersion
  disp <- estimateCommonDisp(DGE) # Estimating the common dispersion
  #tested <- exactTest(disp,pair=c("Normal","Tumor")) # Testing
  tested <- exactTest(disp,pair=c(Cond1type,Cond2type)) # Testing

  # Results visualization
  logFC_table <- tested$table
  logFC_FDR_table <- topTags(tested,n=nrow(tested$table))$table
  return(logFC_FDR_table)

}
