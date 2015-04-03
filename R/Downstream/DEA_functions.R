

DEA_edge5<- function (mat1,mat2,Cond1type,Cond2type) {
  library(edgeR)

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

CreateTabLevel<-function(TF_enriched,FC_FDR_table_mRNA,typeCond1,typeCond2,TableCond1,TableCond2,typeOrder) {

  TableLevel<-matrix(0,nrow(TF_enriched),6)
  TableLevel <- as.data.frame(TableLevel)

  colnames(TableLevel)<-c("mRNA","logFC","FDR",typeCond1,typeCond2,"Delta")


  TableLevel[,"mRNA"]<-TF_enriched
  Tabfilt<-FC_FDR_table_mRNA[which( rownames(FC_FDR_table_mRNA) %in% TF_enriched),]
  TableLevel[,"logFC"]<- as.numeric(Tabfilt[TF_enriched,][,"logFC"])
  TableLevel[,"FDR"]<- as.numeric(Tabfilt[TF_enriched,][,"FDR"])


  MeanTumor<-matrix(0,nrow(TF_enriched),1)
  MeanDiffTumorNormal<-matrix(0,nrow(TF_enriched),1)


  for( i in 1:nrow(TF_enriched)) {
    #print(paste(i, "of", nrow(TF_enriched),TF_enriched[i]))
    TableLevel[i,typeCond1]<- mean(TableCond1[rownames(TableCond1) %in%  TF_enriched[i] , ])
    TableLevel[i,typeCond2]<- mean(TableCond2[rownames(TableCond2) %in%  TF_enriched[i] , ])
  }


  TableLevel[,"Delta"] <- as.numeric(abs(TableLevel[,"logFC"]) * TableLevel[,typeCond1]  )

  TableLevel<-TableLevel[order( as.numeric(TableLevel[,"Delta"]),decreasing=typeOrder),]

  rownames(TableLevel) <-  TableLevel[,"mRNA"]


  return(TableLevel)
}

