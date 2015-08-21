#' @title Visaulize results in format of latex tables.
#' @description Visaulize results in format of latex tables.
#' @param Table write
#' @param rowsForPage write
#' @param TableTitle write
#' @param LabelTitle write
#' @param withrows write
#' @importFrom xtable xtable
#' @examples
#' \dontrun{
#' library(stringr)
#' tabDEGsTFPubmed$PMID <- str_sub(tabDEGsTFPubmed$PMID,0,30)
#' TCGAvisualize_Tables(Table = tabDEGsTFPubmed,
#' rowsForPage = 5,
#' TableTitle = "pip",
#' LabelTitle = "pip2",
#' withrows = FALSE)
#' }
#' @export
#' @return table in latex format to use in beamer presentation or sweave files
TCGAvisualize_Tables <- function(Table, rowsForPage, TableTitle, LabelTitle, withrows){
    require(xtable)

    numberOfprint<-ceiling(nrow(Table)/rowsForPage)
    vectorFirst<-matrix(0,numberOfprint,1)
    vectorLast<-matrix(0,numberOfprint,1)

    i<-0
    for(i in 1:numberOfprint){
        vectorFirst[i]<-rowsForPage*i
        vectorLast[i]<-rowsForPage*i+1
    }

    i<-1
    for(i in 1:numberOfprint){
       if (round(nrow(Table)/rowsForPage) == 1 ) {

            Table_current<-Table[i:nrow(Table),]
            tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size=small)
            print(tablePrint_Table_current,include.rownames = withrows)
        }

        else{
            print(i)
            if( i == 1 ) {
                Table_current<-Table[i:vectorFirst[i],]
                tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size=small)
                print(tablePrint_Table_current,include.rownames = withrows)
            }

            else if (i==numberOfprint) {
                Table_current<-Table[vectorLast[i-1]:nrow(Table),]
                tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size=small)
                print(tablePrint_Table_current,include.rownames = withrows)
            }
            else{
                Table_current<-Table[vectorLast[i-1]:vectorFirst[i],]
                tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size=small)
                print(tablePrint_Table_current,include.rownames = withrows)
            }
        }
    }

}

#' @title Heatmap with more sensible behavior using heatmap.plus
#' @description Heatmap with more sensible behavior using heatmap.plus
#' @param DFfilt write
#' @param DFclin write
#' @param DFsubt write
#' @param data_Hc2 write
#' @param Subtype write
#' @param cbPalette write
#' @param filename write. default = NULL
#' @importFrom heatmap.plus heatmap.plus
#' @examples
#' \dontrun{
#' # from case study n.2 LGG to test the function
#' DFfilt <- datFilt
#' DFclin = dataClin
#' DFsubt = dataSubt
#' data_Hc2 = data_Hc2
# end parameter definition
#' TCGAvisualize_Heatmap(DFfilt,
#' DFclin,
#' DFsubt,
#' data_Hc2)
#' }
#' @export
#' @return Heatmap plotted in pdf or png file.
TCGAvisualize_Heatmap <- function(DFfilt, DFclin, DFsubt, data_Hc2, cbPalette, filename =NULL){

    rownames(DFsubt) <- DFsubt$patient
    rownames(DFclin) <- DFclin$patient
    rownames(DFfilt) <- substr(rownames(DFfilt),1,12)

    ans <- hclust(ddist <- dist(DFfilt), method = "ward.D2")
    hhc <- data_Hc2[[4]]$consensusTree
    consensusClusters<-data_Hc2[[4]]$consensusClass
    sampleOrder <- consensusClusters[hhc$order]

    consensusClusters <- as.factor(data_Hc2[[4]]$clrs[[1]])
    names(consensusClusters) <- attr(ddist, "Labels")
    names(consensusClusters) <- substr(names(consensusClusters),1,12)




    #DFclin <- DFclin[DFclin$bcr_patient_barcode %in% DFsubt$patient,]
    DFclin_merged <- cbind(DFclin, matrix(0,nrow(DFclin),ncol(DFsubt)))
    colnames(DFclin_merged)[((ncol(DFclin_merged)-ncol(DFsubt))+1) :ncol(DFclin_merged)] <- colnames(DFsubt)
    rownames(DFclin_merged) <- DFclin_merged$bcr_patient_barcode

    for( i in 1: ncol(DFsubt)){
        DFsubt[,i] <- as.character(DFsubt[,i])
    }

    for( i in 1: nrow(DFsubt)){
        curSample <- DFsubt$patient[i]
        for( j in 1: ncol(DFsubt)){
            curColumn <- colnames(DFsubt)[j]
            DFclin_merged[curSample,curColumn] <- DFsubt[curSample,curColumn]
        }
    }


    # adding information about gropus from consensus Cluster in clinical data
    DFclin_merged <- cbind(DFclin_merged, groupsHC = matrix(0,nrow(DFclin_merged),1))
    rownames(DFclin_merged) <- DFclin_merged$bcr_patient_barcode

    for( i in 1: nrow(DFclin_merged)){
        currSmp <- DFclin_merged$bcr_patient_barcode[i]
        DFclin_merged[currSmp,"groupsHC"] <- as.character(consensusClusters[currSmp])
    }

    groupsColors <-  levels(as.factor(DFclin_merged$groupsHC))


    for(j in 1:length(table(DFclin_merged$groupsHC))){
        curCol <- groupsColors[j]
        DFclin_merged[DFclin_merged$groupsHC == curCol,"groupsHC"]<-paste0("EC",j)
    }


    DFfilt <- DFfilt[rownames(DFclin_merged),]

    orderCL <- as.character(substr(names(sampleOrder),1,12))
    orderCL <- intersect(orderCL, rownames(DFfilt))
    GE <- t(TCGAbiolinks:::.quantileNormalization(t(DFfilt)))
    rownames(GE) <- substr(rownames(GE),1,12)

    oGE<- GE[orderCL,]  #ordering according cluster
    DFclin_merged <-DFclin_merged[orderCL,]

    # histology
    HISTOLOGY <- DFclin_merged[,"histological_type"]
    names(HISTOLOGY)<-rownames(DFclin_merged)
    HISTOLOGY <- HISTOLOGY[rownames(GE)]
    HISTOLOGY.col <- rep("white",length(HISTOLOGY))
    HISTOLOGY.col[HISTOLOGY=="Astrocytoma"]<-"red"
    HISTOLOGY.col[HISTOLOGY=="glioblastoma"]<-"purple"
    HISTOLOGY.col[HISTOLOGY=="Oligoastrocytoma"]<-"cyan"
    HISTOLOGY.col[HISTOLOGY=="Oligodendroglioma"]<-"green3"
    names(HISTOLOGY.col)<-names(HISTOLOGY)
    oHISTOLOGY.col <- HISTOLOGY.col[orderCL]

    #subtype
    SUBTYPE <- DFclin_merged[,"IDH.1p19q.Subtype"]
    names(SUBTYPE) <- rownames(DFclin_merged)
    SUBTYPE<-SUBTYPE[rownames(GE)]
    SUBTYPE.col <- rep("white",length(SUBTYPE))
    SUBTYPE.col[SUBTYPE=="IDHmut-codel"]<-"cyan"
    SUBTYPE.col[SUBTYPE=="IDHmut-non-codel"]<-"tomato"
    SUBTYPE.col[SUBTYPE=="IDHwt"]<-"gold"
    names(SUBTYPE.col)<-names(SUBTYPE)
    oSUBTYPE.col <- SUBTYPE.col[orderCL]

    #clusters CNCluster
    CNC <- DFclin_merged[,"CNCluster"]
    names(CNC)<- rownames(DFclin_merged)
    CNC <- CNC[rownames(GE)]
    CNC.col <- rep("white",length(CNC))
    names(CNC.col)<-names(CNC)
    CNC.col[CNC=="C1"] <- "green"
    CNC.col[CNC=="C2"] <- "red"
    CNC.col[CNC=="C3"] <- "purple"
    oCNC.col <- CNC.col[orderCL]

    #clusters COCluster
    COC <- DFclin_merged[,"COCCluster"]
    names(COC)<- rownames(DFclin_merged)
    COC <- COC[rownames(GE)]
    COC.col <- rep("white",length(COC))
    names(COC.col)<-names(COC)
    COC.col[COC=="coc1"] <- "green"
    COC.col[COC=="coc2"] <- "red"
    COC.col[COC=="coc3"] <- "purple"
    oCOC.col <- COC.col[orderCL]

    #clusters ONCOluster
    ONCO <- DFclin_merged[,"OncosignCluster"]
    names(ONCO)<- rownames(DFclin_merged)
    ONCO <- ONCO[rownames(GE)]
    ONCO.col <- rep("white",length(ONCO))
    names(ONCO.col)<-names(ONCO)
    ONCO.col[ONCO=="OSC1"] <- "green"
    ONCO.col[ONCO=="OSC2"] <- "red"
    ONCO.col[ONCO=="OSC3"] <- "purple"
    ONCO.col[ONCO=="OSC4"] <- "orange"
    ONCO.col[ONCO=="Unclassified"] <- "gray"
    oONCO.col <- ONCO.col[orderCL]

    oConsensus <- as.character(consensusClusters[hhc$order])
    #oConsensus <- as.character(consensusClusters[orderCL])

    names(consensusClusters[hhc$order])


    #source("heatmap.plus.R")

    cc.col <- matrix(c(oHISTOLOGY.col,
                       oSUBTYPE.col,
                       oCNC.col,
                       oCOC.col,
                       oONCO.col,
                       as.character(oConsensus)),
                     nrow = nrow(oGE), ncol = 6)

    colnames(cc.col)<-c("Histology",
                        "Subtype",
                        "CNCluster",
                        "COCCluster",
                        "OncosignCluster",
                        "Expression Cluster")

    cc.col <- as.data.frame(cc.col)
    rownames(cc.col) <- orderCL
    cc.col <- cc.col[order(cc.col$`Expression Cluster`),]
    cc.col <- as.matrix(cc.col)
    oGE<- oGE[rownames(cc.col),]

    #source('~/Downloads/heatmap.plus2.R')
    #take from here https://github.com/mercutio22/kate/blob/master/heatmap.plus.R

    if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}

    curDate <- as.character(unlist(strsplit(gsub(" ","_h",
                                                 gsub("-","_",as.character(Sys.time()))),":"))[1])


    pdf(file=paste0(curDate,"_",cancer,"_heatmap_with_subtypes_withHeatmapPlus.pdf"))

    heatmap.plus.sm(
        t(oGE),
        na.rm=TRUE,
        scale="none",
        #RowSideColor=probe.cc,
        #ColSideColors=cc.col,
        col=greenred(75),
        key=FALSE,  #changed
        symkey=FALSE,
        density.info="none",
        trace="none",
        Rowv=FALSE,
        Colv=NA,
        cexRow=1,
        cexCol=1.6,
        keysize=2,
        dendrogram = "none",
        main = "Heatmap from consensus cluster",
        labRow=NA,labCol=NA,
        #labCol=NA
    )
    dev.off()
}
