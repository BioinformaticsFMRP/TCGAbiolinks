#' @title GenesCutID
#' @description
#'   GenesCutID
#' @param GeneList GeneList
#' @export
#' @return list of gene symbol without IDs
#' @examples
#' GenesCutID(c("CRKL|1399","TADA2A|6871","KRT76|51350"))
GenesCutID <- function(GeneList){
    GeneListCutID <- as.matrix(matrix(unlist(strsplit(as.character(GeneList),"|",fixed = TRUE)),length(GeneList),2,byrow = TRUE))[,1]
    return(as.matrix(GeneListCutID))
}

#' @title TimeUse
#' @description
#'   TimeUse
#' @param func func
#' @export
#' @examples
#'  TimeUse(print(1))
#' \dontrun{
#'  Genelist <-rownames(dataDEGsFiltLevel)
#'  TimeUse(ansEA <- EAcomplete(TFname="DEA genes Normal Vs Tumor",
#'  Genelist))
#' }
#' @return time execution of a function
TimeUse <- function(func){
    time1 <- proc.time() # mod1 to determine time calculation
    result <- func
    time2 <- proc.time()
    show(timeelapsed <- time2 - time1)
}

#' @title Filtering rnaseq by quantile
#' @description
#'   Filtering rnaseq by quantile
#' @param TableRnaseq TableRnaseq
#' @param QuantileThresh QuantileThresh
#' @export
#' @return table filtered
#' @examples
#' \dontrun{
#' dataNorm <- TCGAbiolinks::RnaSeqNormalization(dataBRCA, geneInfo)
#' dataFilt <- RnaSeqFilt(dataNorm, 0.25)
#'}
RnaSeqFilt <- function(TableRnaseq,QuantileThresh ){
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
#' @importFrom EDASeq newSeqExpressionSet withinLaneNormalization betweenLaneNormalization exprs counts
#' @export
#' @return table normalized
#' @examples
#' \dontrun{
#' dataNorm <- TCGAbiolinks::RnaSeqNormalization(dataBRCA, geneInfo)
#'}
RnaSeqNormalization <- function(TCGA_RnaseqTable,geneInfo){

    TCGA_RnaseqTable <- TCGA_RnaseqTable[ !(GenesCutID(as.matrix(rownames(TCGA_RnaseqTable))) == "?"),]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[ !(GenesCutID(as.matrix(rownames(TCGA_RnaseqTable))) == "SLC35E2"),]
    rownames(TCGA_RnaseqTable) <- GenesCutID(as.matrix(rownames(TCGA_RnaseqTable)))
    TCGA_RnaseqTable <- TCGA_RnaseqTable[rownames(TCGA_RnaseqTable) != "?", ]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[!duplicated(rownames(TCGA_RnaseqTable)), !duplicated(colnames(TCGA_RnaseqTable))]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[, which(substr(colnames(TCGA_RnaseqTable), 14, 15) != "02")]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[rownames(TCGA_RnaseqTable) %in% rownames(geneInfo),]
    TCGA_RnaseqTable <- as.matrix(TCGA_RnaseqTable)

    geneInfo <- geneInfo[rownames(geneInfo) %in% rownames(TCGA_RnaseqTable), ]
    geneInfo <- geneInfo[!duplicated(rownames(geneInfo)), ]
    toKeep <- which(geneInfo[, "geneLength"] != 0)
    geneInfo <- geneInfo[toKeep, ]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[toKeep, ]
    geneInfo <- as.data.frame(geneInfo)
    TCGA_RnaseqTable <- round(TCGA_RnaseqTable)
    commonGenes <- intersect(rownames(TCGA_RnaseqTable),rownames(geneInfo))

    TCGA_RnaseqTable <- TCGA_RnaseqTable[commonGenes,]

    timeEstimated <- format(ncol(TCGA_RnaseqTable)*nrow(TCGA_RnaseqTable)/80000,digits = 2)
    print(messageEstimation <- paste("I Need about ", timeEstimated, "seconds for this Complete Normalization Upper Quantile [Processing 80k elements /s]  "))

    print("Step 1 of 4: newSeqExpressionSet ...")
    TimeUse(TCGA_RnaseqTable_norm <- EDASeq::newSeqExpressionSet(TCGA_RnaseqTable, featureData = geneInfo))
    print("Step 2 of 4: withinLaneNormalization ...")
    TimeUse(TCGA_RnaseqTable_norm <- EDASeq::withinLaneNormalization(TCGA_RnaseqTable_norm, "geneLength", which = "upper", offset = FALSE))
    print("Step 3 of 4: betweenLaneNormalization ...")
    TimeUse(TCGA_RnaseqTable_norm <- EDASeq::betweenLaneNormalization(TCGA_RnaseqTable_norm, which = "upper", offset = FALSE))
    print("Step 4 of 4: exprs ...")

    #TimeUse(TCGA_RnaseqTable_norm <- EDASeq::exprs(TCGA_RnaseqTable_norm))
    TimeUse(TCGA_RnaseqTable_norm <- EDASeq::counts(TCGA_RnaseqTable_norm))

    return(TCGA_RnaseqTable_norm)
}

#' @title Differentially expression analysis (DEA)
#' @description
#'    Perform DEA to identify differentially expressed genes (DEGs). It is possible to do a two-class analysis.
#' @param mat1 numeric matrix, each row represents a gene, each column represents a sample with Cond1type
#' @param mat2 numeric matrix, each row represents a gene, each column represents a sample with Cond2type
#' @param Cond1type a string containing the class label of the samples in mat1  (e.g., control group)
#' @param Cond2type a string containing the class label of the samples in mat2  (e.g., case group)
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags
#' @export
#' @examples
#' \dontrun{
#' library(TCGAbiolinks)
#' dataNorm <- TCGAbiolinks::RnaSeqNormalization(dataBRCA, geneInfo)
#' dataFilt <- RnaSeqFilt(dataNorm, 0.25)
#' samplesNT <- MultiSampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- MultiSampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- DEArnaSEQ(dataFilt[,samplesNT], dataFilt[,samplesTP],"Normal", "Tumor")
#' }
#' @return table containing for each gene logFC, logCPM, pValue,and    FDR
DEArnaSEQ <- function(mat1,mat2,Cond1type,Cond2type) {

    TOC <- cbind(mat1,mat2)
    Cond1num <- ncol(mat1)
    Cond2num <- ncol(mat2)

    print(message1 <- paste( "there are Cond1 type", Cond1type ,"in ",
                             Cond1num, "samples"))
    print(message2 <- paste( "there are Cond2 type", Cond2type ,"in ",
                             Cond2num, "samples"))
    print(message3 <- paste( "there are ", nrow(TOC) ,
                             "features as miRNA or genes "))

    timeEstimated <- format(ncol(TOC)*nrow(TOC)/30000,digits = 2)
    print(messageEstimation <- paste("I Need about ", timeEstimated,
                    "seconds for this DEA. [Processing 30k elements /s]  "))

    # Reading in the data and creating a DGEList object
    colnames(TOC) <- paste0('s',1:ncol(TOC))
    #DGE <- DGEList(TOC,group=rep(c("Normal","Tumor"),c(NormalSample,
    #TumorSample)))
    DGE <- edgeR::DGEList(TOC,group = rep(c(Cond1type,Cond2type),
                                          c(Cond1num,Cond2num)))

    # Analysis using common dispersion
    disp <- edgeR::estimateCommonDisp(DGE) # Estimating the common dispersion
    #tested <- exactTest(disp,pair=c("Normal","Tumor")) # Testing
    tested <- edgeR::exactTest(disp,pair = c(Cond1type,Cond2type)) # Testing

    # Results visualization
    logFC_table <- tested$table
    logFC_FDR_table <- edgeR::topTags(tested,n = nrow(tested$table))$table
    return(logFC_FDR_table)

}

#' @title CreateTabLevel for Differentially expression analysis (DEA)
#' @description
#'    CreateTabLevel for Differentially expression analysis (DEA)
#' @param FC_FDR_table_mRNA Output of dataDEGs filter by abs(LogFC) >=1
#' @param typeCond1 a string containing the class label of the samples in TableCond1  (e.g., control group)
#' @param typeCond2 a string containing the class label of the samples in TableCond2  (e.g., case group)
#' @param TableCond1 numeric matrix, each row represents a gene, each column represents a sample with Cond1type
#' @param TableCond2 numeric matrix, each row represents a gene, each column represents a sample with Cond2type
#' @param typeOrder typeOrder
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags
#' @export
#' @return table with DEGs, log Fold Change (FC), false discovery rate (FDR), the gene expression level
#' for samples in  Cond1type, and Cond2type, and Delta value (the difference of gene expression between the two
#' conditions multiplied logFC)
#' @examples
#' \dontrun{
#' library(TCGAbiolinks)
#' dataNorm <- TCGAbiolinks::RnaSeqNormalization(dataBRCA, geneInfo)
#' dataFilt <- RnaSeqFilt(dataNorm, 0.25)
#' samplesNT <- MultiSampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- MultiSampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- DEArnaSEQ(dataFilt[,samplesNT], dataFilt[,samplesTP],"Normal", "Tumor")
#' dataDEGsFilt <- dataDEGs[abs(dataDEGs$logFC) >= 1,]
#' dataTP <- dataFilt[,samplesTP]
#' dataTN <- dataFilt[,samplesNT]
#' dataDEGsFiltLevel <- CreateTabLevel(dataDEGsFilt,"Tumor","Normal",dataTP,dataTN)
#' }
CreateTabLevel <- function(FC_FDR_table_mRNA,typeCond1,typeCond2,
                           TableCond1,TableCond2,typeOrder = TRUE) {

    TF_enriched <- as.matrix(rownames(FC_FDR_table_mRNA))
    TableLevel <- matrix(0,nrow(TF_enriched),6)
    TableLevel <- as.data.frame(TableLevel)

    colnames(TableLevel) <- c("mRNA","logFC","FDR",typeCond1,typeCond2,"Delta")


    TableLevel[,"mRNA"] <- TF_enriched
    Tabfilt <- FC_FDR_table_mRNA[which( rownames(FC_FDR_table_mRNA) %in%
                                            TF_enriched),]
    TableLevel[,"logFC"] <- as.numeric(Tabfilt[TF_enriched,][,"logFC"])
    TableLevel[,"FDR"] <- as.numeric(Tabfilt[TF_enriched,][,"FDR"])


    MeanTumor <- matrix(0,nrow(TF_enriched),1)
    MeanDiffTumorNormal <- matrix(0,nrow(TF_enriched),1)


    for (i in 1:nrow(TF_enriched)) {
        #print(paste(i, "of", nrow(TF_enriched),TF_enriched[i]))
        TableLevel[i,typeCond1] <- mean(TableCond1[rownames(TableCond1) %in%
                                                       TF_enriched[i] , ])
        TableLevel[i,typeCond2] <- mean(TableCond2[rownames(TableCond2) %in%
                                                       TF_enriched[i] , ])
    }


    TableLevel[,"Delta"] <- as.numeric(abs(TableLevel[,"logFC"]) *
                                           TableLevel[,typeCond1]  )

    TableLevel <- TableLevel[order( as.numeric(TableLevel[,"Delta"]),
                                    decreasing = typeOrder),]

    rownames(TableLevel) <-  TableLevel[,"mRNA"]
    return(TableLevel)
}

#' @title plotPCAforGroups
#' @description
#'   plotPCAforGroups
#' @param dataFilt dataFilt
#' @param dataDEGsFiltLevel dataDEGsFiltLevel
#' @param ntopgenes ntopgenes
#' @import ggplot2
#' @export
#' @return PCA plot
#' @examples
#' \dontrun{
#' # normalization of genes
#' dataNorm <- TCGAbiolinks::RnaSeqNormalization(dataBRCA, geneInfo)
#' # quantile filter of genes
#' dataFilt <- RnaSeqFilt(dataNorm, 0.25)
#' # Principal Component Analysis plot for ntop selected DEGs
#' plotPCAforGroups(dataFilt,dataDEGsFiltLevel, ntopgenes = 200)
#'}
plotPCAforGroups <- function(dataFilt,dataDEGsFiltLevel ,ntopgenes) {
    ComparisonSelected <- "Normal vs Tumor"
    TitlePlot <- paste0("PCA ", "top ", ntopgenes,
                        " Up and down diff.expr genes between ",
                        ComparisonSelected)

    Genelist <- rownames(dataDEGsFiltLevel)[1:ntopgenes]
    commonGenes <- intersect(Genelist, rownames(dataFilt) )
    expr2 <- dataFilt[commonGenes,]
    color1 <- "blue"
    color2 <- "red"

    # selection of normal samples "NT"
    samplesNT <- TCGAbiolinks::MultiSampleTypes(colnames(dataFilt),
                                                typesample = c("NT"))
    # selection of tumor samples "TP"
    samplesTP <- TCGAbiolinks::MultiSampleTypes(colnames(dataFilt),
                                                typesample = c("TP"))

    nsample1 <- length(samplesNT)
    nsample2 <- length(samplesTP)

    #sampleColors <- rep(c(color1,color2), c(nsample1, nsample2))
    #sampleColors <- rep(c("blue","red"), c(length(samplesNT),
    #                     length(samplesTP)))
    sampleColors <- c(rep("blue", length(samplesNT)),
                      rep("red", length(samplesTP)))


    names(sampleColors) <- colnames(expr2)
    cancer.pca <- stats::prcomp(t(expr2),cor = TRUE)

    # print(sample.colors)
    g <- ggbiplot(cancer.pca, obs.scale = 1, var.scale = 1,
                  groups = sampleColors, ellipse = TRUE, circle = FALSE)
    g <- g + scale_colour_manual(name = "",
                                 values = c("blue" = "blue","red" = "red"))
    g <- g + geom_point(aes(colour = sampleColors), size = 3)
    #shape = tabClusterNew$Study)
    g <- g + theme(legend.direction = 'horizontal',  legend.position = 'top')
    g <- g + ggtitle(TitlePlot)
    print(g)
}

#' @title EAcomplete
#' @description
#'   EAcomplete
#' @param TFname TFname
#' @param RegulonList RegulonList
#' @export
#' @return EAcomplete plot
#' @examples
#' \dontrun{
#' Genelist <- rownames(dataDEGsFiltLevel)
#' TimeUse(ansEA <- EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
#' }
EAcomplete <- function(TFname, RegulonList){
    EAGenes <- get("EAGenes")
    DAVID_BP_matrix <- get("DAVID_BP_matrix")
    DAVID_BP_matrix <- get("DAVID_BP_matrix")
    DAVID_MF_matrix <- get("DAVID_MF_matrix")
    DAVID_CC_matrix <- get("DAVID_CC_matrix")
    listEA_pathways <- get("listEA_pathways")

    print(paste("I need about ", "1 minute to finish complete ",
                "Enrichment analysis GO[BP,MF,CC] and Pathways... "))

    ResBP <- EnrichmentAnalysis(TFname,RegulonList,DAVID_BP_matrix,
                                EAGenes,GOtype = "DavidBP")
    print("GO Enrichment Analysis BP completed....done")
    ResMF <- EnrichmentAnalysis(TFname,RegulonList,DAVID_MF_matrix,
                                EAGenes,GOtype = "DavidMF")
    print("GO Enrichment Analysis MF completed....done")
    ResCC <- EnrichmentAnalysis(TFname,RegulonList,DAVID_CC_matrix,
                                EAGenes,GOtype = "DavidCC")
    print("GO Enrichment Analysis CC completed....done")
    ResPat <- EnrichmentAnalysis(TFname,RegulonList,listEA_pathways,
                                 EAGenes,GOtype = "Pathway")
    print("Pathway Enrichment Analysis completed....done")

    ans <- list(ResBP = ResBP, ResMF = ResMF, ResCC = ResCC, ResPat = ResPat)
    return(ans)
}

#' @title EnrichmentAnalysis
#' @description
#'   EnrichmentAnalysis
#' @param GeneName GeneName
#' @param TableEnrichment TableEnrichment
#' @param RegulonList RegulonList
#' @param GOtype GOtype
#' @param FDRThresh FDRThresh
#' @param EAGenes EAGenes
#' @export
#' @return EAcomplete plot
#' @examples
#' \dontrun{
#' EAGenes <- get("EAGenes")
#' DAVID_BP_matrix <- get("DAVID_BP_matrix")
#' ResBP <- EnrichmentAnalysis(TFname,RegulonList,DAVID_BP_matrix,
#'                            EAGenes,GOtype = "DavidBP")
#' }
EnrichmentAnalysis <- function(GeneName,RegulonList,TableEnrichment,
                               EAGenes,GOtype,FDRThresh=0.01) {
    topPathways <- nrow(TableEnrichment)
    topPathways_tab <- matrix(0,1,topPathways)
    topPathways_tab <- as.matrix(topPathways_tab)
    rownames(topPathways_tab) <- GeneName

    rownames(EAGenes) <- toupper(rownames(EAGenes) )
    EAGenes <- EAGenes[!duplicated(EAGenes[,"ID"]),]
    rownames(EAGenes) <- EAGenes[,"ID"]
    allgene <- EAGenes[,"ID"]
    current_pathway_from_EA <- as.matrix(TableEnrichment[,GOtype]) # genes from EA pathways

    TableNames <- gsub("David","",paste("Top ", GOtype, " n. ", 1:topPathways,
                                        " of ", topPathways, sep = ""))
    colnames(topPathways_tab) <- TableNames
    topPathways_tab <- as.data.frame(topPathways_tab)

    table_pathway_enriched <- matrix(1, nrow(current_pathway_from_EA),7)
    colnames(table_pathway_enriched) <- c("Pathway","GenesInPathway","Pvalue",
                                          "FDR","CommonGenesPathway",
                                          "PercentPathway","PercentRegulon")
    table_pathway_enriched <- as.data.frame(table_pathway_enriched)

    for (i in 1:nrow(current_pathway_from_EA)) {
        table_pathway_enriched[i,"Pathway"] <- as.character(current_pathway_from_EA[i,])

        if (nrow(TableEnrichment) == 589) {
            genes_from_current_pathway_from_EA <- GeneSplitRegulon(TableEnrichment[ TableEnrichment[GOtype] == as.character(current_pathway_from_EA[i,]) ,][,"Molecules"], ",")
        }
        else {
            genes_from_current_pathway_from_EA <- GeneSplitRegulon(TableEnrichment[ TableEnrichment[GOtype] == as.character(current_pathway_from_EA[i,]) ,][,"Molecules"], ", ")
        }

        genes_common_pathway_TFregulon <- as.matrix(intersect(toupper(RegulonList),toupper(genes_from_current_pathway_from_EA)))



        if (length(genes_common_pathway_TFregulon) != 0) {
            current_pathway_commongenes_num <- length(genes_common_pathway_TFregulon)
            seta <-  allgene %in% RegulonList
            setb <-  allgene %in% genes_from_current_pathway_from_EA
            ft <- fisher.test(seta,setb)
            FisherpvalueTF <- ft$p.value
            table_pathway_enriched[i,"Pvalue"] <- as.numeric(FisherpvalueTF)
            if (FisherpvalueTF < 0.01) {
                current_pathway_commongenes_percent <- paste("(",format( (current_pathway_commongenes_num/length(genes_from_current_pathway_from_EA)) * 100,digits = 2),"%)")
                current_pathway_commongenes_num_with_percent <- gsub(" ","",paste(current_pathway_commongenes_num, current_pathway_commongenes_percent,"pv=",format(FisherpvalueTF,digits=2)))
                table_pathway_enriched[i,"CommonGenesPathway"] <- length(genes_common_pathway_TFregulon)
                table_pathway_enriched[i,"GenesInPathway"] <- length(genes_from_current_pathway_from_EA)
                table_pathway_enriched[i,"PercentPathway"] <-  as.numeric(table_pathway_enriched[i,"CommonGenesPathway"]) / as.numeric(table_pathway_enriched[i,"GenesInPathway"])  *100
                table_pathway_enriched[i,"PercentRegulon"] <-  as.numeric(table_pathway_enriched[i,"CommonGenesPathway"]) / length(RegulonList)  *100
            } }
    }
    table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[,"Pvalue"],decreasing = FALSE),]
    table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"Pvalue"] < 0.01 ,]
    table_pathway_enriched[,"FDR"] <- p.adjust(table_pathway_enriched[,"Pvalue"],method = "fdr")
    table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"FDR"] < FDRThresh ,]
    table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[,"FDR"],decreasing = FALSE),]

    tmp <- table_pathway_enriched[1:topPathways,]
    tmp <- paste(tmp[,"Pathway"],"; FDR= ", format(tmp[,"FDR"],digits = 3),"; (ng="   ,round(tmp[,"GenesInPathway"]),"); (ncommon=", format(tmp[,"CommonGenesPathway"],digits = 2), ")" ,sep = "")
    tmp <- as.matrix(tmp)
    topPathways_tab[1,] <- tmp
    rm(tmp)

    return(topPathways_tab)
}

#' @title GeneSplitRegulon
#' @description
#'   GeneSplitRegulon
#' @param Genelist Genelist
#' @param Sep Sep
#' @export
#' @return GeneSplitRegulon
#' @examples
#' GeneSplitRegulon("CRKL;TADA2A;KRT76",Sep =";")
GeneSplitRegulon <- function(Genelist,Sep){
    RegSplitted <- as.matrix(unlist(strsplit(as.character(Genelist), Sep)))

    return(RegSplitted)
}

#' @title EAbarplot
#' @description
#'   EAbarplot
#' @param tf tf
#' @param GOBPTab GOBPTab
#' @param GOCCTab GOCCTab
#' @param GOMFTab GOMFTab
#' @param PathTab PathTab
#' @param nBar nBar
#' @param nRGTab nRGTab
#' @export
#' @importFrom EDASeq barplot
#' @return EAbarplot
#' @examples
#' \dontrun{
#' Genelist <- rownames(dataDEGsFiltLevel)
#' TimeUse(ansEA <- EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
#' Enrichment Analysis EA (TCGAVisualize)
#' Gene Ontology (GO) and Pathway enrichment barPlot
#' EAbarplot(tf = rownames(ansEA$ResBP),
#'          GOBPTab = ansEA$ResBP,
#'          GOCCTab = ansEA$ResCC,
#'          GOMFTab = ansEA$ResMF,
#'         PathTab = ansEA$ResPat,
#'          nRGTab = Genelist,
#'          nBar = 10)
#'}
EAbarplot <- function(tf, GOMFTab, GOBPTab, GOCCTab, PathTab, nBar, nRGTab){
    splitFun <- function(tf, Tab, nBar){
        tmp <- lapply(Tab[tf, ], function(x) strsplit(x, ";"))
        names(tmp) <- NULL
        tmp <- matrix(unlist(tmp), ncol = 4, byrow = TRUE)
        if (nrow(tmp) == 0 | tmp[1, 1] == "NA") return(matrix(0, ncol = 2))
        tmp <- tmp[tmp[, 1] != "NA", , drop = FALSE]
        tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
        tmp[, 2] <- as.numeric(sub(" FDR= ", "", tmp[, 2]))
        tmp[, 3] <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(tmp[, 3], "=")), nrow = 2)[2, ], ")")))
        tmp[, 4] <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(tmp[, 4], "=")), nrow = 2)[2, ], ")")))

        if (nrow(tmp) < nBar) nBar <- nrow(tmp)

        tmp[, 2] <- -log10(tmp[, 2])
        o <- order(tmp[, 2], decreasing = TRUE)
        toPlot <- tmp[o[nBar:1], 1:2]
        toPlot[, 1] <- paste(toPlot[, 1], " (n=", tmp[o[nBar:1], 4], ")", sep = "")
        toPlot[, 3] <- tmp[o[nBar:1], 4]/tmp[o[nBar:1], 3]

        return(toPlot)
    }

    par(mfrow = c(2, 2))

    toPlot <- splitFun(tf, GOBPTab, nBar)
    xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = "orange",
                     main = "GO:Biological Process", xlab = "-log10(FDR)")
    labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
    text(x = 1, y = xAxis, labs, pos = 4)
    lines(x = toPlot[, 3], y = xAxis, col = "red")
    points(x = toPlot[, 3], y = xAxis, col = "red")
    axis(side = 3, at = pretty(range(0:1)), col = "red")

    toPlot <- splitFun(tf, GOCCTab, nBar)
    xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = "cyan",
                     main = "GO:Cellular Component", xlab = "-log10(FDR)")
    labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
    text(x = 1, y = xAxis, labs, pos = 4)
    lines(x = toPlot[, 3], y = xAxis, col = "red")
    points(x = toPlot[, 3], y = xAxis, col = "red")
    axis(side = 3, at = pretty(range(0:1)), col = "red")

    toPlot <- splitFun(tf, GOMFTab, nBar)
    xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = "green",
                     main = "GO:Molecular Function", xlab = "-log10(FDR)")
    labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
    text(x = 1, y = xAxis, labs, pos = 4)
    lines(x = toPlot[, 3], y = xAxis, col = "red")
    points(x = toPlot[, 3], y = xAxis, col = "red")
    axis(side = 3, at = pretty(range(0:1)), col = "red")

    toPlot <- splitFun(tf, PathTab, nBar)
    xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = "yellow",
                     main = "Pathways", xlab = "-log10(FDR)")
    labs <- toPlot[, 1]
    text(x = 1, y = xAxis, labs, pos = 4)
    lines(x = toPlot[, 3], y = xAxis, col = "red")
    points(x = toPlot[, 3], y = xAxis, col = "red")
    #axis(side = 1, at = pretty(range(0:1)), col = "red", line = 2.5)
    axis(side = 3, at = pretty(range(0:1)), col = "red")

    #par(new = TRUE)
    #plot(toPlot[, 3], xAxis, axes = FALSE, bty = "n", xlab = "",
    # ylab = "", col = "blue")
    #par(new = TRUE)
    #plot(toPlot[, 3], xAxis, type = "l", axes = FALSE, bty = "n", xlab = "",
    # ylab = "", col = "blue")
    #axis(side = 2, at = pretty(range(xAxis)))
    #axis(side = 1, at = pretty(range(toPlot[, 3])), col = "red", line=2.5)
    #axis(side = 3, at = pretty(range(toPlot[, 3])), col = "red")

    if ( is.character( nRGTab)) {
        nRG <- length(nRGTab)
    } else {
        nRG <- nRGTab[tf, "RegSizeTF"]
    }

    mainLab <- paste(tf, " (nRG = ", nRG, ")", sep = "")
    mtext(mainLab, side = 3, line = -1, outer = TRUE, font = 2)
}

