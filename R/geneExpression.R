#' @title GenesCutID
#' @description
#'   GenesCutID
#' @param GeneList GeneList
# @export
#' @return list of gene symbol without IDs
# @examples
# GenesCutID(c("CRKL|1399","TADA2A|6871","KRT76|51350"))
#' @keywords internal
GenesCutID <- function(GeneList){
    GeneListCutID <- as.matrix(matrix(unlist(strsplit(as.character(GeneList),
                                                      "|",fixed = TRUE)),length(GeneList),2,byrow = TRUE))[,1]
    return(as.matrix(GeneListCutID))
}

#' @title Filtering mRNA transcripts and miRNA selecting a threshold.
#' @description
#'    TCGAanalyze_Filtering allows user to filter mRNA transcripts and miRNA,
#'    selecting a threshold. For istance returns all mRNA or miRNA with mean across all
#'    samples, higher than the threshold defined quantile mean across all samples.
#' @param TableRnaseq is a dataframe or numeric matrix, each row represents a gene,
#' each column represents a sample come from TCGAPrepare
#' @param QuantileThresh is threshold selected as mean for filtering
#' @export
#' @return A filtered dataframe or numeric matrix where each row represents a gene,
#' each column represents a sample
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
TCGAanalyze_Filtering <- function(TableRnaseq,QuantileThresh ){
    GeneThresh <- as.numeric(quantile(rowMeans(TableRnaseq), QuantileThresh))
    geneFiltered <- names(which(rowMeans(TableRnaseq) > GeneThresh))
    Table_Rnaseq_Rawcount_Filt <- TableRnaseq[geneFiltered, ]
    return( Table_Rnaseq_Rawcount_Filt)
}

#' @title normalization mRNA transcripts and miRNA using EDASeq package.
#' @description
#'   TCGAanalyze_Normalization allows user to normalize mRNA transcripts and miRNA,
#'    using EDASeq package. Normalization for RNA-Seq Numerical and graphical summaries of RNA-Seq read data. Within-lane normalization procedures
#'    to adjust for GC-content effect (or other gene-level effects) on read counts:
#'    loess robust local regression, global-scaling, and full-quantile normalization
#'    (Risso et al., 2011). Between-lane normalization procedures to adjust for
#'    distributional differences between lanes (e.g., sequencing depth): global-scaling and full-quantile normalization (Bullard et al., 2010).
#'    For istance returns all mRNA or miRNA with mean across all
#'    samples, higher than the threshold defined quantile mean across all samples.
#'    TCGAanalyze_Normalization performs normalization using following functions from EDASeq
#'    1. EDASeq::newSeqExpressionSet
#'    2. EDASeq::withinLaneNormalization
#'    3. EDASeq::betweenLaneNormalization
#'    4. EDASeq::counts
#' @param TCGA_RnaseqTable Rnaseq numeric matrix, each row represents a gene,
#' each column represents a sample
#' @param geneInfo Information matrix of 20531 genes about geneLength and gcContent
#' @importFrom EDASeq newSeqExpressionSet withinLaneNormalization
#'  betweenLaneNormalization exprs counts
#' @export
#' @return Rnaseq matrix normalized with counts slot holds the count data as a matrix
#' of non-negative integer count values, one row for each observational unit (gene or the like),
#' and one column for each sample.
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
TCGAanalyze_Normalization <- function(TCGA_RnaseqTable,geneInfo){

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
    geneInfo <- geneInfo[commonGenes,]

    timeEstimated <- format(ncol(TCGA_RnaseqTable)*nrow(TCGA_RnaseqTable)/80000,digits = 2)
    print(messageEstimation <- paste("I Need about ", timeEstimated,
                                     "seconds for this Complete Normalization Upper Quantile",
                                     " [Processing 80k elements /s]  "))

    print("Step 1 of 4: newSeqExpressionSet ...")
    system.time(TCGA_RnaseqTable_norm <- EDASeq::newSeqExpressionSet(TCGA_RnaseqTable, featureData = geneInfo))
    print("Step 2 of 4: withinLaneNormalization ...")
    system.time(TCGA_RnaseqTable_norm <- EDASeq::withinLaneNormalization(TCGA_RnaseqTable_norm, "geneLength", which = "upper", offset = FALSE))
    print("Step 3 of 4: betweenLaneNormalization ...")
    system.time(TCGA_RnaseqTable_norm <- EDASeq::betweenLaneNormalization(TCGA_RnaseqTable_norm, which = "upper", offset = FALSE))
    print("Step 4 of 4: exprs ...")

    #system.time(TCGA_RnaseqTable_norm <- EDASeq::exprs(TCGA_RnaseqTable_norm))
    system.time(TCGA_RnaseqTable_norm <- EDASeq::counts(TCGA_RnaseqTable_norm))

    return(TCGA_RnaseqTable_norm)
}

#' @title Differentially expression analysis (DEA) using edgeR package.
#' @description
#'    TCGAanalyze_DEA allows user to perform Differentially expression analysis (DEA),
#'    using edgeR package. DEA to identify differentially expressed genes (DEGs).
#'     It is possible to do a two-class analysis.
#'     TCGAanalyze_DEA performs DEA using following functions from edgeR
#'     1. edgeR::DGEList converts the count matrix into an edgeR object.
#'     2. edgeR::estimateCommonDisp each gene gets assigned the same dispersion estimate.
#'     3. edgeR::exactTest performs pair-wise tests for differential expression between two groups.
#'     4. edgeR::topTags takes the output from exactTest(), adjusts the raw p-values using the
#'     False Discovery Rate (FDR) correction, and returns the top differentially expressed genes.
#' @param mat1 numeric matrix, each row represents a gene,
#' each column represents a sample with Cond1type
#' @param mat2 numeric matrix, each row represents a gene,
#' each column represents a sample with Cond2type
#' @param Cond1type a string containing the class label of the samples in mat1
#'  (e.g., control group)
#' @param Cond2type a string containing the class label of the samples in mat2
#' (e.g., case group)
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags
#' @export
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
#' samplesNT <- MultiSampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- MultiSampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- TCGAanalyze_DEA(dataFilt[,samplesNT],
#'                       dataFilt[,samplesTP],"Normal", "Tumor")
#' @return table with DEGs containing for each gene logFC, logCPM, pValue,and FDR
TCGAanalyze_DEA <- function(mat1,mat2,Cond1type,Cond2type) {

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

#' @title Adding information related to DEGs genes from DEA as mean values in two conditions.
#' @description
#'    TCGAanalyze_LevelTab allows user to add information related to DEGs genes from
#'    Differentially expression analysis (DEA) such as mean values and in two conditions.
#' @param FC_FDR_table_mRNA Output of dataDEGs filter by abs(LogFC) >=1
#' @param typeCond1 a string containing the class label of the samples
#'  in TableCond1  (e.g., control group)
#' @param typeCond2 a string containing the class label of the samples
#' in TableCond2  (e.g., case group)
#' @param TableCond1 numeric matrix, each row represents a gene, each column
#'  represents a sample with Cond1type
#' @param TableCond2 numeric matrix, each row represents a gene, each column
#' represents a sample with Cond2type
#' @param typeOrder typeOrder
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags
#' @export
#' @return table with DEGs, log Fold Change (FC), false discovery rate (FDR),
#' the gene expression level
#' for samples in  Cond1type, and Cond2type, and Delta value (the difference
#' of gene expression between the two
#' conditions multiplied logFC)
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
#' samplesNT <- MultiSampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- MultiSampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- TCGAanalyze_DEA(dataFilt[,samplesNT], dataFilt[,samplesTP],
#' "Normal", "Tumor")
#' dataDEGsFilt <- dataDEGs[abs(dataDEGs$logFC) >= 1,]
#' dataTP <- dataFilt[,samplesTP]
#' dataTN <- dataFilt[,samplesNT]
#' dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGsFilt,"Tumor","Normal",
#' dataTP,dataTN)
TCGAanalyze_LevelTab <- function(FC_FDR_table_mRNA,typeCond1,typeCond2,
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

#' @title Principal components analysis (PCA) plot
#' @description
#'   TCGAvisualize_PCA performs a principal components analysis (PCA) on the given data matrix
#'   and returns the results as an object of class prcomp, and shows results in PCA level.
#' @param dataFilt A filtered dataframe or numeric matrix where each row represents a gene,
#' each column represents a sample from function TCGAanalyze_Filtering
#' @param dataDEGsFiltLevel table with DEGs, log Fold Change (FC), false discovery rate (FDR),
#' the gene expression level, etc, from function TCGAanalyze_LevelTab.
#' @param ntopgenes number of DEGs genes to plot in PCA
#' @import ggplot2
#' @export
#' @return principal components analysis (PCA) plot of PC1 and PC2
#' @examples
#' # normalization of genes
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' # quantile filter of genes
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
#' # Principal Component Analysis plot for ntop selected DEGs
#' TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel, ntopgenes = 200)
#'
TCGAvisualize_PCA <- function(dataFilt,dataDEGsFiltLevel ,ntopgenes) {
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


    g <- ggbiplot(cancer.pca, obs.scale = 1, var.scale = 1,
                  groups = sampleColors, ellipse = TRUE, circle = FALSE)
    g <- g + scale_colour_manual(name = "",
                                 values = c("blue" = "blue","red" = "red"))
    with(g,
         g <- g + geom_point(aes(colour = sampleColors), size = 3)
    )
    #shape = tabClusterNew$Study)
    g <- g + theme(legend.direction = 'horizontal',  legend.position = 'top')
    g <- g + ggtitle(TitlePlot)
    print(g)
    return(cancer.pca)
}

#' @title Enrichment analysis for Gene Ontology (GO) [BP,MF,CC] and Pathways
#' @description
#'   Researchers, in order to better understand the underlying biological
#'   processes, often want to retrieve a functional profile of a set of genes
#'   that might have an important role. This can be done by performing an
#'   enrichment analysis.
#'
#'We will perform an enrichment analysis on gene sets using the TCGAanalyze_EAcomplete
#'function. Given a set of genes that are
#'up-regulated under certain conditions, an enrichment analysis will find
#'identify classes of genes or proteins that are #'over-represented using
#'annotations for that gene set.
#' @param TFname is the name of the list of genes or TF's regulon.
#' @param RegulonList List of genes such as TF's regulon or DEGs where to find enrichment.
#' @export
#' @return Enrichment analysis GO[BP,MF,CC] and Pathways complete table enriched by genelist.
#' @examples
#' Genelist <- c("FN1","COL1A1")
#' ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
#' \dontrun{
#' Genelist <- rownames(dataDEGsFiltLevel)
#' system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
#' }
TCGAanalyze_EAcomplete <- function(TFname, RegulonList){

    print(paste("I need about ", "1 minute to finish complete ",
                "Enrichment analysis GO[BP,MF,CC] and Pathways... "))

    ResBP <- TCGAanalyze_EA(TFname,RegulonList,DAVID_BP_matrix,
                                EAGenes,GOtype = "DavidBP")
    print("GO Enrichment Analysis BP completed....done")
    ResMF <- TCGAanalyze_EA(TFname,RegulonList,DAVID_MF_matrix,
                                EAGenes,GOtype = "DavidMF")
    print("GO Enrichment Analysis MF completed....done")
    ResCC <- TCGAanalyze_EA(TFname,RegulonList,DAVID_CC_matrix,
                                EAGenes,GOtype = "DavidCC")
    print("GO Enrichment Analysis CC completed....done")
    ResPat <- TCGAanalyze_EA(TFname,RegulonList,listEA_pathways,
                                 EAGenes,GOtype = "Pathway")
    print("Pathway Enrichment Analysis completed....done")

    ans <- list(ResBP = ResBP, ResMF = ResMF, ResCC = ResCC, ResPat = ResPat)
    return(ans)
}

#' @title Enrichment analysis of a gene-set with GO [BP,MF,CC]  and pathways.
#' @description
#' The rational behind a enrichment analysis ( gene-set, pathway etc) is to compute
#' statistics of whether the overlap between the focus list (signature) and the gene-set
#' is significant. ie the confidence that overlap between the list is not due to chance.
#'  The Gene Ontology project describes genes (gene products) using terms from
#'  three structured vocabularies: biological process, cellular component and molecular function.
#'  The Gene Ontology Enrichment component, also referred to as the GO Terms" component, allows
#'  the genes in any such "changed-gene" list to be characterized using the Gene Ontology terms
#'  annotated to them. It asks, whether for any particular GO term, the fraction of genes
#'  assigned to it in the "changed-gene" list is higher than expected by chance
#'  (is over-represented), relative to the fraction of genes assigned to that term in the
#'  reference set.
#'  In statistical terms it peform the analysis tests the null hypothesis that,
#'  for any particular ontology term, there is no diffeerence in the proportion of genes
#'  annotated to it in the reference list and the proportion annotated to it in the test list.
#'  We adopted a Fisher Exact Test to perform the EA.
#' @param GeneName is the name of gene signatures list
#' @param TableEnrichment is a table related to annotations of gene symbols such as
#' GO[BP,MF,CC] and Pathways. It was created from DAVID gene ontology on-line.
#' @param RegulonList is a gene signature (lisf of genes) in which perform EA.
#' @param GOtype is type of gene ontology Biological process (BP), Molecular Function (MF),
#' Cellular componet (CC)
#' @param FDRThresh pvalue corrected (FDR) as threshold to selected significant
#' BP, MF,CC, or pathways. (default FDR < 0.01)
#' @param EAGenes is a table with informations about genes
#' such as ID, Gene, Description, Location and Family.
# @export
#' @import stats
#' @return Table with enriched GO or pathways by selected gene signature.
#' @examples
#' \dontrun{
#' EAGenes <- get("EAGenes")
#' RegulonList <- rownames(dataDEGsFiltLevel)
#' ResBP <- TCGAanalyze_EA(GeneName="DEA genes Normal Vs Tumor",
#'                            RegulonList,DAVID_BP_matrix,
#'                            EAGenes,GOtype = "DavidBP")
#'}
TCGAanalyze_EA <- function(GeneName,RegulonList,TableEnrichment,
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
# @export
#' @return GeneSplitRegulon
# @examples
# GeneSplitRegulon("CRKL;TADA2A;KRT76",Sep =";")
GeneSplitRegulon <- function(Genelist,Sep){
    RegSplitted <- as.matrix(unlist(strsplit(as.character(Genelist), Sep)))

    return(RegSplitted)
}

#' @title barPlot for a complete Enrichment Analysis
#' @description
#'   TCGAvisualize_EAbarplot plots the result from TCGAanalyze_EAcomplete in a complete barPlot
#' @param tf is a list of gene symbols
#' @param GOBPTab is results from TCGAanalyze_EAcomplete related to Biological Process (BP)
#' @param GOCCTab is results from TCGAanalyze_EAcomplete related to Cellular Component (CC)
#' @param GOMFTab is results from TCGAanalyze_EAcomplete related to Molecular Function (MF)
#' @param PathTab is results from TCGAanalyze_EAcomplete related to Pathways EA
#' @param nBar is the number of bar histogram selected to show (default = 10)
#' @param nRGTab is the gene signature list with gene symbols.
#' @export
#' @importFrom EDASeq barplot
#' @import graphics
#' @return Complete barPlot from Enrichment Analysis showing significant (default FDR < 0.01)
#' BP,CC,MF and pathways enriched by list of genes.
#' @examples
#' Genelist <- c("FN1","COL1A1")
#' ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
#' TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
#'          GOBPTab = ansEA$ResBP,
#'          GOCCTab = ansEA$ResCC,
#'          GOMFTab = ansEA$ResMF,
#'         PathTab = ansEA$ResPat,
#'          nRGTab = Genelist,
#'          nBar = 10)
#' \dontrun{
#' Genelist <- rownames(dataDEGsFiltLevel)
#' system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
#' # Enrichment Analysis EA (TCGAVisualize)
#' # Gene Ontology (GO) and Pathway enrichment barPlot
#' TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
#'          GOBPTab = ansEA$ResBP,
#'          GOCCTab = ansEA$ResCC,
#'          GOMFTab = ansEA$ResMF,
#'         PathTab = ansEA$ResPat,
#'          nRGTab = Genelist,
#'          nBar = 10)
#'}
TCGAvisualize_EAbarplot <- function(tf, GOMFTab, GOBPTab, GOCCTab, PathTab, nBar, nRGTab){
    splitFun <- function(tf, Tab, nBar){
        tmp <- lapply(Tab[tf, ], function(x) strsplit(x, ";"))
        names(tmp) <- NULL
        tmp <- matrix(unlist(tmp), ncol = 4, byrow = TRUE)
        if (nrow(tmp) == 0 | tmp[1, 1] == "NA") return(matrix(0, ncol = 2))
        tmp <- tmp[tmp[, 1] != "NA", , drop = FALSE]
        tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
        tmp[, 2] <- as.numeric(sub(" FDR= ", "", tmp[, 2]))
        tmp[, 3] <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(tmp[, 3],
                                                                      "=")), nrow = 2)[2, ], ")")))
        tmp[, 4] <- as.numeric(unlist(strsplit(matrix(unlist(strsplit(tmp[, 4],
                                                                      "=")), nrow = 2)[2, ], ")")))

        if (nrow(tmp) < nBar) nBar <- nrow(tmp)

        tmp[, 2] <- -log10(tmp[, 2])
        o <- order(tmp[, 2], decreasing = TRUE)
        toPlot <- tmp[o[nBar:1], 1:2]
        toPlot[, 1] <- paste0(toPlot[, 1], " (n=", tmp[o[nBar:1], 4], ")")
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

