#' @title Survival analysis with univariate Cox regression package (dnet)
#' @description TCGAvisualize_SurvivalCoxNET can help an user to identify a group of survival genes that are
#' significant from univariate Kaplan Meier Analysis and also for Cox Regression.
#' It shows in the end a network build with community of genes with similar range of pvalues from
#' Cox regression (same color) and that interaction among those genes is already validated in
#' literatures using the STRING database (version 9.1).
#' TCGAvisualize_SurvivalCoxNET perform survival analysis with univariate Cox regression
#' and package (dnet) using following functions wrapping from these packages:
#' \enumerate{
#' \item survival::coxph
#' \item igraph::subgraph.edges
#' \item igraph::layout.fruchterman.reingold
#' \item igraph::spinglass.community
#' \item igraph::communities
#' \item dnet::dRDataLoader
#' \item dnet::dNetInduce
#' \item dnet::dNetPipeline
#' \item dnet::visNet
#' \item dnet::dCommSignif
#' }
#' @importFrom igraph subgraph.edges layout.fruchterman.reingold
#'             spinglass.community degree E communities crossing V V<-
#' @importFrom dnet dRDataLoader dNetInduce dNetPipeline
#'             visNet dCommSignif
#' @importFrom supraHex visColormap visColoralpha
#' @importFrom grDevices dev.list
#' @details TCGAvisualize_SurvivalCoxNET allow user to perform the complete workflow using coxph
#' and dnet package related to survival analysis with an identification of gene-active networks from
#' high-throughput omics data using gene expression and clinical data.
#' \enumerate{
#' \item Cox regression survival analysis to obtain hazard ratio (HR) and pvaules
#' \item fit a Cox proportional hazards model and ANOVA (Chisq test)
#' \item Network comunites
#' \item An igraph object that contains a functional protein association network in human.
#' The network is extracted from the STRING database (version 9.1).
#' Only those associations with medium confidence (score>=400) are retained.
#' \item restrict to those edges with high confidence (score>=700)
#' \item extract network that only contains genes in pvals
#' \item Identification of gene-active network
#' \item visualisation of the gene-active network itself
#' \item the layout of the network visualisation (fixed in different visuals)
#' \item color nodes according to communities (identified via a spin-glass model and simulated annealing)
#' \item node sizes according to degrees
#' \item highlight different communities
#' \item visualise the subnetwork
#' }
#' @param clinical_patient is a data.frame using function 'clinic' with information
#' related to barcode / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_followup , vital_status, etc
#' @param dataGE is a matrix of Gene expression (genes in rows, samples in cols) from TCGAprepare
#' @param Genelist is a list of gene symbols where perform survival KM.
#' @param org.Hs.string an igraph object that contains a functional protein association network
#' in human. The network is extracted from the STRING database (version 10).
#' @param scoreConfidence restrict to those edges with high confidence (eg. score>=700)
#' @param titlePlot is the title to show in the final plot.
#' @importFrom survival coxph
#' @importFrom igraph subgraph.edges layout.fruchterman.reingold
#'             spinglass.community degree E communities crossing V V<-
#' @importFrom dnet dRDataLoader dNetInduce dNetPipeline
#'             visNet dCommSignif
#' @importFrom supraHex visColormap visColoralpha
#' @importFrom grDevices dev.list
#' @export
#' @return net IGRAPH with related Cox survival genes in community (same pval and color) and with
#' interactions from STRING database.
#' query <- TCGAquery(tumor = "lgg")
TCGAvisualize_SurvivalCoxNET <- function(clinical_patient,
                                         dataGE,
                                         Genelist,
                                         org.Hs.string,
                                         scoreConfidence = 700,
                                         titlePlot = "TCGAvisualize_SurvivalCoxNET Example"){

    #clinical_patient<- dataClin
    #dataGE <- dataFilt
    #Genelist <- rownames(dataSurv)
    #scoreConfidence = 700

    combined_score <- NULL
    if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}

    pdf("SurvivalCoxNETOutput.pdf", width = 15, height = 10)

    ## fit a Cox proportional hazards model for age, gender, tumor type
    cfu<-clinical_patient[clinical_patient[,"bcr_patient_barcode"] %in% substr(colnames(dataGE),1,12),]
    rownames(cfu)<- cfu$bcr_patient_barcode
    cfu <- as.data.frame(subset(cfu, select=c("bcr_patient_barcode",
                                              "days_to_last_followup",
                                              "days_to_death",
                                              "vital_status",
                                              "age_at_initial_pathologic_diagnosis",
                                              "gender")
    )
    )

    rownames(cfu)<- cfu$bcr_patient_barcode
    cfu<-cfu[,-1]
    colnames(cfu)   <- c("time","timeDead","status","Age","Gender")

    cfu[which(cfu$status=="Alive"),"status"]<-0
    cfu[which(cfu$status=="Dead"),"time"]<- cfu[which(cfu$status=="Dead"),"timeDead"]
    cfu[which(cfu$status=="Dead"),"status"]<-1
    cfu$Gender <-tolower(cfu$Gender)
    cfu<-cfu[,-2]
    cfu <- as.data.frame(subset(cfu, select=c("time","status")))
    cfu$time <- as.numeric(cfu$time)
    cfu$status <- as.numeric(cfu$status)

    data <-as.data.frame(dataGE)
    colnames(data) <- substr(colnames(data),1,12)
    commonSamples <- intersect(rownames(cfu), colnames(data))

    data <- t(data)
    data <- data[commonSamples,Genelist]

    #rownames(tabSurvKMfilt) <- tabSurvKMfilt$mRNA
    #colnames(Cancer_rnaseqv2) <- substr(colnames(Cancer_rnaseqv2),1,12)
    #md_selected<-log2(Cancer_rnaseqv2[rownames(tabSurvKMfilt),rownames(cfu)])

    md_selected<- data
    pd <- cfu

    ## survival analysis to obtain hazard ratio (HR) and pvaules
    HR <- rep(1, ncol(md_selected))
    pvals <- rep(1, ncol(md_selected))
    for(i in 1:ncol(md_selected)){
        ## fit a Cox proportional hazards model
        data <- cbind(pd, gene = md_selected[,i])
        data <- data [which(data$gene !="-Inf"),]

        fit <- coxph(formula=Surv(time,status) ~., data=data)
        ## ANOVA (Chisq test)
        cox <- summary(fit)
        cat(paste( (ncol(md_selected)-i),".",sep=""))
        # res <- as.matrix(anova(fit, test="Chisq"))
        HR[i] <- as.numeric(cox$coefficients[2])
        pvals[i] <- as.numeric(cox$logtest[3])
    }

    names(HR) <- colnames(md_selected)
    names(pvals) <- colnames(md_selected)
    # Network comunites >>=

    # An igraph object that contains a functional protein association network
    # in human. The network is extracted from the STRING database (version 9.1).
    # Only those associations with medium confidence (score>=400) are retained.
    #  org.Hs.string <- dRDataLoader(RData='org.Hs.string')
    # restrict to those edges with high confidence (score>=700)
    # with(org.Hs.string,{
    #    network <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=scoreConfidence])})
    #network
    network <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=scoreConfidence])


    # extract network that only contains genes in pvals
    ind <- match(V(network)$symbol, names(pvals))
    ## for extracted graph
    nodes_mapped <- V(network)$name[!is.na(ind)]
    network <- dNetInduce(g=network, nodes_query=nodes_mapped, knn=0,
                          remove.loops=FALSE, largest.comp=TRUE)
    V(network)$name <- V(network)$symbol

    # Identification of gene-active network
    net <- dNetPipeline(g=network, pval=pvals, method="customised",
                        significance.threshold=5e-02)
    # visualisation of the gene-active network itself
    ## the layout of the network visualisation (fixed in different visuals)
    glayout <- layout.fruchterman.reingold(net)
    ## color nodes according to communities (identified via a spin-glass model and simulated annealing)
    com <- spinglass.community(net, spins=25)
    com$csize <- sapply(1:length(com),function(x) sum(com$membership==x))
    vgroups <- com$membership
    colormap <- "yellow-darkorange"
    palette.name <- visColormap(colormap=colormap)
    mcolors <- palette.name(length(com))
    vcolors <- mcolors[vgroups]
    com$significance <- dCommSignif(net, com)
    ## node sizes according to degrees
    vdegrees <- degree(net)
    ## highlight different communities
    mark.groups <- communities(com)
    mark.col <- visColoralpha(mcolors, alpha=0.2)
    mark.border <- visColoralpha(mcolors, alpha=0.2)
    edge.color <- c("#C0C0C0", "#000000")[crossing(com,net)+1]
    edge.color <- visColoralpha(edge.color, alpha=0.5)
    ## visualise the subnetwrok
    visNet(g=net, glayout=glayout, vertex.label=V(net)$geneSymbol,
           vertex.color=vcolors, vertex.frame.color=vcolors,
           vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col,
           mark.border=mark.border, mark.shape=1, mark.expand=10,
           edge.color=edge.color, newpage=FALSE, vertex.label.color="blue",
           vertex.label.dist=0.4, vertex.label.font=2, main = titlePlot)
    legend_name <- paste("C",1:length(mcolors)," (n=",com$csize,", pval=",signif(com$significance,digits=2),")",sep='')
    legend("topleft", legend=legend_name, fill=mcolors, bty="n", cex=1.4)

    dev.off()

    return(net)

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
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo = geneInfo,
#' method = "geneLength")
#' # quantile filter of genes
#' dataFilt <- TCGAanalyze_Filtering(tabDF = dataBRCA, method = "quantile", qnt.cut =  0.25)
#' # Principal Component Analysis plot for ntop selected DEGs
#' TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel, ntopgenes = 200)
TCGAvisualize_PCA <- function(dataFilt,dataDEGsFiltLevel ,ntopgenes) {
    ComparisonSelected <- "Normal vs Tumor"
    TitlePlot <- paste0("PCA ", "top ", ntopgenes,
                        " Up and down diff.expr genes between ",
                        ComparisonSelected)


    dataFilt <- dataFilt[!duplicated(GenesCutID(rownames(dataFilt))),]
    rownames(dataFilt) <- GenesCutID(rownames(dataFilt))

    Genelist <- rownames(dataDEGsFiltLevel)[1:ntopgenes]
    commonGenes <- intersect(Genelist, rownames(dataFilt) )
    expr2 <- dataFilt[commonGenes,]
    color1 <- "blue"
    color2 <- "red"

    # selection of normal samples "NT"
    samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt),
                                       typesample = c("NT"))
    # selection of tumor samples "TP"
    samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt),
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

    pdf("TCGAvisualize_EAbarplot_Output.pdf", width = 15, height = 15)

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

    dev.off()
}

#' @title Barplot of subtypes and clinical info in groups of gene expression clustered.
#' @description
#'   Barplot of subtypes and clinical info in groups of gene expression clustered.
#' @param DFfilt write
#' @param DFclin write
#' @param DFsubt write
#' @param data_Hc2 write
#' @param Subtype write
#' @param cbPalette Define the colors of the bar.
#' @param filename The name of the pdf file
#' @param width Image width
#' @param height Image height
#' @param dpi Image dpi
#' @import ggplot2
#' @export
#' @return barplot image in pdf or png file
#' @examples
#' query <- TCGAquery(tumor = "lgg")
TCGAvisualize_BarPlot <- function(DFfilt, DFclin, DFsubt, data_Hc2,
                                  Subtype, cbPalette, filename, width, height,dpi){

    if(Subtype =="AGE"){
        dataClinNew <- dataClin
        colnames(dataClin)[which(colnames(dataClin) == "age_at_initial_pathologic_diagnosis")] <- "AGE"
        dataClin$AGE <- as.numeric(as.character(dataClin$AGE))
        dataClin <- cbind(dataClin, AGE2 = matrix(0,nrow(dataClin),1))
        dataClin[ dataClin$AGE <= 32, "AGE2"] <- "<=32 yr"
        dataClin[ dataClin$AGE >= 41, "AGE2"] <- ">=41 yr"
        dataClin[ dataClin$AGE2 == 0, "AGE2"] <- "33-40 yr"
        dataClin$AGE <- as.character(dataClin$AGE2)
        DFclin <- dataClin
    }

    ans <- hclust(ddist <- dist(DFfilt), method = "ward.D2")
    hhc <- data_Hc2[[4]]$consensusTree
    consensusClusters<-data_Hc2[[4]]$consensusClass
    sampleOrder <- consensusClusters[hhc$order]

    consensusClusters <- as.factor(data_Hc2[[4]]$clrs[[1]])
    names(consensusClusters) <- attr(ddist, "Labels")
    names(consensusClusters) <- substr(names(consensusClusters),1,12)

    # adding information about gropus from consensus Cluster in clinical data
    DFclin <- cbind(DFclin, groupsHC = matrix(0,nrow(DFclin),1))
    rownames(DFclin) <- DFclin$bcr_patient_barcode

    for( i in 1: nrow(DFclin)){
        currSmp <- DFclin$bcr_patient_barcode[i]
        DFclin[currSmp,"groupsHC"] <- as.character(consensusClusters[currSmp])
    }

    DFclin_filt <- DFclin[DFclin$bcr_patient_barcode %in% DFsubt$patient,]
    DFclin_filt <- DFclin_filt[order(DFclin_filt$bcr_patient_barcode),]
    DFsubt <- DFsubt[order(DFsubt$patient),]
    DFclin_merged <- cbind(DFclin_filt,DFsubt)

    subtype_sel <- colnames(DFclin_merged) == Subtype
    DFclin_merged[,Subtype] <- as.character(DFclin_merged[,Subtype])
    DFclin_merged <- DFclin_merged[!(is.na(DFclin_merged[,Subtype])),]
    DFclin_merged <- DFclin_merged[DFclin_merged[,Subtype] !="NA",]

    groupsColors <-  levels(as.factor(DFclin_merged$groupsHC))
    for(j in 1:length(table(DFclin_merged$groupsHC))){
        curCol <- groupsColors[j]
        DFclin_merged[DFclin_merged$groupsHC == curCol,"groupsHC"]<-paste0("EC",j)
    }

    subtypeCluster <- factor(DFclin_merged[,Subtype])
    pplot <- qplot(factor(DFclin_merged$groupsHC),
                   data=DFclin_merged, geom="bar",
                   fill=subtypeCluster , xlab="Expression group") +
        theme_bw()  +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) +
        guides(fill=guide_legend(title=Subtype))  +
        theme(legend.position="top",
              legend.text = element_text(size = 18),
              legend.title = element_text(size=18, face="bold")) +
        theme( legend.text = element_text(size = 18),
               legend.title = element_text(size = 18),
               axis.text= element_text(size = 30),
               axis.title.x= element_text(size = 22),
               axis.title.y= element_text(size = 30))   +
        scale_fill_manual(values=cbPalette)

    ggsave(pplot, filename = filename, width = width, height = height, dpi = dpi)
    message(paste("Plot saved in: ", file.path(getwd(),filename)))
}

#' @title Visaulize results in format of latex tables.
#' @description Visaulize results in format of latex tables.
#' @param Table write
#' @param rowsForPage write
#' @param TableTitle write
#' @param LabelTitle write
#' @param withrows write
#' @param size size selected for font, 'small', 'tiny'
#' @importFrom xtable xtable
#' @importFrom gplots greenred
#' @examples
#' library(stringr)
#' tabDEGsTFPubmed$PMID <- str_sub(tabDEGsTFPubmed$PMID,0,30)
#' TCGAvisualize_Tables(Table = tabDEGsTFPubmed,
#' rowsForPage = 5,
#' TableTitle = "pip",
#' LabelTitle = "pip2",
#' withrows = FALSE,
#' size = "small")
#' @export
#' @return table in latex format to use in beamer presentation or sweave files
TCGAvisualize_Tables <- function(Table, rowsForPage, TableTitle, LabelTitle, withrows, size){
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
            tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size= size)
            print(tablePrint_Table_current,include.rownames = withrows)
        }

        else{
            print(i)
            if( i == 1 ) {
                Table_current<-Table[i:vectorFirst[i],]
                tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size=size)
                print(tablePrint_Table_current,include.rownames = withrows)
            }

            else if (i==numberOfprint) {
                Table_current<-Table[vectorLast[i-1]:nrow(Table),]
                tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size= size)
                print(tablePrint_Table_current,include.rownames = withrows)
            }
            else{
                Table_current<-Table[vectorLast[i-1]:vectorFirst[i],]
                tablePrint_Table_current<-xtable(Table_current, caption = paste(TableTitle,"(",i,")"),label = gsub(" ","",paste(LabelTitle,".",i)) , size= size)
                print(tablePrint_Table_current,include.rownames = withrows)
            }
        }
    }

}

#' @title Heatmap with more sensible behavior using heatmap.plus
#' @description Heatmap with more sensible behavior using heatmap.plus
#' @param cancer tumor selected for the analysis
#' @param DFfilt write
#' @param DFclin write
#' @param DFsubt write
#' @param data_Hc2 write
#' @param cbPalette write
#' @param filename write. default = NULL
#' @importFrom heatmap.plus heatmap.plus
#' @examples
#' query <- TCGAquery(tumor = "lgg")
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
TCGAvisualize_Heatmap <- function(cancer, DFfilt, DFclin, DFsubt, data_Hc2, cbPalette, filename =NULL){

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
    GE <- t(.quantileNormalization(t(DFfilt)))
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

    if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}

    curDate <- as.character(unlist(strsplit(gsub(" ","_h",
                                                 gsub("-","_",as.character(Sys.time()))),":"))[1])


    pdf(file=paste0(curDate,"_",cancer,"_heatmap_with_subtypes_withHeatmapPlus.pdf"))

    .heatmap.plus.sm(
        t(oGE),
        na.rm=TRUE,
        scale="none",
        #RowSideColor=probe.cc,
        #ColSideColors=cc.col,
        col=gplots::greenred(75),
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



#' @title Profile plot
#' @description Displaty the association between cancer subtypes and any kind of clustering.
#' @param data A data frame with the cluters and subytpe of cancers
#' @param subtypeCol Name of the column with the subtype information
#' @param groupCol Names of tre columns with the cluster information
#' @param filename Name of the file to save the plot, can be pdf, png, svg etc..
#' @param na.rm Remove NA groups? Default = FALSE
#' @importFrom sjPlot sjp.stackfrq
#' @examples
#' query <- TCGAquery(tumor = "lgg")
#' \dontrun{
#' clin <- TCGAquery_clinic("lgg","clinical_patient")
#' TCGAvisualize_profilePlot ()
#' }
#' @export
#' @return A plot
TCGAvisualize_profilePlot <- function (data = NULL,
                                       groupCol = NULL,
                                       subtypeCol = NULL,
                                       colors = NULL,
                                       filename = NULL,
                                       na.rm = FALSE) {

    if (is.null(groupCol)) stop("Please provide the groupCol argument")
    if (is.null(subtypeCol)) stop("Please provide the subtypeCol argument")
    if (is.null(data)) stop("Please provide the data argument")
    if (is.null(filename)) filename <- paste0(groupCol,subtypeCol,".pdf")

    if(na.rm){
        data <- data[!is.na(data[,groupCol]),]
        data <- data[which(data[,groupCol] != "NA"),]
    }

    # use https://github.com/cttobin/ggthemr
    # when it is in cran
    if(is.null(colors)) colors <- c("#34495e",
                                    "#3498db",
                                    "#2ecc71",
                                    "#f1c40f",
                                    "#e74c3c",
                                    "#9b59b6",
                                    "#1abc9c",
                                    "#f39c12",
                                    "#d35400")

    # The ideia is: we have a data frame like this:
    #----------------------
    # GROUP COL  SUBTYPE
    #---------------------
    # group 1     WT
    # group 2     NC
    # group 1     WT
    # group 3     NC
    #---------------------
    # And we need
    #-------------------
    # Group 1 Group 2 Group 3
    #   1       2       2
    #   1       NA      NA
    # ---------------------
    # where 1 is WT and NC 2

    df <- as.data.frame(data)
    df <- dcast(df, as.formula(paste0(subtypeCol, " ~ ", groupCol)))
    var.labels <- unique(df[,1]) # get the cluster names
    m <- max(apply(df[,-1],2,sum)) # get the max number of subtypes in the clusters

    # create a data frame with all values
    for(i in 2:ncol(df)){
        x <- c()
        for(j in 1:nrow(df)){
            x <- c(x,rep(j,(df[j,i])))
        }
        missing <- m-length(x)
        x <- c(x,rep(NA,missing))
        if(i == 2) data <- x
        if(i > 2) data <- cbind(data,x)
    }
    colnames(data) <- colnames(df)[-1]

    # create a collumn for all values
    all <- as.numeric(unlist(data))
    idx <- length(all) - nrow(data)
    for( i in 1:idx) {
        data <- rbind(data, rep(NA,ncol(data)))
    }

    data <- cbind(all,data)
    colnames(data)[1] <- subtypeCol
    data <- as.data.frame(data)

    p <- sjp.stackfrq(data,
                      legendTitle = subtypeCol,
                      axisTitle.x = groupCol,
                      #sort.frq = "last.desc",
                      expand.grid = FALSE,
                      legendLabels = as.character(var.labels),
                      jitterValueLabels = TRUE,
                      showSeparatorLine = TRUE,
                      showValueLabels = FALSE,
                      geom.colors = colors[1:length(var.labels)],
                      #separatorLineColor = "#6699cc"
                      printPlot = TRUE)
    ggsave(p$plot, filename = filename, width = 10, height = 10, dpi = 600)
}




#' @title Visualize mutation
#' @description See % of genes mutated
#' @param data A data frame with the cluters and subytpe of cancers
#' @param samples Samples to consider mutation
#' @param geneList List of genes to plot
#' @param filename Name of the file to save the plot, can be pdf, png, svg etc..
#' @param threshold plot only if number of mutated genes are higher than threshold
#' @importFrom sjPlot sjp.stackfrq
#' @export
#' @return A plot
TCGAvisualize_mutation <- function (data = NULL,
                                    colors = NULL,
                                    samples = NULL,
                                    geneList = NULL,
                                    filename = NULL,
                                    threshold = 0) {

    if (is.null(data)) stop("Please provide the data argument")
    if (is.null(filename)) filename <- "mutation_summary.pdf"

    # subset samples
    if(!is.null(samples)){
        idx <- which(samples == data$patient)
        data <- data[idx,]
    }


    # subset samples
    if(!is.null(geneList)){
        genes <- geneList
    } else {
        genes <- na.omit(unique(unlist(data$genes)))
    }

    df <-  data.frame(matrix(NA, ncol = length(genes), nrow = nrow(data)))
    colnames(df) <- genes
    summary <- table(unlist(data$genes))

    for( i in genes){
        print(i)
        count <- summary[i]
        count <- as.numeric(count)

        if( count < threshold){
            df[,i] <- NULL
        } else {
        print(count)
            df[,i] <- c(rep(1,as.numeric(count)), rep(2,nrow(data)-count))
        }
    }

    if(ncol(df) == 0){
        stop("No genes found")
    }
    if(ncol(df) == 1){
        df <- as.data.frame(df)
    }

    print(df)

    # use https://github.com/cttobin/ggthemr
    # when it is in cran
    if(is.null(colors)) colors <- c("#34495e",
                                    "#3498db",
                                    "#2ecc71",
                                    "#f1c40f",
                                    "#e74c3c",
                                    "#9b59b6",
                                    "#1abc9c",
                                    "#f39c12",
                                    "#d35400")

    # The ideia is: we have a data frame like this:
    #----------------------
    # GROUP COL  SUBTYPE
    #---------------------
    # group 1     WT
    # group 2     NC
    # group 1     WT
    # group 3     NC
    #---------------------
    # And we need
    #-------------------
    # Group 1 Group 2 Group 3
    #   1       2       2
    #   1       NA      NA
    # ---------------------
    # where 1 is WT and NC 2


    p <- sjp.stackfrq(df,
                      #legendTitle = subtypeCol,
                      #axisTitle.x = groupCol,
                      #sort.frq = "last.desc",
                      expand.grid = FALSE,
                      legendLabels = c("Mutated","Not mutated"),
                      jitterValueLabels = TRUE,
                      showSeparatorLine = TRUE,
                      showValueLabels = FALSE,
                      geom.colors = colors[1:2],
                      #separatorLineColor = "#6699cc"
                      printPlot = TRUE)
    ggsave(p$plot, filename = filename, width = 10, height = 10, dpi = 600)
}
