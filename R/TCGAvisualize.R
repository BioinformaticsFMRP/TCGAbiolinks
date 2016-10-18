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
#' pca <- TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel, ntopgenes = 200)
#' if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
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
#' The figure shows canonical pathways significantly overrepresented (enriched) by the DEGs
#' (differentially expressed genes).
#' The most statistically significant canonical pathways identified
#' in DEGs list are listed according to their p value corrected FDR (-Log) (colored bars)
#' and the ratio of list genes found in each pathway over the total number of
#' genes in that pathway (Ratio, red line).
#' @param tf is a list of gene symbols
#' @param GOBPTab is results from TCGAanalyze_EAcomplete related to Biological Process (BP)
#' @param GOCCTab is results from TCGAanalyze_EAcomplete related to Cellular Component (CC)
#' @param GOMFTab is results from TCGAanalyze_EAcomplete related to Molecular Function (MF)
#' @param PathTab is results from TCGAanalyze_EAcomplete related to Pathways EA
#' @param nBar is the number of bar histogram selected to show (default = 10)
#' @param nRGTab is the gene signature list with gene symbols.
#' @param filename Name for the pdf. If null it will return the plot.
#' @param color A vector of colors for each barplot. Deafult:  c("orange", "cyan","green","yellow")
#' @param text.size Text size
#' @param xlim Upper limit of the x-axis.
#' @param mfrow Vector with number of rows/columns of the plot. Default  2 rows/2 columns "c(2,2)"
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
#'          nBar = 10,
#'          filename="a.pdf")
#' while (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
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
TCGAvisualize_EAbarplot <- function(tf, GOMFTab, GOBPTab, GOCCTab, PathTab, nBar, nRGTab,
                                    filename = "TCGAvisualize_EAbarplot_Output.pdf",
                                    text.size = 1.0, mfrow = c(2, 2), xlim = NULL,
                                    color = c("orange", "cyan","green","yellow") ){

    if(!is.null(filename)) pdf(filename, width = 30, height = 15)

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
    par(mfrow = mfrow)

    if(!missing(GOBPTab)){
        if(!is.null(GOBPTab) & !is.na(GOBPTab)){
            # Plotting GOBPTab
            toPlot <- splitFun(tf, GOBPTab, nBar)
            xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[1],
                             main = "GO:Biological Process", xlab = "-log10(FDR)",xlim = xlim)
            labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
            text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
            lines(x = toPlot[, 3], y = xAxis, col = "red")
            points(x = toPlot[, 3], y = xAxis, col = "red")
            axis(side = 3, at = pretty(range(0:1)), col = "red")
        }
    }
    if(!missing(GOCCTab)){
        if(!is.null(GOCCTab) & !is.na(GOCCTab)){
            # Plotting GOCCTab
            toPlot <- splitFun(tf, GOCCTab, nBar)
            xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[2],
                             main = "GO:Cellular Component", xlab = "-log10(FDR)",xlim = xlim)
            labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
            text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
            lines(x = toPlot[, 3], y = xAxis, col = "red")
            points(x = toPlot[, 3], y = xAxis, col = "red")
            axis(side = 3, at = pretty(range(0:1)), col = "red")
        }
    }
    if(!missing(GOMFTab)){
        if(!is.null(GOMFTab) & !is.na(GOMFTab)){
            # Plotting GOMFTab
            toPlot <- splitFun(tf, GOMFTab, nBar)
            xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[3],
                             main = "GO:Molecular Function", xlab = "-log10(FDR)",xlim = xlim)
            labs <- matrix(unlist(strsplit(toPlot[, 1], "~")), nrow = 2)[2, ]
            text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
            lines(x = toPlot[, 3], y = xAxis, col = "red")
            points(x = toPlot[, 3], y = xAxis, col = "red")
            axis(side = 3, at = pretty(range(0:1)), col = "red")
        }
    }
    if(!missing(PathTab)){
        if(!is.null(PathTab) & !is.na(PathTab)){
            # Plotting PathTab
            toPlot <- splitFun(tf, PathTab, nBar)
            xAxis <- barplot(toPlot[, 2], horiz = TRUE, col = color[4],
                             main = "Pathways", xlab = "-log10(FDR)",xlim = xlim)
            labs <- toPlot[, 1]
            text(x = 1, y = xAxis, labs, pos = 4, cex = text.size)
            lines(x = toPlot[, 3], y = xAxis, col = "red")
            points(x = toPlot[, 3], y = xAxis, col = "red")
            #axis(side = 1, at = pretty(range(0:1)), col = "red", line = 2.5)
            axis(side = 3, at = pretty(range(0:1)), col = "red")
        }
    }
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

    if(!is.null(filename)) dev.off()
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
TCGAvisualize_BarPlot <- function(DFfilt,
                                  DFclin,
                                  DFsubt,
                                  data_Hc2,
                                  Subtype,
                                  cbPalette,
                                  filename,
                                  width,
                                  height,
                                  dpi){

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
#' @param data The object to with the heatmap data (expression, methylation)
#' @param col.metadata Metadata for the columns (samples). It should have on of the following columns:
#' barcode (28 characters)  column to match with the samples. It will also work with
#' "bcr_patient_barcode"(12 chars),"patient"(12 chars),"sample"(16 chars) columns but as one patient might
#' have more than one sample, this coul lead to errors in the annotation.
#' The code will throw a warning in case two samples are from the same patient.
#' @param row.metadata  Metadata for the rows  genes (expression) or probes (methylation)
#' @param col.colors A list of names colors
#' @param row.colors A list of named colors
#' @param type Select the colors of the heatmap values. Possible values are
#'  "expression" (default), "methylation"
#' @param show_column_names  Show column names names? Dafault: FALSE
#' @param show_row_names Show row names? Dafault: FALSE
#' @param cluster_rows Cluster rows ? Dafault: FALSE
#' @param cluster_columns Cluster columns ? Dafault: FALSE
#' @param filename Filename to save the heatmap. Default: heatmap.png
#' @param width figure width
#' @param height figure height
#' @param sortCol Name of the column to be used to sort the columns
#' @param title Title of the plot
#' @param rownames.size Rownames size
#' @param color.levels A vector with the colors (low level, middle level, high level)
#' @param extrems Extrems of colors (vector of 3 values)
#' @param values.label Text of the levels in the heatmap
#' @param heatmap.legend.color.bar Heatmap legends values type.
#' Options: "continuous", "disctrete
#' @param scale Use z-score to make the heatmap?
#' If we want to show differences between genes, it is good to make Z-score by samples
#' (force each sample to have zero mean and standard deviation=1).
#' If we want to show differences between samples, it is good to make Z-score by genes
#' (force each gene to have zero mean and standard deviation=1).
#' Possibilities: "row", "col". Default "none"
#' @examples
#'  row.mdat <- matrix(c("FALSE","FALSE",
#'                      "TRUE","TRUE",
#'                      "FALSE","FALSE",
#'                      "TRUE","FALSE",
#'                      "FALSE","TRUE"
#'                 ),
#'               nrow = 5, ncol = 2, byrow = TRUE,
#'               dimnames = list(
#'                   c("probe1", "probe2","probe3","probe4","probe5"),
#'                   c("duplicated", "Enhancer region")))
#' dat <- matrix(c(0.3,0.2,0.3,1,1,0.1,1,1,0, 0.8,1,0.7,0.7,0.3,1),
#'              nrow = 5, ncol = 3, byrow = TRUE,
#'                dimnames = list(
#'                c("probe1", "probe2","probe3","probe4","probe5"),
#'                c("TCGA-DU-6410",
#'                  "TCGA-DU-A5TS",
#'                  "TCGA-HT-7688")))
#'
#' mdat <- data.frame(patient=c("TCGA-DU-6410","TCGA-DU-A5TS","TCGA-HT-7688"),
#'                    Sex=c("Male","Female","Male"),
#'                    COCCluster=c("coc1","coc1","coc1"),
#'                    IDHtype=c("IDHwt","IDHMut-cod","IDHMut-noncod"))
#'
#'TCGAvisualize_Heatmap(dat,
#'                     col.metadata = mdat,
#'                     row.metadata = row.mdat,
#'                     row.colors = list(duplicated = c("FALSE" = "pink",
#'                                                      "TRUE"="green"),
#'                                      "Enhancer region" = c("FALSE" = "purple",
#'                                                             "TRUE"="grey")),
#'                     col.colors = list(Sex = c("Male" = "blue", "Female"="red"),
#'                                       COCCluster=c("coc1"="grey"),
#'                                       IDHtype=c("IDHwt"="cyan",
#'                                       "IDHMut-cod"="tomato"
#'                                       ,"IDHMut-noncod"="gold")),
#'                     type = "methylation",
#'                     show_row_names=TRUE)
#' if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
#' @export
#' @importFrom matlab jet.colors
#' @importFrom circlize colorRamp2
#' @importFrom tools file_ext
#' @import ComplexHeatmap
#' @return Heatmap plotted in the device
TCGAvisualize_Heatmap <- function(data,
                                  col.metadata,
                                  row.metadata,
                                  col.colors=NULL,
                                  row.colors=NULL,
                                  show_column_names = FALSE,
                                  show_row_names = FALSE,
                                  cluster_rows = FALSE,
                                  cluster_columns = FALSE,
                                  sortCol,
                                  extrems = NULL,
                                  rownames.size = 12,
                                  title = NULL,
                                  color.levels = NULL,
                                  values.label = NULL,
                                  filename = "heatmap.pdf",
                                  width = 10,
                                  height = 10,
                                  type = "expression",
                                  scale = "none",
                                  heatmap.legend.color.bar = "continuous"){

    # STEP 1 add columns labels (top of heatmap)

    ha <-  NULL
    if(!missing(col.metadata)) {
        if(!is.null(col.metadata)) {
            id <- NULL
            if("patient"  %in% colnames(col.metadata)) {
                id <- "patient"
                size <- 12
            }
            if("barcode"  %in% colnames(col.metadata)) {
                id <- "barcode"
                stopifnot(nchar(col.metadata[,id])[1] == 28)
                size <- 28
            }
            if("bcr_patient_barcode"  %in% colnames(col.metadata)) {
                id <- "bcr_patient_barcode"
                stopifnot(nchar(col.metadata[,id])[1] == 12)
                size <- 12
            }
            if("sample"  %in% colnames(col.metadata)) {
                id <- "sample"
                stopifnot(nchar(col.metadata[,id])[1] == 16)
                size <- 16
            }

            if(is.null(id)) {
                message("=============== INNPUT ERROR =================")
                message("I'm expecting one of these columns:")
                message(" => barcode")
                message("    Has the complete barcode (TCGA-AA-3833-01A-01D-0904-05)")
                message(" => bcr_patient_barcode")
                message("    Has the patient barcode (TCGA-AA-3833)")
                message(" => patient")
                message("    Has the patient barcode (TCGA-AA-3833)")
                message(" => sample")
                message("    Has the sample barcode (TCGA-AA-3833-01A)")
                message("-----------------------------------------------")
                message("Obs: The complete barcode is the recommended one, as the others might lead to errors")
                return(NULL)
            }
            stopifnot(nchar(as.character(col.metadata[,id])[1]) == size)
            message(paste0("Reorganizing: col.metadata order should be the same of the data object"))
            df <- col.metadata[match(substr(colnames(data),1,size), col.metadata[,id]),]
            df[,id] <- NULL

            duplicated.samples <- any(sapply(col.metadata[,id],
                                             function(x) {length(grep(x,col.metadata[,id])) > 1 }))
            if(duplicated.samples){
                warning("Some samples are from the same patient, this might lead to the wrong upper annotation")
            }

            if (!missing(sortCol)) {
                message(paste0("Sorting columns based on column: ",
                               sortCol))
                column_order <- order(df[,sortCol])
            }
            if(is.null(col.colors)) {
                ha <- HeatmapAnnotation(df = df)
            } else {
                ha <- HeatmapAnnotation(df = df,
                                        col = col.colors)
            }
        }
    }
    # STEP 2 Create heatmap
    if(is.null(color.levels)) {
        if (type == "expression") color.levels <-  c("green", "white", "red")
        if (type == "methylation") color.levels <-  c("blue", "white", "red")
    }


    # If we want to show differences between genes, it is good to make Z-score by samples
    # (force each sample to have zero mean and standard deviation=1).
    # If we want to show differences between samples, it is good to make Z-score by genes
    # (force each gene to have zero mean and standard deviation=1).
    if(scale == "row"){
        message("Calculating z-scores for the rows....")
        data <- t(scale(t(data)))
        all.na <- apply(data,1, function(x) all(is.na(x)))
        data <- data[!all.na,]
    } else if(scale == "col"){
        message("Calculiating z-scores for the columns....")
        data <- scale(data)
    }

    if(is.null(extrems)) {
        if(min(data) < 0) {
            extrems <- c(min(data), 0, max(data))
        } else {
            extrems <- c(0, max(data)/2, max(data))
        }
    }
    if (type == "expression") color <- circlize::colorRamp2(extrems, color.levels)
    if (type == "methylation") color <- circlize::colorRamp2(extrems, color.levels)

    # Creating plot title
    if(is.null(title)) {
        if(type == "methylation") title <- "DNA methylation heatmap"
        if(type == "expression") title <- "Expression heatmap"
    }
    # Change label type
    heatmap_legend_param <- list()
    if(heatmap.legend.color.bar == "continuous" && type == "methylation"){
        heatmap_legend_param <- c(list(color_bar = "continuous"),heatmap_legend_param)
        if(!scale %in% c("row","col")) heatmap_legend_param <- list(color_bar = "continuous", at = c(0,0.2,0.4,0.6,0.8, 1), legend_height = unit(3, "cm"), labels = c("0.0 (hypomethylated)",0.2,0.4,0.6,0.8,"1.0 (hypermethylated)"))
    }
    if(heatmap.legend.color.bar == "continuous" && type == "expression"){
        heatmap_legend_param <- c(list(color_bar = "continuous"),heatmap_legend_param)
    }
    # Change label reference
    if(is.null(values.label)){
        if(type == "methylation") values.label <- "DNA methylation level"
        if(type == "expression") values.label <- "Expression"
    }
    if(!missing(sortCol) & heatmap.legend.color.bar == "continuous"){
        heatmap  <- Heatmap(data, name = values.label,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            row_names_gp =  gpar(fontsize = rownames.size),
                            show_row_names = show_row_names,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_order = column_order,
                            column_title = title,
                            heatmap_legend_param = heatmap_legend_param)
    } else if(missing(sortCol) & heatmap.legend.color.bar == "continuous"){
        heatmap  <- Heatmap(data, name = values.label,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            show_row_names = show_row_names,
                            row_names_gp =  gpar(fontsize = rownames.size),
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_title = title,
                            heatmap_legend_param = heatmap_legend_param)
    }  else if(!missing(sortCol)){
        heatmap  <- Heatmap(data, name = values.label,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            row_names_gp =  gpar(fontsize = rownames.size),
                            show_row_names = show_row_names,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_order = column_order,
                            column_title = title)
    } else {
        heatmap  <- Heatmap(data, name = values.label,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            row_names_gp =  gpar(fontsize = rownames.size),
                            show_row_names = show_row_names,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_title = title,
                            heatmap_legend_param = heatmap_legend_param)
    }

    # STEP 3 row labels (right side)
    if (!missing(row.metadata)) {
        if (!is.null(row.metadata)) {
            for (i in 1:ncol(row.metadata)) {
                if (!missing(row.colors) && !is.null(row.colors[[colnames(row.metadata)[i]]])) {
                    color <- row.colors[[colnames(row.metadata)[i]]]
                    x = Heatmap(row.metadata[,i] ,
                                name = colnames(row.metadata)[i],
                                width = unit(0.5, "cm"),
                                show_row_names = FALSE, col = color )
                } else {
                    x = Heatmap(row.metadata[,i] ,
                                name = colnames(row.metadata)[i],
                                width = unit(0.5, "cm"),
                                show_row_names = FALSE)
                }
                heatmap <- add_heatmap(heatmap,x)
            }
        }
    }
    if(!is.null(filename)){
        if(file_ext(filename) == "png") png(filename, width = width, height = height )
        if(file_ext(filename) == "pdf") pdf(filename, width = width, height = height )
        draw(heatmap)
        dev.off()
    } else {
        draw(heatmap)
    }
}


# unlist labels
# Help function that unlists a list into a vector
unlistlabels <- function(lab) {
    dummy <- unlist(lab)
    labels <- c()
    labels <- c(labels, as.character(dummy))
    return(labels)
}


#' Creating a oncoprint
#' @param mut A dataframe from the mutation annotation file (see TCGAquery_maf from TCGAbiolinks)
#' @param genes Gene list
#' @param filename name of the pdf
#' @param color named vector for the plot
#' @param height pdf height
#' @param width pdf width
#' @param rm.empty.columns If there is no alteration in that sample, whether remove it on the oncoprint
#' @param show.row.barplot  Show barplot annotation on rows?
#' @param show.column.names Show column names? Default: FALSE
#' @param rows.font.size Size of the fonts
#' @param dist.col distance between columns in the plot
#' @param dist.row distance between rows in the plot
#' @param label.font.size Size of the fonts
#' @param row.order Order the genes (rows). Genes with more mutations will be in the first rows
#' @param annotation Matrix or data frame with the annotation.
#' Should have a column bcr_patient_barcode with the same ID of the mutation object
#' @param annotation.position Position of the annotation "bottom" or "top"
#' @param label.title Title of the label
#' @param annotation.legend.side Position of the annotation legend
#' @param heatmap.legend.side Position of the heatmap legend
#' @param information Which column to use as informastion from MAF.
#' Options: 1) "Variant_Classification" (The information will be "Frame_Shift_Del", "Frame_Shift_Ins",
#'         "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",  "Nonsense_Mutation",
#'              "Nonstop_Mutation",  "RNA",  "Silent" ,  "Splice_Site",  "Targeted_Region",  "Translation_Start_Site")
#' 2) "Variant_Type" (The information will be INS,DEL,SNP)
#' @importFrom ComplexHeatmap oncoPrint draw HeatmapAnnotation
#' @importFrom grid gpar grid.rect
#' @importFrom data.table dcast setDT setDF :=
#' @examples
#' mut <- GDCquery_Maf(tumor = "ACC")
#' TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10], rm.empty.columns = TRUE)
#' TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:10],
#'                  filename = "onco.pdf",
#'                  color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"))
#' clin <- GDCquery_clinic("TCGA-ACC","clinical")
#' clin <- clin[,c("bcr_patient_barcode","disease","gender","tumor_stage","race","vital_status")]
#' TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:20],
#'                 filename = "onco.pdf",
#'                 annotation = clin,
#'                 color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
#'                 rows.font.size=10,
#'                 heatmap.legend.side = "right",
#'                 dist.col = 0,
#'                 label.font.size = 10)
#'
#' @export
#' @return A oncoprint plot
TCGAvisualize_oncoprint <- function (mut,
                                     genes,
                                     filename,
                                     color,
                                     annotation.position = "bottom",
                                     annotation,
                                     height,
                                     width = 10,
                                     rm.empty.columns = FALSE,
                                     show.column.names = FALSE,
                                     show.row.barplot = TRUE,
                                     label.title = "Mutation",
                                     label.font.size = 16,
                                     rows.font.size = 16,
                                     dist.col = 0.5,
                                     dist.row = 0.5,
                                     information = "Variant_Type",
                                     row.order = FALSE,
                                     heatmap.legend.side = "bottom",
                                     annotation.legend.side = "bottom"){


    if(missing(mut))   stop("Missing mut argument")
    mut <- setDT(mut)
    mut$value <- 1

    mut$Hugo_Symbol <- as.character(mut$Hugo_Symbol)
    if(!missing(genes) & !is.null(genes)) mut <- subset(mut, mut$Hugo_Symbol %in% genes)

    if(!rm.empty.columns){
        formula <- paste0("Tumor_Sample_Barcode + Hugo_Symbol ~ ", information)
        mat <- dcast(mut, as.formula(formula),value.var = "value",fill = 0,drop = FALSE)
    } else {
        formula <- paste0("Tumor_Sample_Barcode + Hugo_Symbol ~ ", information)
        mat <- dcast(mut, as.formula(formula),value.var = "value",fill = 0,drop = TRUE)
    }

    # mutation in the file
    columns <- colnames(mat)[-c(1:2)]

    # value will be a collum with all the mutations
    mat$value <- ""

    for ( i in columns){
        mat[,i] <-  replace(mat[,i,with = FALSE],mat[,i,with = FALSE]>0,paste0(i,";"))
        mat[,i] <-  replace(mat[,i,with = FALSE],mat[,i,with = FALSE]==0,"")
        mat[,value:=paste0(value,get(i))]
    }

    # After the gene selection, some of the mutation might not exist
    # we will remove them to make the oncoprint work
    mutation.type <- c()
    for (i in columns){
        if(length(grep(i,mat$value)) > 0) mutation.type <- c(mutation.type,i)
    }

    # now we have a matrix with pairs samples/genes mutations
    # we want a matrix with samples vs genes mutations with the content being the value
    mat <- setDF(dcast(mat, Tumor_Sample_Barcode~Hugo_Symbol, value.var="value",fill=""))
    rownames(mat) <- mat[,1]
    mat <- mat[,-1]


    alter_fun = function(x, y, w, h, v) {
        n = sum(v)
        h = h*0.9
        # use `names(which(v))` to correctly map between `v` and `col`
        if(n) {
            grid.rect(x, y - h*0.5 + 1:n/n*h,  w-unit(dist.col, "mm"), 1/n*h,
                      gp = gpar(fill = color[names(which(v))], col = NA), just = "top")
        } else {
            grid.rect(x, y, w-unit(dist.col, "mm"), h-unit(dist.row, "mm"), gp = gpar(fill = color["background"], col = NA))
        }
    }

    # get only the colors to the mutations
    # otherwise it gives errors

    if(missing(color)){
        color <- c(rainbow(length(mutation.type)), "#CCCCCC")
        names(color) <- c(mutation.type,"background")
    } else{
        if("background" %in% names(color)) {
            color <- color[c(mutation.type,"background")]
        } else {
            color <- c(color[mutation.type],"background"= "#CCCCCC")
        }
    }
    # header are samples, rows genes
    mat <- t(mat)

    if(!missing(height)) height <- length(genes)/2
    if(!missing(filename)) pdf(filename,width = width,height = height)

    if(missing(annotation)) annotation <- NULL
    if(!is.null(annotation)){
        idx <- match(substr(colnames(mat),1,12),annotation$bcr_patient_barcode)

        annotation <- annotation[idx,]

        annotation$bcr_patient_barcode <- NULL

        n.col <- sum(sapply(colnames(annotation), function(x) {
            length(unique(annotation[,x]))
        }))

        # add automatic colors: not working

        get.color <- function(df,col){
            idx <- which(colnames(df) == col)
            start <- 1
            if(idx != 1) start <- length(unique(unlist(c(df[,1:(idx-1)])))) + 1
            end <- start + length(unique(df[,col])) -1
            diff.colors <- c("dimgray","thistle","deeppink3","magenta4","lightsteelblue1","black",
                             "chartreuse","lightgreen","maroon4","darkslategray",
                             "lightyellow3","darkslateblue","firebrick1","aquamarine",
                             "dodgerblue4","bisque4","moccasin","indianred1",
                             "yellow","gray93","cyan","darkseagreen4",
                             "lightgoldenrodyellow","lightpink","sienna1",
                             "darkred","palevioletred","tomato4","blue",
                             "mediumorchid4","royalblue1","magenta2","darkgoldenrod1")
            return(diff.colors[start:end])
        }
        col.annot <- lapply(colnames(annotation), function(x) {
            #idx <- which(colnames(annotation) == x) - 1
            #print(idx/n.col)
            ret <- get.color(annotation,x)
            #ret <- rainbow(length(unique(annotation[,x])),start = idx/n.col,alpha=0.5)
            names(ret) <- as.character(unique(annotation[,x]))
            return(ret)
        })
        names(col.annot) <-  colnames(annotation)

        annotHeatmap <- HeatmapAnnotation(df=annotation,
                                          col=col.annot,
                                          annotation_legend_param=list(title_gp=gpar(fontsize=label.font.size,
                                                                                     fontface="bold"),
                                                                       labels_gp=gpar(fontsize=label.font.size),#sizelabels
                                                                       grid_height=unit(8,"mm"))
        )
    }
    if(heatmap.legend.side == "bottom") {
        nrow <- 1
        title_position <- "leftcenter"
    } else {
        nrow <- 10
        title_position <- "topcenter"
    }
    if(is.null(annotation) & !row.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom" & !row.order){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )

    } else if(!is.null(annotation) & annotation.position == "top" & !row.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       row_order = NULL,
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"),  # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    }  else if(is.null(annotation) & row.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    } else if(!is.null(annotation) & annotation.position == "bottom" & row.order){

        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_row_barplot = show.row.barplot,
                       show_column_names = show.column.names,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       bottom_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"), # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )

    } else if(!is.null(annotation) & annotation.position == "top" & row.order){
        p <- oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                       remove_empty_columns = FALSE,
                       show_column_names = show.column.names,
                       show_row_barplot = show.row.barplot,
                       column_order = NULL, # Do not sort the columns
                       alter_fun = alter_fun, col = color,
                       row_names_gp = gpar(fontsize = rows.font.size),  # set size for row names
                       pct_gp = gpar(fontsize = rows.font.size), # set size for percentage labels
                       axis_gp = gpar(fontsize = rows.font.size),# size of axis
                       #column_title = "OncoPrint for TCGA LGG, genes in Glioma signaling",
                       #column_title_gp = gpar(fontsize = 11),
                       row_barplot_width = unit(2, "cm"), #size barplot
                       top_annotation = annotHeatmap,
                       heatmap_legend_param = list(title = label.title, at = names(color),
                                                   labels = names(color),
                                                   title_gp = gpar(fontsize = label.font.size, fontface = "bold"),
                                                   labels_gp = gpar(fontsize = label.font.size), # size labels
                                                   grid_height = unit(8, "mm"),  # vertical distance labels
                                                   nrow = nrow, title_position = title_position
                       )
        )
    }

    draw(p, heatmap_legend_side = heatmap.legend.side, annotation_legend_side = annotation.legend.side)
    if(!missing(filename)) dev.off()

}
