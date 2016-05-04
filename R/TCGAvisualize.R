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
        if(!is.null(GOBPTab)){
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
        if(!is.null(GOCCTab)){
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
        if(!is.null(GOMFTab)){
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
        if(!is.null(PathTab)){
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
#' @param sortCol Name of the column to be used to sort the columns
#' @param title Title of the plot
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
                                  title,
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
            ha <- HeatmapAnnotation(df = df,
                                    col = col.colors)
        }
    }
    # STEP 2 Create heatmap

    # If we want to show differences between genes, it is good to make Z-score by samples
    # (force each sample to have zero mean and standard deviation=1).
    # If we want to show differences between samples, it is good to make Z-score by genes
    # (force each gene to have zero mean and standard deviation=1).
    if(scale == "row"){
        message("Calculiating z-scores for the rows....")
        data <- t(scale(t(data)))
        if (type == "expression") color <- colorRamp2(seq(-4,4,0.1), gplots::greenred(length(seq(-4,4,0.1))))
        if (type == "methylation") color <- colorRamp2(seq(-4,4,0.1), matlab::jet.colors(length(seq(-4,4,0.1))))
    } else if(scale == "col"){
        message("Calculiating z-scores for the columns....")
        data <- scale(data)
        if (type == "expression") color <- colorRamp2(seq(-4,4,0.1), gplots::greenred(length(seq(-4,4,0.1))))
        if (type == "methylation") color <- colorRamp2(seq(-4,4,0.1), matlab::jet.colors(length(seq(-4,4,0.1))))
    } else {
        if (type == "expression") color <- gplots::greenred(200)
        if (type == "methylation") color <- matlab::jet.colors(200)
    }

    # Creating plot title
    if(missing(title)) {
        if(type == "methylation") title <- "Methylation heatmap"
        if(type == "expression") title <- "Expression heatmap"
    }

    # Change label type
    if(heatmap.legend.color.bar == "continuous" && type == "methylation"){
        heatmap_legend_param <- list(color_bar = "continuous", at = c(0,0.2,0.4,0.6,0.8, 1), legend_height = unit(3, "cm"), labels = c("0.0 (hypomethylated)",0.2,0.4,0.6,0.8,"1.0 (hypermethylated)"))
    }
    if(heatmap.legend.color.bar == "continuous" && type == "expression"){
        heatmap_legend_param <- list(color_bar = "continuous")
    }

    # Change label reference
    if(type == "methylation") type <- "Methylation level"
    if(type == "expression") type <- "Expression"

    if(!missing(sortCol) & heatmap.legend.color.bar == "continuous"){
        heatmap  <- Heatmap(data, name = type,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            show_row_names = show_row_names,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_order = column_order,
                            column_title = title,
                            heatmap_legend_param = heatmap_legend_param)
    } else if(missing(sortCol) & heatmap.legend.color.bar == "continuous"){
        heatmap  <- Heatmap(data, name = type,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            show_row_names = show_row_names,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_title = title,
                            heatmap_legend_param = heatmap_legend_param)
    }  else if(!missing(sortCol)){
        heatmap  <- Heatmap(data, name = type,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
                            show_row_names = show_row_names,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_column_names = show_column_names,
                            column_order = column_order,
                            column_title = title)
    } else {
        heatmap  <- Heatmap(data, name = type,
                            top_annotation = ha,
                            bottom_annotation_height = unit(3, "cm"),
                            col = color,
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
    return(heatmap)
}

#' @title Profile plot
#' @description Displaty the association between cancer subtypes and any kind of clustering.
#' @param data A data frame with the cluters and subytpe of cancers
#' @param subtypeCol Name of the column with the subtype information
#' @param groupCol Names of tre columns with the cluster information
#' @param filename Name of the file to save the plot, can be pdf, png, svg etc..
#' @param na.rm.groups Remove NA groups? Default = FALSE
#' @param na.rm.subtypes Remove NA subtypes? Default = FALSE
#' @param colors Vector of colors to be used in the bars
#' @param plot.margin Plot margin for cluster distribution. This can control the size
#' of the bar if the output is not aligned
#' @param axis.title.size axis.title.size
#' @param axis.textsize axis.textsize
#' @param legend.size Size of the legend
#' @param legend.title.size Size of the legend title
#' @param geom.label.size Size of percentage in the left barplot
#' @param geom.label.color Color of percentage in the left barplot
#' @importFrom sjPlot sjp.stackfrq sjp.setTheme
#' @importFrom cowplot ggdraw switch_axis_position plot_grid
#' @importFrom reshape2 dcast
#' @importFrom grDevices gray.colors
#' @import gtable
#' @export
#' @examples
#' while (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
#' cluster <- c(rep("cluster1",30),
#'              rep("cluster2",30),
#'              rep("cluster3",30))
#' subtype <- rep(c(rep("subtype1",10),
#'            rep("subtype2",10),
#'            rep("subtype3",10)),3)
#' df <- data.frame(cluster,subtype)
#' TCGAvisualize_profilePlot(data = df, groupCol = "cluster", subtypeCol = "subtype",
#'                           plot.margin=c(-4.2,-2.5,-0.0,2))
#' while (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
#' cluster <- c(rep("cluster1",10),
#'              rep("cluster2",20),
#'              rep("cluster3",30),
#'              rep("cluster4",40))
#' subtype <- rep(c(rep("subtype1",5),
#'            rep("subtype2",10),
#'            rep("subtype3",10)),4)
#' df <- data.frame(cluster,subtype)
#' plot <- TCGAvisualize_profilePlot(data = df, groupCol = "cluster", subtypeCol = "subtype",
#'                           plot.margin=c(-4.2,-2.5,-0.5,2))
#' @return A plot
TCGAvisualize_profilePlot <- function(data = NULL,
                                      groupCol = NULL,
                                      subtypeCol = NULL,
                                      colors = NULL,
                                      filename = NULL,
                                      na.rm.groups = FALSE,
                                      na.rm.subtypes= FALSE,
                                      plot.margin=c(-2.5,-2.5,-0.5,2),
                                      axis.title.size=1.5,
                                      axis.textsize=1.3,
                                      legend.size=1.5,
                                      legend.title.size=1.5,
                                      geom.label.size = 6.0,
                                      geom.label.color = "black") {

    sjp.setTheme(theme = "scatterw",
                 axis.title.size = axis.title.size,
                 axis.textsize = axis.textsize,
                 legend.size = legend.size,
                 legend.title.size = legend.title.size,
                 geom.label.size = geom.label.size,
                 geom.label.color = geom.label.color)

    if (is.null(groupCol)) stop("Please provide the groupCol argument")
    if (is.null(subtypeCol)) stop("Please provide the subtypeCol argument")
    if (is.null(data)) stop("Please provide the data argument")
    if (is.null(filename)) filename <- paste0(groupCol,subtypeCol,".pdf")

    if(na.rm.groups){
        data <- data[!is.na(data[,groupCol]),]
        data <- data[which(data[,groupCol] != "NA"),]
    }
    if(na.rm.subtypes){
        data <- data[!is.na(data[,subtypeCol]),]
        data <- data[which(data[,subtypeCol] != "NA"),]
    }

    # use https://github.com/cttobin/ggthemr
    # when it is in cran
    if (is.null(colors)) colors <- c("#34495e",
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
    groups <- df[,groupCol]

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

    ngroups <- length(unique(groups))
    nsbutype <- length(unique(all))

    if(NA %in% unique(all)) nsbutype <- nsbutype - 1
    max <- length(na.omit(data[,1]))
    message(paste0("Number of subtypes: ", nsbutype))
    message(paste0("Number of ngroups: ", ngroups))

    for(i in 1:ncol(data)){
        if(i == 1) {
            width = rep(c(0.5),nsbutype)
        }  else {
            width <- c(width, rep((ngroups/max) * length(na.omit(data[,i])),nsbutype))
        }
    }
    var.labels <- as.character(var.labels)
    var.labels[which(is.na(var.labels)==TRUE)] <- "NA"
    # Create the horizontal barplot
    p <- .mysjp.stackfrq(data,
                         legendTitle = subtypeCol,
                         #axisTitle.x = groupCol,
                         axisTitle.x = "",
                         #sort.frq = "last.desc",
                         expand.grid = FALSE,
                         geom.size = width,
                         legendLabels = as.character(var.labels),
                         jitterValueLabels = TRUE,
                         #showSeparatorLine = TRUE,
                         showValueLabels = FALSE,
                         geom.colors = colors[1:length(var.labels)])$plot
    p <-  ggdraw(switch_axis_position(p , axis = 'y'))

    j <- 1
    groups <- as.data.frame(groups)
    groups$x <- 1
    for (i in sort(unique(groups[,1]))){
        idx <- which(groups[,1] == i)
        groups[idx,"x"] <- as.numeric(j)
        j <- j + 1
    }

    # Create the vertical barplot with the percentage of element in each group
    p2 <- sjp.stackfrq(groups[,2],
                       #legendTitle = subtypeCol,
                       axisTitle.y = "Cluster distribution",
                       #sort.frq = "first.asc",
                       #expand.grid = TRUE,
                       legendLabels = as.character(unique(groups[,2])),
                       coord.flip = FALSE,
                       hideLegend = TRUE,
                       showSeparatorLine = FALSE,
                       #showValueLabels = FALSE,
                       includeN = FALSE,
                       #geom.size = c(0.1,0.2,0.4,0.3),
                       geom.colors = gray.colors(length(unique(groups[,2])),
                                                 start = 0.6,
                                                 end = 0.9,
                                                 gamma = 2.2,
                                                 alpha = 0.1))$plot

    p2 <- p2 +
        theme(#legend.position="none",
            #plot.margin=unit(c(-3.0,4.5,-0.75,5.5), "cm"),
            #plot.margin=unit(c(-3.3,-2.5,-1.0,2), "cm"),
            plot.margin=unit(plot.margin, "cm"),
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            #panel.background=element_blank(),
            #panel.border=element_blank(),
            #panel.grid.major=element_blank(),
            #panel.grid.minor=element_blank(),
            plot.background=element_blank())

    # put the plots together
    p <-  plot_grid(p2,
                    p,
                    ncol = 2,
                    scale = c(0.5,1),
                    rel_heights = c(2,8),
                    rel_widths = c(0.6,2))
    plot(p)
    if(!is.null(filename))
        ggsave(p, filename = filename, width = 20, height = 10, dpi = 600)
}




#' @title Visualize mutation
#' @description See % of genes mutated
#' @param data A data frame with the cluters and subytpe of cancers
#' @param samples Samples to consider mutation
#' @param colors Vector of colors for the barplot
#' @param geneList List of genes to plot
#' @param filename Name of the file to save the plot, can be pdf, png, svg etc..
#' @param threshold plot only if number of mutated genes are higher than threshold
#' @param by type of plot. Options: "genes" and "cluster".
#'  By "genes" the axis will receive the genes
#' and how many samples are mutated. By "cluster" will show the mutations of
#' only one gene in the cluster.
#' @param groupCol Must be provided if by is set to "cluster"
#' @param na.rm Remove NA groups
#' @importFrom sjPlot sjp.stackfrq sjp.setTheme
#' @export
#' @return A plot
TCGAvisualize_mutation <- function (data = NULL,
                                    colors = NULL,
                                    samples = NULL,
                                    geneList = NULL,
                                    filename = NULL,
                                    threshold = 0,
                                    by="genes",
                                    groupCol=NULL,
                                    na.rm = TRUE) {


    if(na.rm){
        data <- data[!is.na(data[,groupCol]),]
        data <- data[which(data[,groupCol] != "NA"),]
    }

    if (is.null(data)) stop("Please provide the data argument")
    if (is.null(filename)) filename <- "mutation_summary.pdf"

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

    if(by == "cluster" & !is.null(groupCol)){
        if(is.null(groupCol)){
            stop("Please provide groupCol argument")
        }

        if(length(geneList) != 1 ){
            stop("Please geneList in this case must be only one gene")
        }

        clusters <- unique(data[,groupCol])
        df <-  data.frame(matrix(NA, ncol = length(clusters), nrow = nrow(data)))
        colnames(df) <- clusters
        all <- c()
        for(j in clusters){
            idx <-  which(data[,groupCol] == j)
            aux <- data[idx,]
            summary <- table(unlist(aux$genes))
            count <- summary[geneList]
            df[,j] <- c(rep(1,as.numeric(count)),
                        rep(2,nrow(aux)-count),
                        rep(NA,nrow(data)- nrow(aux)))

            all <- c(all,c(rep(1,as.numeric(count))))
        }

        df[,groupCol] <- c(all,rep(2,nrow(data)-length(all)))
        df <- df[,rev(colnames(df))]

        # add theme when this library is in CRAN
        # https://github.com/cttobin/ggthemr
        #sjp.setTheme(theme = "539")
        p <- sjp.stackfrq(df,
                          title = paste0(geneList),
                          legendTitle = "Status",
                          #axisTitle.x = groupCol,
                          #sort.frq = "last.desc",
                          expand.grid = FALSE,
                          legendLabels = c("Mutated","Not mutated"),
                          showSeparatorLine = TRUE,
                          showValueLabels = FALSE,
                          geom.colors = colors[1:2],
                          #separatorLineColor = "#6699cc"
                          printPlot = TRUE)


    } else {
        df <-  data.frame(matrix(NA, ncol = length(genes), nrow = nrow(data)))
        colnames(df) <- genes
        summary <- table(unlist(data$genes))

        for( i in genes){

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
        p <- sjp.stackfrq(df,
                          #legendTitle = subtypeCol,
                          #axisTitle.x = groupCol,
                          #sort.frq = "last.desc",
                          expand.grid = FALSE,
                          legendLabels = c("Mutated","Not mutated"),
                          showSeparatorLine = TRUE,
                          showValueLabels = FALSE,
                          geom.colors = colors[1:2],
                          #separatorLineColor = "#6699cc"
                          printPlot = TRUE)

    }

    ggsave(p$plot, filename = filename, width = 10, height = 10, dpi = 600)
}

#' @import magrittr
#' @keywords internal
#' @import sjmisc
.mysjp.stackfrq <- function(items,
                            legendLabels = NULL,
                            sort.frq = NULL,
                            weightBy = NULL,
                            weightByTitleString = NULL,
                            hideLegend = FALSE,
                            title = NULL,
                            legendTitle = NULL,
                            includeN = TRUE,
                            axisLabels.y = NULL,
                            breakTitleAt = 50,
                            breakLabelsAt = 30,
                            breakLegendTitleAt = 30,
                            breakLegendLabelsAt = 28,
                            gridBreaksAt = 0.2,
                            expand.grid = FALSE,
                            geom.size = 0.5,
                            geom.colors = "Blues",
                            axisTitle.x = NULL,
                            axisTitle.y = NULL,
                            showValueLabels = TRUE,
                            labelDigits = 1,
                            showPercentageAxis = TRUE,
                            jitterValueLabels = FALSE,
                            showItemLabels = TRUE,
                            showSeparatorLine = FALSE,
                            separatorLineColor = "grey80",
                            separatorLineSize = 0.3,
                            coord.flip = TRUE,
                            printPlot = TRUE) {
    # --------------------------------------------------------
    # check param. if we have a single vector instead of
    # a data frame with several items, convert vector to data frame
    # --------------------------------------------------------
    if (!is.data.frame(items) && !is.matrix(items)) items <- as.data.frame(items)
    # --------------------------------------------------------
    # check sorting
    # --------------------------------------------------------
    if (!is.null(sort.frq)) {
        if (sort.frq == "first.asc") {
            sort.frq  <- "first"
            reverseOrder <- FALSE
        } else if (sort.frq == "first.desc") {
            sort.frq  <- "first"
            reverseOrder <- TRUE
        } else if (sort.frq == "last.asc") {
            sort.frq  <- "last"
            reverseOrder <- TRUE
        } else if (sort.frq == "last.desc") {
            sort.frq  <- "last"
            reverseOrder <- FALSE
        } else {
            sort.frq  <- NULL
            reverseOrder <- FALSE
        }
    } else {
        reverseOrder <- FALSE
    }
    # --------------------------------------------------------
    # try to automatically set labels is not passed as parameter
    # --------------------------------------------------------
    if (is.null(legendLabels)) legendLabels <- sjmisc::get_labels(items[[1]],
                                                                  attr.only = F,
                                                                  include.values = NULL,
                                                                  include.non.labelled = T)
    if (is.null(axisLabels.y)) {
        axisLabels.y <- c()
        # if yes, iterate each variable
        for (i in 1:ncol(items)) {
            # retrieve variable name attribute
            vn <- sjmisc::get_label(items[[i]], def.value = colnames(items)[i])
            # if variable has attribute, add to variableLabel list
            if (!is.null(vn)) {
                axisLabels.y <- c(axisLabels.y, vn)
            } else {
                # else break out of loop
                axisLabels.y <- NULL
                break
            }
        }
    }
    # --------------------------------------------------------
    # If axisLabels.y were not defined, simply use column names
    # --------------------------------------------------------
    if (is.null(axisLabels.y)) axisLabels.y <- colnames(items)
    # --------------------------------------------------------
    # unlist/ unname axis labels
    # --------------------------------------------------------
    if (!is.null(axisLabels.y)) {
        # unlist labels, if necessary, so we have a simple
        # character vector
        if (is.list(axisLabels.y)) axisLabels.y <- unlistlabels(axisLabels.y)
        # unname labels, if necessary, so we have a simple
        # character vector
        if (!is.null(names(axisLabels.y))) axisLabels.y <- as.vector(axisLabels.y)
    }
    # --------------------------------------------------------
    # unlist/ unname axis labels
    # --------------------------------------------------------
    if (!is.null(legendLabels)) {
        # unlist labels, if necessary, so we have a simple
        # character vector
        if (is.list(legendLabels)) legendLabels <- unlistlabels(legendLabels)
        # unname labels, if necessary, so we have a simple
        # character vector
        if (!is.null(names(legendLabels))) legendLabels <- as.vector(legendLabels)
    }
    if (is.null(legendLabels)) {
        # if we have no legend labels, we iterate all data frame's
        # columns to find all unique items of the data frame.
        # In case one item has missing categories, this may be
        # "compensated" by looking at all items, so we have the
        # actual values of all items.
        legendLabels <- as.character(sort(unique(unlist(
            apply(items, 2, function(x) unique(stats::na.omit(x)))))))
    }
    # --------------------------------------------------------
    # Check whether N of each item should be included into
    # axis labels
    # --------------------------------------------------------
    if (includeN && !is.null(axisLabels.y)) {
        for (i in 1:length(axisLabels.y)) {
            axisLabels.y[i] <- paste(axisLabels.y[i],
                                     sprintf(" (n=%i)", length(stats::na.omit(items[[i]]))),
                                     sep = "")
        }
    }
    # -----------------------------------------------
    # if we have legend labels, we know the exact
    # amount of groups
    # -----------------------------------------------
    countlen <- length(legendLabels)
    # -----------------------------------------------
    # create cross table for stats, summary etc.
    # and weight variable. do this for each item that was
    # passed as parameter
    #---------------------------------------------------
    mydat <- c()
    # ----------------------------
    # determine minimum value. if 0, add one, because
    # vector indexing starts with 1
    # ----------------------------
    if (any(apply(items, c(1, 2), is.factor)) || any(apply(items, c(1, 2), is.character))) {
        diff <- ifelse(min(apply(items, c(1, 2), as.numeric), na.rm = TRUE) == 0, 1, 0)
    } else {
        diff <- ifelse(min(items, na.rm = TRUE) == 0, 1, 0)
    }
    # iterate item-list
    for (i in 1:ncol(items)) {
        # get each single items
        variable <- items[[i]]
        # -----------------------------------------------
        # create proportional table so we have the percentage
        # values that should be used as y-value for the bar charts
        # We now have a data frame with categories, group-association
        # and percentage values (i.e. each cell as separate row in the
        # data frame)
        # -----------------------------------------------
        # check whether counts should be weighted or not
        if (is.null(weightBy)) {
            df <- as.data.frame(prop.table(table(variable)))
        } else {
            df <- as.data.frame(prop.table(round(stats::xtabs(weightBy ~ variable), 0)))
        }

        grp <- NULL
        ypos <- NULL

        # give columns names
        names(df) <- c("var", "prc")
        # need to be numeric, so percentage values (see below) are
        # correctly assigned, i.e. missing categories are considered
        df$var <- sjmisc::to_value(df$var, keep.labels = FALSE) + diff # if categories start with zero, fix this here
        # Create a vector of zeros
        prc <- rep(0, countlen)
        # Replace the values in prc for those indices which equal df$var
        prc[df$var] <- df$prc
        # create new data frame. We now have a data frame with all
        # variable categories abd their related percentages, including
        # zero counts, but no(!) missings!
        mydf <- data.frame(grp = i,
                           cat = 1:countlen,
                           prc)
        # now, append data frames
        mydat <- data.frame(rbind(mydat, mydf))
    }
    # ----------------------------
    # make sure group and count variable
    # are factor values
    # ----------------------------
    mydat$grp <- as.factor(mydat$grp)
    mydat$cat <- as.factor(mydat$cat)
    # add half of Percentage values as new y-position for stacked bars
    mydat <- mydat %>%
        dplyr::group_by(grp) %>%
        dplyr::mutate(ypos = cumsum(prc) - 0.5 * prc) %>%
        dplyr::arrange(grp)
    # --------------------------------------------------------
    # Caculate vertical adjustment to avoid overlapping labels
    # --------------------------------------------------------
    jvert <- rep(c(1.1, -0.1), length.out = length(unique(mydat$cat)))
    jvert <- rep(jvert, length.out = nrow(mydat))
    # --------------------------------------------------------
    # Prepare and trim legend labels to appropriate size
    # --------------------------------------------------------
    # wrap legend text lines
    legendLabels <- sjmisc::word_wrap(legendLabels, breakLegendLabelsAt)
    # check whether we have a title for the legend
    # if yes, wrap legend title line
    if (!is.null(legendTitle)) legendTitle <- sjmisc::word_wrap(legendTitle, breakLegendTitleAt)
    # check length of diagram title and split longer string at into new lines
    # every 50 chars
    if (!is.null(title)) {
        # if we have weighted values, say that in diagram's title
        if (!is.null(weightByTitleString)) title <- paste0(title, weightByTitleString)
        title <- sjmisc::word_wrap(title, breakTitleAt)
    }
    # check length of x-axis-labels and split longer strings at into new lines
    # every 10 chars, so labels don't overlap
    if (!is.null(axisLabels.y)) axisLabels.y <- sjmisc::word_wrap(axisLabels.y, breakLabelsAt)
    # ----------------------------
    # Check if ordering was requested
    # ----------------------------
    if (!is.null(sort.frq)) {
        # order by first cat
        if (sort.frq == "first") {
            facord <- order(mydat$prc[which(mydat$cat == 1)])
        } else {
            # order by last cat
            facord <- order(mydat$prc[which(mydat$cat == countlen)])
        }
        # create dummy vectors from 1 to itemlength
        dummy1 <- dummy2 <- c(1:length(facord))
        # facords holds the ordered item indices! we now need to
        # change the original item-index with its ordered position index.
        # example:
        # we have 4 items, and they may be ordered like this:
        # 1 3 4 2
        # so the first item is the one with the lowest count , item 3 is on second postion,
        # item 4 is on third position and item 2 is the last item (with highest count)
        # we now need their order as subsequent vector: 1 4 2 3
        # (i.e. item 1 is on first pos, item 2 is on fourth pos, item 3 is on
        # second pos and item 4 is on third pos in order)
        if (reverseOrder) {
            dummy2[rev(facord)] <- dummy1
        } else {
            dummy2[facord] <- dummy1
        }
        # now we have the order of either lowest to highest counts of first
        # or last category of "items". We now need to repeat these values as
        # often as we have answer categories
        orderedrow <- unlist(tapply(dummy2, 1:length(dummy2), function(x) rep(x, countlen)))
        # replace old grp-order by new order
        mydat$grp <- as.factor(orderedrow)
        # reorder axis labels as well
        axisLabels.y <- axisLabels.y[order(dummy2)]
    }
    # --------------------------------------------------------
    # check if category-oder on x-axis should be reversed
    # change category label order then
    # --------------------------------------------------------
    if (reverseOrder && is.null(sort.frq)) axisLabels.y <- rev(axisLabels.y)
    # --------------------------------------------------------
    # define vertical position for labels
    # --------------------------------------------------------
    if (coord.flip) {
        # if we flip coordinates, we have to use other parameters
        # than for the default layout
        vert <- 0.35
    } else {
        vert <- waiver()
    }
    # --------------------------------------------------------
    # set diagram margins
    # --------------------------------------------------------
    if (expand.grid) {
        expgrid <- waiver()
    } else {
        expgrid <- c(0, 0)
    }
    # --------------------------------------------------------
    # Set value labels and label digits
    # --------------------------------------------------------
    mydat$labelDigits <- labelDigits
    if (showValueLabels) {
        if (jitterValueLabels) {
            ggvaluelabels <-  geom_text(aes(y = ypos, label = sprintf("%.*f%%", labelDigits, 100 * prc)),
                                        vjust = jvert)
        } else {
            ggvaluelabels <-  geom_text(aes(y = ypos, label = sprintf("%.*f%%", labelDigits, 100 * prc)),
                                        vjust = vert)
        }
    } else {
        ggvaluelabels <-  geom_text(label = "")
    }
    # --------------------------------------------------------
    # Set up grid breaks
    # --------------------------------------------------------
    if (is.null(gridBreaksAt)) {
        gridbreaks <- waiver()
    } else {
        gridbreaks <- c(seq(0, 1, by = gridBreaksAt))
    }
    # --------------------------------------------------------
    # check if category-oder on x-axis should be reversed
    # change x axis order then
    # --------------------------------------------------------
    #geom.size <- c(1,geom.size)
    l <- length(unique(mydat$cat))
    w <- geom.size[seq(l +1, length(geom.size),l)]
    #w <- c(1.5,w)
    pos <- 0.5 * (cumsum(w) + cumsum(c(0, w[-length(w)])))
    pos <- pos + 1.4
    pos <- c(1,pos)
    mydat$grp <- as.numeric(sort(rep(pos,l)))


    if (reverseOrder && is.null(sort.frq)) {
        baseplot <- ggplot(mydat, aes(x = rev(grp), y = prc, fill = cat))
    } else {
        baseplot <- ggplot(mydat, aes(x = grp, y = prc, fill = cat))
    }
    baseplot <- baseplot +
        # plot bar chart
        geom_bar(aes(x = grp, y = prc, fill = cat),stat = "identity",
                 position = "stack", width = 0.98 *geom.size, colour="black")

    # --------------------------------------------------------
    # check whether bars should be visually separated by an
    # additional separator line
    # --------------------------------------------------------
    if (showSeparatorLine) {
        baseplot <- baseplot +
            geom_vline(xintercept = c(seq(1.5, length(items), by = 1)),
                       size = separatorLineSize,
                       colour = separatorLineColor)
    }
    # -----------------
    # show/hide percentage values on x axis
    # ----------------------------
    if (!showPercentageAxis) percent <- NULL
    baseplot <- baseplot +
        # show absolute and percentage value of each bar.
        ggvaluelabels +
        # no additional labels for the x- and y-axis, only diagram title
        labs(title = title, x = axisTitle.x, y = axisTitle.y, fill = legendTitle) +
        # print value labels to the x-axis.
        # If parameter "axisLabels.y" is NULL, the category numbers (1 to ...)
        # appear on the x-axis
        scale_x_continuous(labels = axisLabels.y, breaks = pos) +
        # set Y-axis, depending on the calculated upper y-range.
        # It either corresponds to the maximum amount of cases in the data set
        # (length of var) or to the highest count of var's categories.
        scale_y_continuous(breaks = gridbreaks,
                           limits = c(0, 1),
                           expand = expgrid,
                           labels = percent)
    # check whether coordinates should be flipped, i.e.
    # swap x and y axis
    if (coord.flip) baseplot <- baseplot + coord_flip()
    # ---------------------------------------------------------
    # set geom colors
    # ---------------------------------------------------------
    baseplot <- sjPlot:::sj.setGeomColors(baseplot,
                                          geom.colors,
                                          length(legendLabels),
                                          ifelse(isTRUE(hideLegend), FALSE, TRUE),
                                          legendLabels)
    # ---------------------------------------------------------
    # Check whether ggplot object should be returned or plotted
    # ---------------------------------------------------------
    if (printPlot) plot(baseplot)
    # -------------------------------------
    # return results
    # -------------------------------------
    invisible(structure(class = "sjpstackfrq",
                        list(plot = baseplot,
                             df = mydat)))
}


# unlist labels
# Help function that unlists a list into a vector
unlistlabels <- function(lab) {
    dummy <- unlist(lab)
    labels <- c()
    labels <- c(labels, as.character(dummy))
    return(labels)
}
