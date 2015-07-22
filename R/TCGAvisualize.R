#' @title Survival analysis with univariate Cox regression package (dnet)
#' @description Survival analysis with univariate Cox regression and package (dnet)
#' @param clinical_patient clinical_patient
#' @param dataGE dataGE
#' @param Genelist Genelist
#' @param scoreConfidence restrict to those edges with high confidence (eg. score>=700)
#' @param titlePlot titlePlot
#' @importFrom survival coxph
#' @importFrom igraph subgraph.edges layout.fruchterman.reingold
#'             spinglass.community degree E communities crossing V V<-
#' @importFrom dnet dRDataLoader dNetInduce dNetPipeline
#'             visNet dCommSignif
#' @importFrom supraHex visColormap visColoralpha
#' @importFrom grDevices dev.list
#' @export
#' @return net IGRAPH with attr: name (v/c), seqid (v/c), geneid (v/n), symbol (v/c), description (v/c) ...
TCGAvisualize_SurvivalCoxNET <- function(clinical_patient,dataGE,Genelist,
                                         scoreConfidence = 700,
                                         titlePlot = "TCGAvisualize_SurvivalCoxNET Example"){

    combined_score <- NULL
    if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}

    png("SurvivalCoxNETOutput.png", width = 800, height = 800)

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
    org.Hs.string <- dRDataLoader(RData='org.Hs.string')
    # restrict to those edges with high confidence (score>=700)
    with(org.Hs.string,{
        network <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=scoreConfidence])})
    network

    # extract network that only contains genes in pvals
    ind <- match(V(network)$symbol, names(pvals))
    ## for extracted graph
    nodes_mapped <- V(network)$name[!is.na(ind)]
    network <- dNetInduce(g=network, nodes_query=nodes_mapped, knn=0,
                          remove.loops=FALSE, largest.comp=TRUE)
    V(network)$name <- V(network)$symbol
    network

    # Identification of gene-active network
    net <- dNetPipeline(g=network, pval=pvals, method="customised",
                        significance.threshold=5e-02)
    net

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
    legend("topright", legend=legend_name, fill=mcolors, bty="n", cex=1.4)

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


