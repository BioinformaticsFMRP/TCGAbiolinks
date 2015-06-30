#' @title SurvivalCoxNET
#' @description Survival analysis with univariate Cox regression package (dnet)
#' @param clinical_patient clinical_patient
#' @param dataGE dataGE
#' @param Genelist Genelist
#' @param scoreConfidence restrict to those edges with high confidence (eg. score>=700)
#' @param titlePlot titlePlot
#' @importFrom survival coxph
#' @importFrom igraph subgraph.edges layout.fruchterman.reingold
#'             spinglass.community degree E communities crossing V
#' @importFrom dnet dRDataLoader dNetInduce dNetPipeline visNet dCommSignif
#' @importFrom supraHex visColormap visColoralpha
#' @importFrom grDevices dev.list
#' @export
#' @return net IGRAPH with attr: name (v/c), seqid (v/c), geneid (v/n), symbol (v/c), description (v/c) ...
SurvivalCoxNET <- function(clinical_patient,dataGE,Genelist,
                           scoreConfidence = 700,
                           titlePlot = "SurvivalCoxNET Example"){


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


    network <- subgraph.edges(org.Hs.string, eids=E(org.Hs.string)[combined_score>=scoreConfidence])
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


