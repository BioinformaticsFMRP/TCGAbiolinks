#' @title Calculate diffmean methylation between two groups
#' @description
#'    Calculate diffmean methylation of probes between two groups removing lines
#'    that has NA values.
#' @param data SummarizedExperiment object obtained from TCGAPrepare
#' @param groupCol Columns in colData(data) that defines the groups.
#' @param group1 Name of group1 to be used in the analysis
#' @param group2 Name of group2  to be used in the analysis
#' @import ggplot2
#' @import graphics
#' @importFrom grDevices png dev.off
#' @importFrom S4Vectors values
#' @importFrom SummarizedExperiment colData rowRanges assay rowRanges<- values<-
#' @return Saves in the rowRages(data) the columns: mean.group1, mean.group2
#'        diffmean.group1.group2; Where group1 and group2 are the names of the
#'        groups.
#' @examples
#' \dontrun{
#' nrows <- 200; ncols <- 20
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                    IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                     strand=sample(c("+", "-"), 200, TRUE),
#'                     feature_id=sprintf("ID%03d", 1:200))
#'colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 10),
#'                     row.names=LETTERS[1:20],
#'                     group=rep(c("group1","group2"),c(10,10)))
#'data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=counts),
#'          rowRanges=rowRanges,
#'          colData=colData)
#' data <- diffmean(data)
#' }
diffmean <- function(data, groupCol = NULL, group1 = NULL, group2 = NULL) {

    if (is.null(groupCol)) {
        message("Please, set the groupCol parameter")
        return(NULL)
    }
    if ( length(unique(colData(data)[,groupCol])) != 2 &&
         is.null(group1) && is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    } else if (length(unique(colData(data)[,groupCol])) == 2) {
        group1 <- unique(colData(data)[,groupCol])[1]
        group2 <- unique(colData(data)[,groupCol])[2]
    }
    message("Calculating the diference between the mean methylation of the groups...")

    m <- assay(data)
    idx1 <- which(colData(data)[,groupCol] == group1)
    idx2 <- which(colData(data)[,groupCol] == group2)
    mean.g1 <- rowMeans(m[,idx1], na.rm = TRUE)
    mean.g2 <- rowMeans(m[,idx2], na.rm = TRUE)
    diffmean <- mean.g1 - mean.g2

    # Saves the result
    values(rowRanges(data))[,paste0("mean.",group1)] <-  mean.g1
    values(rowRanges(data))[,paste0("mean.",group2)] <-  mean.g2
    values(rowRanges(data))[,paste0("diffmean.",group1,".",group2)] <-  diffmean

    # Ploting a histogram to evaluate the data
    message("Saved histogram_diffmean.png...")
    png(filename = "histogram_diffmean.png")
    hist(diffmean)
    dev.off()

    return(data)
}

#' @title Creates survival analysis
#' @description Creates a survival plot from TCGA patient clinical data
#' using survival library. It uses the fields days_to_death and vital, plus a
#' columns for groups.
#'
#' @param data TCGA Clinical patient with the information days_to_death
#' @param clusterCol Column with groups to plot. This is a mandatory field, the
#' caption will be based in this column
#' @param legend Legend title of the figure
#' @param cutoff xlim This parameter will be a limit in the x-axis. That means, that
#' patients with days_to_deth > cutoff will be set to Alive.
#' @param main main title of the plot
#' @param ylab y axis text of the plot
#' @param xlab x axis text of the plot
#' @param filename The name of the pdf file
#' @param color Define the colors of the lines.
#' @param width Image width
#' @param height Image height
#' @importFrom GGally ggsurv
#' @importFrom survival survfit Surv
#' @importFrom scales percent
#' @export
#' @return Survival plot
#' @examples
#' days_to_death <- floor(runif(200, 1, 1000))
#' vital_status <- c(rep("Dead",200))
#' groups <- c(rep(c("G1","G2"),c(100,100)))
#' df <- data.frame(days_to_death,vital_status,groups)
#' survivalAnalysis(df,clusterCol="groups")
#' \dontrun{
#' clinical <- clinic("gbm","clinical_patient")
#' survivalAnalysis(clinical,"gender", filename = "surv.pdf", legend="Gender")
#' }
survivalAnalysis <- function(data,
                             clusterCol=NULL,
                             legend = "Legend", cutoff = 0,
                             main = "Kaplan-Meier Overall Survival Curves",
                             ylab = "PROBABILITY OF SURVIVAL",
                             xlab = "TIME SINCE DIAGNOSIS (DAYS)",
                             filename = "survival.pdf",
                             color = c("green", "firebrick4", "orange3", "blue"),
                             height=10,
                             width=5
) {
    .e <- environment()
    group <- NULL
    if (is.null(clusterCol)) {
        message("Please provide the clusterCol argument")
        return(NULL)
    }
    notDead <- which(data$days_to_death == "[Not Applicable]")

    if (length(notDead) > 0) {
        data[notDead,]$days_to_death <- data[notDead,]$days_to_last_followup
    }

    if (cutoff != 0) {
        # Axis-x cut-off
        aux <- subset(data, data$days_to_death > cutoff)
        # 'cut' info bigger than cutoff
        data[rownames(aux), "days_to_death"] <- cutoff
        # Pacient is alive (0:alive,1:dead)
        data[rownames(aux), "vital_status"] <- "Alive"
    }
    # create a column to be used with survival package, info need
    # to be TRUE(DEAD)/FALSE (ALIVE)
    data$s <- !(data$vital_status == "Dead")

    # Column with groups
    data$type <- as.factor(data[,clusterCol])
    # create the formula for survival analysis
    f.m <- formula(Surv(as.numeric(data$days_to_death),
                        data$s) ~ data$type)
    fit <- survfit(f.m, data = data)

    surv <- ggsurv(fit, CI = "def", plot.cens = FALSE,
                   surv.col = "gg.def",
                   cens.col = "red", lty.est = 1,
                   lty.ci = 2, cens.shape = 3,
                   back.white = TRUE,
                   xlab = xlab, ylab = ylab, main = main)
    if (cutoff != 0) {
        surv <- surv + ggplot2::coord_cartesian(xlim = c(0, cutoff))
    }
    label.add.n <- function(x) {
        paste0(x, " (n = ",
               nrow(subset(data,data[,clusterCol] == x)), ")")
    }
    with(data,{
        surv <- surv + scale_colour_discrete(name = legend,
                                             labels = sapply(levels(data$type),label.add.n)
        )
        with(surv,{
            surv <- surv + geom_point(aes(colour = group),
                                      shape = 3,size = 2)
            surv <- surv + guides(linetype = FALSE) +
                scale_y_continuous(labels = scales::percent)
            ggsave(surv, filename = filename,
                   width = width, height = height)
        })
    })

}
#' @title Mean methylation boxplot
#' @description
#'   Creates a mean methylation boxplot divided in by groups
#' @param data SummarizedExperiment object obtained from TCGAPrepare
#' @param groupCol Columns in colData(data) that defines the groups. If no
#' columns defined a columns called "Patients" will be used
#' @param filename The name of the pdf that will be saved
#' @param legend Caption title
#' @param color vector of colors to be used in graph
#' @param title main title in the plot
#' @param ylab y axis text in the plot
#' @param xlab x axis text in the plot
#' @param sort Sort by mean methylation? False by default
#' @import ggplot2 stats
#' @importFrom SummarizedExperiment colData rowRanges assay
#' @export
#' @return Save the pdf survival plot
#' @examples
#' nrows <- 200; ncols <- 21
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                    IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                     strand=sample(c("+", "-"), 200, TRUE),
#'                     feature_id=sprintf("ID%03d", 1:200))
#'colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input","Other"), 7),
#'                     row.names=LETTERS[1:21],
#'                     group=rep(c("group1","group2","group3"),c(7,7,7)))
#'data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=counts),
#'          rowRanges=rowRanges,
#'          colData=colData)
#' meanMethylationAnalysis(data,groupCol  = "group",sort=TRUE)
meanMethylationAnalysis <- function(data,
                                    groupCol=NULL,
                                    sort = FALSE,
                                    filename = "G-CIMP-mean.methylation.pdf",
                                    ylab = "Mean DNA methylation",
                                    xlab = "DNA Methylation Clusters",
                                    title = "Mean DNA methylation by cluster",
                                    legend = "Legend",
                                    color = c("green", "red", "purple",
                                              "orange", "salmon", "grey")) {
    .e <- environment()
    mean <- apply(assay(data), 2, mean,na.rm = TRUE)

    if (is.null(groupCol)){
        groups <- rep("Patient",length(mean))
    } else {
        groups <- colData(data)[,groupCol]
    }
    df <- data.frame(mean = mean, groups = groups)

    # Plot for methylation analysis Axis x: LGm clusters Axis y:
    # mean methylation
    if (sort) {
        p <- ggplot(df, aes(reorder(factor(df$groups), df$mean),
                            df$mean),  environment = .e) +
            geom_boxplot(aes(fill = reorder(factor(df$groups),
                                            df$mean)),
                         notchwidth = 0.25) +
            geom_jitter(height = 0,position = position_jitter(width = 0.1),
                        size = 3) +
            scale_fill_manual(values = color,
                              labels = levels(reorder(factor(df$groups),
                                                      df$mean)),
                              name = legend)
    } else {
        p <- ggplot(df, aes(factor(df$groups), df$mean),
                    environment = .e) +
            geom_boxplot(aes(fill = factor(df$groups)),
                         notchwidth = 0.25) +
            geom_jitter(height = 0,
                        position = position_jitter(width = 0.1),
                        size = 3) +
            scale_fill_manual(values = color,
                              labels = levels(factor(df$groups)),
                              name = legend)
    }
    p <- p + ylab(ylab) + xlab(xlab) + labs(title = title) +
        theme(axis.title.x = element_text(face = "bold", size = 20),
              axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         size = 16),
              axis.title.y = element_text(face = "bold",
                                          size = 20),
              axis.text.y = element_text(size = 16),
              plot.title = element_text(face = "bold", size = 16))

    # saving box plot to analyse it
    ggsave(p, filename = filename, width = 10, height = 10, dpi = 600)
    message(paste("Plot saved in: ", file.path(getwd(),filename)))
}

#' @title Calculate pvalues
#' @details
#'    Verify if the data is significant between two groups. For the methylation
#'    we search for probes that have a difference in the mean methylation and
#'    also a significant value.
#'    Input: A SummarizedExperiment object that will be used to
#'    compared two groups with wilcoxon test, a boolean value to do a
#'    paired or non-paired test
#'    Output: p-values (non-adj/adj) histograms, p-values (non-adj/adj)
#' @param data  SummarizedExperiment obtained from the TCGAPrepare
#' @param groupCol  Columns with the groups inside the SummarizedExperiment
#'  object. (This will be obtained by the function colData(data))
#' @param group1 In case our object has more than 2 groups, you should set the
#'  groups
#' @param group2 In case our object has more than 2 groups, you should set the
#'  groups
#' @param paired  Do a paired wilcoxon test? Default: True
#' @param exact  Do a exact wilcoxon test? Default: True
#' @param  method P-value adjustment method. Default:"BH" Benjamini-Hochberg
#' @return Data frame with cols p values/p values adjusted
#' @import graphics
#' @importFrom grDevices png dev.off pdf
#' @import stats
#' @importFrom coin wilcox_test wilcoxsign_test pvalue
#' @importFrom SummarizedExperiment colData rowRanges rowRanges<- colData<-
#' @return Data frame with two cols
#'         p-values/p-values adjusted
#' @examples
#' \dontrun{
#' nrows <- 200; ncols <- 20
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                    IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                     strand=sample(c("+", "-"), 200, TRUE),
#'                     feature_id=sprintf("ID%03d", 1:200))
#'colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 10),
#'                     row.names=LETTERS[1:20],
#'                     group=rep(c("group1","group2"),c(10,10)))
#'data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=counts),
#'          rowRanges=rowRanges,
#'          colData=colData)
#' data <- calculate.pvalues(data,"group")
#' }
#' @keywords internal
calculate.pvalues <- function(data,
                              groupCol = NULL,
                              group1 = NULL,
                              group2 = NULL,
                              paired = FALSE,
                              method = "BH",
                              exact = TRUE) {

    if (is.null(groupCol)) {
        message("Please, set the groupCol parameter")
        return(NULL)
    }
    if ( length(unique(colData(data)[,groupCol])) != 2 &&
         is.null(group1) && is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    } else if (length(unique(colData(data)[,groupCol])) == 2) {
        group1 <- unique(colData(data)[,groupCol])[1]
        group2 <- unique(colData(data)[,groupCol])[2]
    }
    message("Calculating the p-values of each probe...")
    # Apply Wilcoxon test in order to calculate the p-values
    idx1 <- which(colData(data)[,groupCol] == group1)
    idx2 <- which(colData(data)[,groupCol] == group2)

    if(!paired){
        p.value <- apply(assay(data),1,
                         function(x) {
                             if(!is.factor(colData(data)[,groupCol])) {
                                 colData(data)[,groupCol] <- factor(
                                     colData(data)[,groupCol]
                                     )
                             }
                             aux <-data.frame(beta=x[c(idx1,idx2)],
                                              cluster=droplevels(colData(data)[c(idx1,idx2),groupCol]))
                             pvalue(wilcox_test(beta ~ cluster, data=aux, distribution = "exact"))
                         }
        )
    } else {
        p.value <- apply(assay(data),1,
                         function(x) {
                             aux <-data.frame(beta=x[c(idx1,idx2)],
                                              cluster=droplevels(colData(data)[c(idx1,idx2),groupCol]))
                             pvalue(wilcoxsign_test(beta ~ cluster, data=aux, distribution = exact()))
                         }
        )
    }
    ## Plot a histogram
    message("Saved histogram_pvalues.png...")
    png(filename = "histogram_pvalues.png")
    hist(p.value)
    dev.off()

    ## Calculate the adjusted p-values by using Benjamini-Hochberg
    ## (BH) method
    p.value.adj <- p.adjust(p.value, method = method)

    ## Plot a histogram
    message("Saved histogram_pvalues_adj.png")
    png(filename = "histogram_pvalues_adj.png")
    hist(p.value.adj)
    dev.off()

    #Saving the values into the object
    colp <- paste("p.value", group1, group2, sep = ".")
    values(rowRanges(data))[,colp] <-  p.value
    coladj <- paste("p.value.adj",group1,group2, sep = ".")
    values(rowRanges(data))[,coladj] <-  p.value.adj

    return(data)
}

#' @title Differentially methylated regions Analysis
#' @description
#'   This function will search for differentially methylated CpG sites,
#'   which are regarded as possible functional regions involved
#'   in gene transcriptional regulation.
#'   In order to find these regions we use the beta-values (methylation values
#'   ranging from 0.0 to 1.0) to compare two groups.
#'   Firstly, it calculates the difference between the mean methylation of each
#'   group for each probes. Secondly, it calculates the p-value using the
#'   wilcoxon test using the Benjamini-Hochberg adjustment method.
#'   The default parameters will require a minimum absolute beta values delta
#'   of 0.2 and a false discovery rate (FDR)-adjusted Wilcoxon rank-sum P-value
#'   of <0.01 for the difference.
#'   After these analysis, we save a volcano plot (x-axis:diff mean methylation,
#'   y-axis: significance) that will help the user identify the differentially
#'   methylated CpG sites and return the object with the calculus in the rowRanges.
#' @param data  SummarizedExperiment obtained from the TCGAPrepare
#' @param groupCol  Columns with the groups inside the SummarizedExperiment
#'  object. (This will be obtained by the function colData(data))
#' @param group1 In case our object has more than 2 groups, you should set
#' the name of the group
#' @param group2 In case our object has more than 2 groups, you should set
#' the name of the group
#' @param filename pdf filename
#' @param legend legend title
#' @param color vector of colors to be used in graph
#' @param title main title
#' @param ylab y axis text
#' @param xlab x axis text
#' @param xlim x limits to cut image
#' @param ylim y limits to cut image
#' @param label vector of labels to be used in the figure
#' @param p.cut p values threshold
#' @param diffmean.cut diffmean threshold
#' @param adj.method Adjusted method for the p-value calculation
#' @param paired Wilcoxon paired parameter
#' @param overwrite Overwrite the pvalues and diffmean values if already in the object
#' for both groups? Default: FALSE
#' @import ggplot2
#' @importFrom SummarizedExperiment colData rowRanges assay rowRanges<- values<-
#' @export
#' @return Volcano plot saved and the given data with the results
#' (diffmean.group1.group2,p.value.group1.group2,
#' p.value.adj.group1.group2,status.group1.group2)
#' in the rowRanges where group1 and group2 are the names of the groups
#' @examples
#' nrows <- 200; ncols <- 20
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                    IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                     strand=sample(c("+", "-"), 200, TRUE),
#'                     feature_id=sprintf("ID%03d", 1:200))
#'colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
#'                     row.names=LETTERS[1:20],
#'                     group=rep(c("group1","group2"),c(10,10)))
#'data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=counts),
#'          rowRanges=rowRanges,
#'          colData=colData)
#' SummarizedExperiment::colData(data)$group <- c(rep("group1",ncol(data)/2),
#'                          rep("group2",ncol(data)/2))
#' hypo.hyper <- DMRAnalysis(data, p.cut = 0.85,"group","group1","group2")
DMRAnalysis <- function(data,
                        groupCol=NULL,
                        group1=NULL,
                        group2=NULL,
                        filename = "volcano.pdf",
                        ylab =  expression(paste(-Log[10],
                                                 " (FDR corrected -P values)")),
                        xlab = "DNA Methylation\n difference",
                        title = "Volcano plot",
                        legend = "Legend",
                        color = c("1" = "black", "2" = "green",
                                  "3" = "red"),
                        label = c("1" = "Not Significant",
                                  "2" = "Hypermethylated",
                                  "3" = "Hypomethylated"),
                        xlim = NULL,
                        ylim = NULL,
                        p.cut = 0.01,
                        diffmean.cut = 0.2,
                        paired = FALSE,
                        adj.method="BH",
                        overwrite=FALSE) {
    .e <- environment()

    if (is.null(groupCol)) {
        message("Please, set the groupCol parameter")
        return(NULL)
    }
    if ( length(unique(colData(data)[,groupCol])) != 2 &&
         is.null(group1) && is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    } else if (length(unique(colData(data)[,groupCol])) == 2) {
        group1 <- unique(colData(data)[,groupCol])[1]
        group2 <- unique(colData(data)[,groupCol])[2]
    }

    diffcol <- paste("diffmean",group1,group2,sep = ".")
    if (!(diffcol %in% colnames(values(rowRanges(data)))) || overwrite) {
        data <- diffmean(data,groupCol, group1 = group1, group2 = group2)
        if (!(diffcol %in% colnames(values(rowRanges(data))))) stop("Error!")
    }
    pcol <- paste("p.value.adj",group1,group2,sep = ".")
    if (!(pcol %in% colnames(values(rowRanges(data)))) || overwrite) {
        data <- calculate.pvalues(data,groupCol, group1, group2,
                                  paired = paired,method = adj.method)
        # An error should not happen, if it happens (probably due to an incorret
        # user input) we will stop
        if (!(pcol %in% colnames(values(rowRanges(data))))) stop("Error!")
    }

    statuscol <- paste("status",group1,group2,sep = ".")
    values(rowRanges(data))[,statuscol] <-  "Not Significant"
    rowRanges(data)$threshold <- "1"

    # get significant data
    sig <-  values(rowRanges(data))[,pcol] < p.cut

    # hypermethylated samples compared to old state
    hyper <- values(rowRanges(data))[,diffcol] < (-diffmean.cut)

    if (any(hyper & sig)) rowRanges(data)[hyper & sig,]$threshold <- "2"
    if (any(hyper & sig)) values(rowRanges(data))[hyper & sig,statuscol] <- "Hypermethylated"
    # hypomethylated samples compared to old state
    hypo <-  values(rowRanges(data))[,diffcol]  > diffmean.cut
    if (any(hypo & sig)) rowRanges(data)[hypo & sig,]$threshold <- "3"
    if (any(hypo & sig)) values(rowRanges(data))[hypo & sig,statuscol] <- "Hypomethylated"

    # Plot a volcano plot
    p <- ggplot(data = as.data.frame(rowRanges(data)),
                aes(x = values(rowRanges(data))[,diffcol] ,
                    y = -1 * log10(values(rowRanges(data))[,pcol]),
                    colour = rowRanges(data)$threshold ),
                environment = .e) + geom_point() +
        labs(title = title) + ylab(ylab) + xlab(xlab) +
        geom_vline(aes(xintercept = -diffmean.cut),
                   colour = "black",linetype = "dashed") +
        geom_vline(aes(xintercept = diffmean.cut),
                   colour = "black", linetype = "dashed") +
        geom_hline(aes(yintercept = -1 * log10(p.cut)),
                   colour = "black", linetype = "dashed") +
        scale_color_manual(breaks = c("1", "2", "3"),
                           values = color,
                           labels = label,
                           name = legend)
    # saving box plot to analyse it
    ggsave(p, filename = filename, width = 10, height = 5, dpi = 600)
    rowRanges(data)$threshold <- NULL
    return(data)
}

# Get latest Genome Reference Consortium Human Build And save
# it as Genomic Ranges
#' @importFrom biomaRt useMart getBM
get.GRCh.bioMart <- function(genome="hg19") {

    if (genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "grch37.ensembl.org",
                           path = "/biomart/martservice" ,
                           dataset = "hsapiens_gene_ensembl")
    } else {
        # for hg39
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    }

    chrom <- c(1:22, "X", "Y")
    gene.location <- getBM(attributes = c("chromosome_name",
                                          "start_position",
                                          "end_position", "strand",
                                          "external_gene_name",
                                          "entrezgene"),
                           filters = c("chromosome_name"),
                           values = list(chrom), mart = ensembl)

    gene.location <- gene.location[!is.na(gene.location$entrezgene),]
    # the duplicates are different transcripts, not different
    # coordinates
    gene.location <- gene.location[!duplicated(gene.location$entrezgene),]
    save(gene.location, file = "inst/extdata/GRCh.rda")
}

#' @title Create starburst plot
#'
#' @description
#'   Create Starburst plot for comparison of DNA methylation and gene expression.
#'    The log10 (FDR-corrected P value) is plotted for beta value for DNA
#'    methylation (x axis) and gene expression (y axis) for each gene.
#'    The black dashed line shows the FDR-adjusted P value of 0.01.
#'
#' @details
#'    Input: data with gene expression/methylation expression
#'    Output: starburst plot
#'
#' @param met SummarizedExperiment with methylation data obtained from the
#' TCGAPrepare. Expected colData columns: diffmean,  p.value.adj  and p.value
#' Execute volcanoPlot function in order to obtain these values for the object.
#' @param exp Object obtained by DEArnaSEQ function
#' @param filename pdf filename
#' @param legend legend title
#' @param color vector of colors to be used in graph
#' @param label vector of labels to be used in graph
#' @param title main title
#' @param ylab y axis text
#' @param xlab x axis text
#' @param xlim x limits to cut image
#' @param ylim y limits to cut image
#' @param p.cut p value cut-off
#' @param group1 In case our object has more than 2 groups, you should set
#' the name of the group
#' @param group2 In case our object has more than 2 groups, you should set
#' the name of the group
#' @import ggplot2
#' @importFrom SummarizedExperiment subsetByOverlaps rowRanges rowRanges<-
#'             values<-
#' @export
#' @return Save a starburst plot
#' @examples
#' nrows <- 20000; ncols <- 20
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' ranges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(5000, 15000)),
#'                    IRanges::IRanges(floor(runif(20000, 1e5, 1e6)), width=100),
#'                     strand=sample(c("+", "-"), 20000, TRUE),
#'                     probeID=sprintf("ID%03d", 1:20000),
#'                     Gene_Symbol=sprintf("ID%03d", 1:20000))
#'colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
#'                     row.names=LETTERS[1:20],
#'                     group=rep(c("group1","group2"),c(10,10)))
#'data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=counts),
#'          rowRanges=ranges,
#'          colData=colData)
#' met <- data
#' exp <- data.frame(row.names=sprintf("ID%03d", 1:20000),
#'                   logFC=runif(20000, -0.2, 0.2),
#'                   FDR=runif(20000, 0.01, 1))
#' SummarizedExperiment::rowRanges(met)$diffmean.g1.g2 <- c(runif(20000, -0.1, 0.1))
#' SummarizedExperiment::rowRanges(met)$p.value.g1.g2 <- c(runif(20000, 0, 1))
#' SummarizedExperiment::rowRanges(met)$p.value.adj.g1.g2 <- c(runif(20000, 0, 1))
#' result <- starburstAnalysis(met,exp,p.cut = 0.05,"g1","g2")
starburstAnalysis <- function(met,
                              exp,
                              group1=NULL,
                              group2=NULL,
                              filename = "volcano.pdf",
                              ylab = expression(atop("Gene Expression",
                                                     paste(Log[10],
                                                           " (FDR corrected P values)"))),
                              xlab = expression(atop("DNA Methylation",
                                                     paste(Log[10],
                                                           " (FDR corrected P values)"))),
                              title = "Starburst Plot",
                              legend = "Methylation/Expression Relation",
                              color = c("1" = "black",
                                        "2" = "purple",
                                        "3" = "darkgreen",
                                        "4" = "blue",
                                        "5" = "darkred",
                                        "6" = "red",
                                        "7" = "green",
                                        "8" = "yellow",
                                        "9" = "orange"),
                              label = c("1" = "Not Significant",
                                        "2" = "Up regulated & Hypo methylated",
                                        "3" = "Down regulated & Hypo methylated",
                                        "4" = "hypo methylated",
                                        "5" = "hyper methylated",
                                        "6" = "Up regulated",
                                        "7" = "Down regulated",
                                        "8" = "Up regulated & Hyper methylated",
                                        "9" = "Down regulated & Hyper methylated"),
                              xlim = NULL, ylim = NULL, p.cut = 0.01
)
{
    .e <- environment()

    if ( is.null(group1) || is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    }

    met <- as.data.frame(rowRanges(met))
    exp$Gene_Symbol <- row.names(exp)
    volcano <- merge(met, exp, by = "Gene_Symbol")

    volcano$ID <- paste(volcano$Gene_Symbol,
                        volcano$probeID, sep = ".")

    # Preparing gene expression
    volcano$geFDR <- log10(volcano$FDR)
    volcano$geFDR2 <- volcano$geFDR
    volcano[volcano$logFC > 0, "geFDR2"] <-
        -1 * volcano[volcano$logFC > 0, "geFDR"]

    # Preparing methylation
    pcol <- paste("p.value.adj",group1,group2,sep = ".")
    diffcol <- paste("diffmean",group1,group2,sep = ".")
    volcano$meFDR <- log10(volcano[,pcol])
    volcano$meFDR2 <- volcano$meFDR
    volcano[volcano[,diffcol] > 0, "meFDR2"] <-
        -1 * volcano[volcano[,diffcol] > 0, "meFDR"]

    volcano$threshold.starburst <- "1"
    volcano$threshold.size <- "1"

    # subseting by regulation (geFDR) and methylation level
    # (meFDR) down regulated up regulated lowerthr
    # |||||||||||||||| upperthr hypomethylated hipermethylated
    lowerthr <- log10(p.cut)  # - 1.30103
    upperthr <- (-lowerthr)  # +1.30103

    # Group 2:up regulated and hypomethylated
    a <- subset(volcano,
                volcano$geFDR2 > upperthr &
                    volcano$meFDR2 < lowerthr)

    # Group 3: down regulated and hypomethylated
    b <- subset(volcano,
                volcano$geFDR2 < lowerthr &
                    volcano$meFDR2 < lowerthr)

    # Group 4: hypomethylated
    c <- subset(volcano,
                volcano$geFDR2 > lowerthr &
                    volcano$geFDR2 < upperthr &
                    volcano$meFDR2 < lowerthr)

    # Group 5: hypermethylated
    d <- subset(volcano,
                volcano$geFDR2 > lowerthr &
                    volcano$geFDR2 < upperthr &
                    volcano$meFDR2 > upperthr)

    # Group 6: upregulated
    e <- subset(volcano,
                volcano$geFDR2 > upperthr &
                    volcano$meFDR2 < upperthr &
                    volcano$meFDR2 > lowerthr)

    # Group 7: downregulated
    f <- subset(volcano,
                volcano$geFDR2 < lowerthr &
                    volcano$meFDR2 < upperthr &
                    volcano$meFDR2 > lowerthr)

    # Group 8: upregulated and hypermethylated
    g <- subset(volcano,
                volcano$geFDR2 > upperthr &
                    volcano$meFDR2 > upperthr)

    # Group 9: downregulated and hypermethylated
    h <- subset(volcano,
                volcano$geFDR2 < lowerthr &
                    volcano$meFDR2 > upperthr)

    size <- c("1", "1", "1", "1", "1", "1", "1","1")
    groups <- c("2", "3", "4", "5", "6", "7","8","9")
    s <- list(a, b, c, d, e, f,g,h)
    for (i in seq_along(s)) {
        idx <- rownames(s[[i]])
        if (length(idx) > 0) {
            volcano[idx, "threshold.starburst"] <- groups[i]
            volcano[idx, "threshold.size"] <- size[i]
        }
    }

    ## starburst plot
    p <- ggplot(data = volcano, environment = .e,
                aes(x = volcano$meFDR2,
                    y = volcano$geFDR2,
                    colour = volcano$threshold.starburst,
                    size = volcano$threshold.size)) +
        geom_point()
    if (!is.null(xlim)) {
        p <- p + xlim(xlim)
    }
    if (!is.null(ylim)) {
        p <- p + ylim(ylim)
    }
    p <- p + labs(title) + ylab(ylab) + xlab(xlab)
    p <- p + scale_color_manual(values = color, labels = label, name = legend)
    p <-  p + geom_hline(aes(yintercept = lowerthr), colour = "black",
                         linetype = "dashed") +
        geom_hline(aes(yintercept = upperthr), colour = "black",
                   linetype = "dashed") +
        geom_vline(aes(xintercept = lowerthr), colour = "black",
                   linetype = "dashed") +
        geom_vline(aes(xintercept = upperthr), colour = "black",
                   linetype = "dashed")
    ggsave(filename = "starbust.gcimp.pdf", width = 14, height = 10)

    # return methylation < 0, expressao >0
    return (volcano)
}
