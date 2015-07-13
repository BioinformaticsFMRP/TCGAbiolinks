#' @title Calculate diffmean methylation between two groups
#' @description
#'    Calculate diffmean methylation between two groups
#' @param data SummarizedExperiment object obtained from TCGAPrepare
#' @param groupCol Columns in colData(data) that defines the groups.
#' @param group2 Name of group2 to be used in the analysis
#' @param group1 Name of group1  to be used in the analysis
#' @import ggplot2
#' @import graphics
#' @importFrom grDevices png dev.off
#' @importFrom SummarizedExperiment colData rowRanges assay rowRanges<-
#' @return dataframe with diffmean values
#' @examples
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
diffmean <- function(data, groupCol = NULL, group2 = NULL, group1 = NULL) {

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
    mean.g1 <- apply(m[,idx1], 1, mean, na.rm = TRUE)
    mean.g2 <- apply(m[,idx2], 1, mean, na.rm = TRUE)
    diffmean <- mean.g1 - mean.g2

    # Saves the result into diffmean column
    rowRanges(data)$diffmean <-  diffmean

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
#' @param clinical_patient TCGA Clinical patient with the information days_to_death
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
#' @importFrom GGally ggsurv
#' @importFrom survival survfit Surv
#' @importFrom scales percent
#' @export
#' @return Survival plot
#' @examples
#' clinical <- clinic("gbm","clinical_patient")
#' survivalPlot(clinical,"gender", filename = "surv.pdf", legend="Gender")
#' survivalPlot(met.md)
survivalPlot <- function(clinical_patient,
                         clusterCol=NULL,
                         legend = "Legend", cutoff = 0,
                         main = "Kaplan-Meier Overall Survival Curves",
                         ylab = "PROBABILITY OF SURVIVAL",
                         xlab = "TIME SINCE DIAGNOSIS (DAYS)",
                         filename = "survival.pdf",
                         color = c("green", "firebrick4", "orange3", "blue")
) {
    .e <- environment()
    group <- NULL
    if (is.null(clusterCol)) {
        message("Please provide the clusterCol argument")
        return(NULL)
    }
    notDead <- which(clinical_patient$days_to_death == "[Not Applicable]")

    if (length(notDead) > 0) {
        clinical_patient[notDead,]$days_to_death <- Inf
    }

    if (cutoff != 0) {
        # Axis-x cut-off
        aux <- subset(clinical_patient, clinical_patient$days_to_death > cutoff)
        # 'cut' info bigger than cutoff
        clinical_patient[rownames(aux), "days_to_death"] <- cutoff
        # Pacient is alive (0:alive,1:dead)
        clinical_patient[rownames(aux), "vital_status"] <- "Alive"
    }
    # create a column to be used with survival package, info need
    # to be TRUE(DEAD)/FALSE (ALIVE)
    clinical_patient$s <- clinical_patient$vital == "Dead"

    # Column with groups
    clinical_patient$type <- as.factor(clinical_patient[,clusterCol])

    # create the formula for survival analysis
    f.m <- formula(Surv(as.numeric(clinical_patient$days_to_death),
                        clinical_patient$s) ~ clinical_patient$type)
    fit <- survfit(f.m, data = clinical_patient)

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
               nrow(subset(clinical_patient,clinical_patient[,clusterCol] == x)), ")")
    }
    with(clinical_patient,{
        surv <- surv + scale_colour_discrete(name = legend,
                    labels = sapply(levels(clinical_patient$type),label.add.n)
        )
        with(surv,{
            surv <- surv + geom_point(aes(colour = group),
                                      shape = 3,size = 2)
            surv <- surv + guides(linetype = FALSE) +
                scale_y_continuous(labels = scales::percent)
            ggsave(surv, filename = filename)
        })
    })

}
#' @title Mean methylation boxplot
#' @description
#'   Creates a mean methylation boxplot based for patients divided in by groups
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
#' metMeanBoxplot(data,groupCol  = "group",sort=TRUE)
metMeanBoxplot <- function(data,
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
#' @importFrom exactRankTests wilcox.exact
#' @import graphics
#' @importFrom grDevices png dev.off pdf
#' @import stats
#' @importFrom SummarizedExperiment colData rowRanges rowRanges<-
#' @return Data frame with two cols
#'         p-values/p-values adjusted
#' @examples
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
#' data <- calculate.pvalues(data)
#' @keywords internal
calculate.pvalues <- function(data,
                              groupCol = NULL,
                              group1 = NULL,
                              group2 = NULL,
                              paired = TRUE,
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
    print(idx1)
    print(idx2)
    p.value <- apply(assay(data),1,
                    function(x) {
                        print(x[idx1])
                        wilcox.test(x[idx1], x[idx2],
                        exact = exact, paired = paired)$p.value
                        }
    )

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
    rowRanges(data)$p.value <-  p.value
    rowRanges(data)$p.value.adj <-  p.value.adj

    return(data)
}

#' @title Volcano plot
#' @description
#'   In order to searches for the probes that are different methylated and
#'   are significant, we use a volcano plot (x-axis:diff mean methylation,
#'   y-axis: significance) that compares the methylation data between the two groups.
#'   Firstly, it calculates the difference between the mean methylation of each group
#' for each probes. After, it calculates the p-value using the wilcoxon test
#' using the Benjamini-Hochberg adjustment method. With both values, it is possible
#' to analyse the data.
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
#' @import ggplot2
#' @importFrom SummarizedExperiment colData rowRanges assay rowRanges<-
#' @export
#' @return A dataframe with the Composite.Element.REF and
#'         the group it was classified
#'          group 1 = Not Significant
#'          group 2 = Hypermethylated
#'          group 3 = Hypomethylated
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
#' hypo.hyper <- volcanoPlot(data, p.cut = 0.85)
volcanoPlot <- function(data,
                        groupCol=NULL,
                        group1=NULL,
                        group2=NULL,
                        filename = "volcano.pdf",
                        ylab = "- 1*log10 of the Significance",
                        xlab = "DNA Methylation",
                        title = "Volcano Plot TCGA GBM Tumors",
                        legend = "Legend",
                        color = c("1" = "black", "2" = "green",
                                  "3" = "red"),
                        label = c("1" = "Not Significant",
                                  "2" = "Hypermethylated",
                                  "3" = "Hypomethylated"),
                        xlim = NULL,
                        ylim = NULL,
                        p.cut = 0.05,
                        diffmean.cut = 0) {
    .e <- environment()


    data <- calculate.pvalues(data,groupCol, group1, group2)
    data <- diffmean(data,groupCol, group1 = group1, group2 = group2)
    if (is.null(data$diffmean)) return(NULL)

    rowRanges(data)$threshold <- "1"

    # get significant data
    sig <-  rowRanges(data)$p.value.adj < p.cut

    # hypermethylated samples compared to old state
    hyper <- rowRanges(data)$diffmean < (-diffmean.cut)
    if (any(hyper & sig)) rowRanges(data)[hyper & sig,]$threshold <- "2"

    # hypomethylated samples compared to old state
    hypo <-  rowRanges(data)$diffmean  > diffmean.cut
    if (any(hypo & sig))rowRanges(data)[hypo & sig,]$threshold <- "3"

    # Plot a volcano plot
    p <- ggplot(data = as.data.frame(rowRanges(data)), aes(x = rowRanges(data)$diffmean ,
                                 y = -1 * log10(rowRanges(data)$p.value.adj),
                                 colour = rowRanges(data)$threshold),
                environment = .e) + geom_point() +
        labs(title = title) + ylab(ylab) + xlab(xlab) +
        geom_vline(aes(xintercept = -diffmean.cut), colour = "black",
                   linetype = "dashed") +
        geom_hline(aes(yintercept = -1 * log10(p.cut)),
                   colour = "black",
                   linetype = "dashed") +
        scale_color_manual(breaks = c("1", "2", "3"),
                           values = color,
                           labels = label,
                           name = legend)
    # saving box plot to analyse it
    ggsave(p, filename = filename, width = 10, height = 5, dpi = 600)
    data <- subset(data, rowRanges(data)$threshold  == "2" | rowRanges(data)$threshold == "3")
    return(data)
}

# Get latest Genome Reference Consortium Human Build And save
# it as Genomic Ranges
#' @importFrom biomaRt useMart getBM
get.GRCh.bioMart <- function(genome="hg19") {

    if(genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                           host="grch37.ensembl.org",
                           path="/biomart/martservice" ,
                           dataset="hsapiens_gene_ensembl")
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

#' @title Preparing methylation data and expression data
#'
#' @description
#'  Preparing methylation data and expression data
#'
#' @param met methylation data
#' @param expression expression data
#' @return Dataframe with methylation and expression data
#' @importFrom GenomicRanges GRanges distanceToNearest
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' sample <- "TCGA-06-0939-01A-01D-1228-05"
#' query <- TCGAQuery(tumor = "GBM",samples = sample, level = 3)
#' TCGADownload(query,path = "exampleData",samples = sample, quiet = TRUE)
#' met <- TCGAPrepare(query, dir="exampleData",summarizedExperiment = FALSE)
#' met <-  met[ met$Chromosome != "X" &  met$Chromosome != "Y", ]
#' met <- na.omit(met)
#' met$diffmean <- runif(nrow(met),-3,3)
#' FDR <- runif(nrow(met),0,1)
#' p.value.adj <- runif(nrow(met),0,1)
#' GeneSymbol <- met$Gene_Symbol
#' dm <-  runif(nrow(met),-3,3)
#' expression <- data.frame(GeneSymbol,FDR,p.value.adj,dm)
#' gene.met <- starbursAnalysis(met,expression)
starbursAnalysis <- function(met, expression) {
    #### fix methylation gene names before merging.  map gene ID to
    #### genomic coordinates
    gene.GR <- GRanges(seqnames = paste0("chr", gene.location$chromosome_name),
                       ranges = IRanges(start = gene.location$start_position,
                                        end = gene.location$end_position),
                       strand = gene.location$strand,
                       symbol = gene.location$external_gene_name,
                       EntrezID = gene.location$entrezgene)

    probe.info <- GRanges(seqnames = paste0("chr", met$Chromosome),
                          ranges = IRanges(start = met$Genomic_Coordinate,
                                           end = met$Genomic_Coordinate),
                          probeID = met$Composite.Element.REF)

    # closest gene to each 450k probe ##your data
    distance <- as.data.frame(distanceToNearest(probe.info, gene.GR))
    rownames(gene.location) <- NULL
    gene.order.by.distance <- gene.location[distance$subjectHits,]
    gene.order.by.distance$distance <- as.matrix(distance$distance)

    data.nearest.gene <- merge(met,
                               gene.order.by.distance[, c("external_gene_name",
                                                          "entrezgene",
                                                          "distance")],
                               by.x="Gene_Symbol", by.y="external_gene_name")

    volcano <- merge(met, expression,
                     by.x = "Gene_Symbol",
                     by.y = "GeneSymbol",
                     incomparables = NA)

    volcano$ID <- paste(volcano$Gene_Symbol,
                        volcano$probeID,
                        volcano$cgID, sep = ".")

    # Preparing gene expression
    volcano$geFDR <- log10(volcano$FDR)
    volcano$geFDR2 <- volcano$geFDR
    volcano[volcano$dm > 0, "geFDR2"] <- -1 * volcano[volcano$dm > 0, "geFDR"]

    # Preparing methylation
    volcano$meFDR <- log10(volcano$p.value.adj)
    volcano$meFDR2 <- volcano$meFDR
    volcano[volcano$diffmean > 0, "meFDR2"] <-
        -1 * volcano[volcano$diffmean > 0, "meFDR"]

    return(volcano)
}

#' @title Create starburst plot
#'
#' @description
#'   Create starburst plot
#'
#' @details
#'    Input: data with gene expression/methylation expression
#'    Output: starburst plot
#'
#' @param data data frame columns: diffmean, p.value.adj
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
#' @import ggplot2
#' @export
#' @return Save a starburst plot
#' @examples
#' sample <- "TCGA-06-0939-01A-01D-1228-05"
#' query <- TCGAQuery(tumor = "GBM",samples = sample, level = 3)
#' TCGADownload(query,path = "exampleData",samples = sample, quiet = TRUE)
#' met <- TCGAPrepare(query, dir="exampleData", summarizedExperiment = FALSE)
#' met <-  met[ met$Chromosome != "X" &  met$Chromosome != "Y", ]
#' met <- na.omit(met)
#' met$diffmean <- runif(nrow(met),-3,3)
#' FDR <- runif(nrow(met),0,1)
#' p.value.adj <- runif(nrow(met),0,1)
#' GeneSymbol <- met$Gene_Symbol
#' dm <-  runif(nrow(met),-3,3)
#' expression <- data.frame(GeneSymbol,FDR,p.value.adj,dm)
#' gene.met <- starbursAnalysis(met,expression)
#' starburstPlot(gene.met)
starburstPlot <- function(data,
                          filename = "volcano.pdf",
                          ylab = paste0("Gene Expression\nlog10 of the",
                                        "adjusted Significance (FDR)"),
                          xlab = paste0("DNA Methylation\nlog10 of the",
                                        " adjusted Significance (FDR)"),
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
                                    "2" = "Up & Hypo",
                                    "3" = "Down & Hypo",
                                    "4" = "hypo",
                                    "5" = "hyper",
                                    "6" = "Up",
                                    "7" = "Down",
                                    "8" = "Up & Hyper",
                                    "9" = "Down & Hyper"),
                          xlim = NULL, ylim = NULL, p.cut = 0.05
)
{
    .e <- environment()
    volcano.m <- data
    volcano.m$threshold.starburst <- "1"
    volcano.m$threshold.size <- "1"

    # subseting by regulation (geFDR) and methylation level
    # (meFDR) down regulated up regulated lowerthr
    # |||||||||||||||| upperthr hypomethylated hipermethylated
    lowerthr <- log10(p.cut)  # - 1.30103
    upperthr <- (-lowerthr)  # +1.30103

    # Group 2:up regulated and hypomethylated
    a <- subset(volcano.m,
                volcano.m$geFDR2 > upperthr &
                    volcano.m$meFDR2 < lowerthr)

    # Group 3: down regulated and hypomethylated
    b <- subset(volcano.m,
                volcano.m$geFDR2 < lowerthr &
                    volcano.m$meFDR2 < lowerthr)

    # Group 4: hypomethylated
    c <- subset(volcano.m,
                volcano.m$geFDR2 > lowerthr &
                    volcano.m$geFDR2 < upperthr &
                    volcano.m$meFDR2 < lowerthr)

    # Group 5: hypermethylated
    d <- subset(volcano.m,
                volcano.m$geFDR2 > lowerthr &
                    volcano.m$geFDR2 < upperthr &
                    volcano.m$meFDR2 > upperthr)

    # Group 6: upregulated
    e <- subset(volcano.m,
                volcano.m$geFDR2 > upperthr &
                    volcano.m$meFDR2 < upperthr &
                    volcano.m$meFDR2 > lowerthr)

    # Group 7: downregulated
    f <- subset(volcano.m,
                volcano.m$geFDR2 < lowerthr &
                    volcano.m$meFDR2 < upperthr &
                    volcano.m$meFDR2 > lowerthr)

    # Group 8: upregulated and hypermethylated
    g <- subset(volcano.m,
                volcano.m$geFDR2 > upperthr &
                    volcano.m$meFDR2 > upperthr)

    # Group 9: downregulated and hypermethylated
    h <- subset(volcano.m,
                volcano.m$geFDR2 < lowerthr &
                    volcano.m$meFDR2 > upperthr)

    size <- c("1", "1", "1", "1", "1", "1", "1","1")
    groups <- c("2", "3", "4", "5", "6", "7","8","9")
    s <- list(a, b, c, d, e, f,g,h)
    for (i in seq_along(s)) {
        idx <- rownames(s[[i]])
        if (length(idx) > 0) {
            volcano.m[idx, "threshold.starburst"] <- groups[i]
            volcano.m[idx, "threshold.size"] <- size[i]
        }
    }

    ## starburst plot
    p <- ggplot(data = volcano.m, environment = .e,
                aes(x = volcano.m$meFDR2,
                    y = volcano.m$geFDR2,
                    colour = volcano.m$threshold.starburst,
                    size = volcano.m$threshold.size)) +
        geom_point()
    if (!is.null(xlim)) {
        p <- p + xlim(xlim)
    }
    if (!is.null(ylim)) {
        p <- p + ylim(ylim)
    }
    p <- p + labs(title) + ylab(ylab) + xlab(xlab)
    p <- p + scale_color_manual(values = color, labels = label,
                                name = legend)
    p + geom_hline(aes(yintercept = lowerthr)) +
        geom_hline(aes(yintercept = upperthr)) +
        geom_vline(aes(xintercept = lowerthr)) +
        geom_vline(aes(xintercept = upperthr))

    ggsave(filename = "starbust.gcimp.pdf", width = 14, height = 10)

    # return methylation < 0, expressao >0
}
