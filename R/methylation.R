#' @title Calculate diffmean methylation between two groups
#' @description
#'    Calculate diffmean methylation of probes between two groups removing lines
#'    that has NA values.
#' @param data SummarizedExperiment object obtained from TCGAPrepare
#' @param groupCol Columns in colData(data) that defines the groups.
#' @param group1 Name of group1 to be used in the analysis
#' @param group2 Name of group2  to be used in the analysis
#' @param save Save histogram of diffmean
#' @import ggplot2
#' @import graphics
#' @importFrom grDevices png dev.off
#' @importFrom S4Vectors values
#' @importFrom SummarizedExperiment colData rowRanges assay rowRanges<- values<-
#' @return Saves in the rowRages(data) the columns: mean.group1, mean.group2
#'        diffmean.group1.group2; Where group1 and group2 are the names of the
#'        groups.
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
#'  diff.mean <- TCGAbiolinks:::diffmean(data,groupCol = "group")
#' @keywords internal
diffmean <- function(data, groupCol = NULL, group1 = NULL, group2 = NULL, save = FALSE) {

    if (is.null(groupCol)) {
        message("Please, set the groupCol parameter")
        return(NULL)
    }
    if ( length(unique(colData(data)[,groupCol])) != 2 &&
         is.null(group1) && is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    } else if (length(unique(colData(data)[,groupCol])) == 2 &&
               is.null(group1) && is.null(group2)) {
        group1 <- unique(colData(data)[,groupCol])[1]
        group2 <- unique(colData(data)[,groupCol])[2]
    }
    message("Calculating the diference between the mean methylation of the groups...")

    m <- assay(data)
    idx1 <- which(colData(data)[,groupCol] == group1)
    idx2 <- which(colData(data)[,groupCol] == group2)
    mean.g1 <- rowMeans(m[,idx1], na.rm = TRUE)
    mean.g2 <- rowMeans(m[,idx2], na.rm = TRUE)
    diffmean <- mean.g2 - mean.g1

    # Saves the result
    group1.col <- gsub("[[:punct:]]| ", ".", group1)
    group2.col <- gsub("[[:punct:]]| ", ".", group2)
    values(rowRanges(data))[,paste0("mean.", group1.col)] <-  mean.g1
    values(rowRanges(data))[,paste0("mean.", group2.col)] <-  mean.g2
    values(rowRanges(data))[,paste0("diffmean.",group1.col,".", group2.col)] <-  diffmean
    values(rowRanges(data))[,paste0("diffmean.",group2.col,".", group1.col)] <-  -diffmean
    # Ploting a histogram to evaluate the data
    tryCatch({
        if(save) {
            fhist <- paste0("histogram_diffmean.",group1.col,group2.col,".png")
            message("Saving histogram of diffmean values: ", fhist)
            p <- qplot(diffmean,
                       geom="histogram",
                       binwidth = 0.5,
                       main = "Histogram for diffmeans",
                       xlab = "Diffmean",
                       fill=I("blue"))
            png(filename = fhist)
            print(p)
            dev.off()
        }
    })
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
#' @param labels labels of the plot
#' @param ylab y axis text of the plot
#' @param xlab x axis text of the plot
#' @param filename The name of the pdf file.
#' @param color Define the colors of the lines.
#' @param width Image width
#' @param height Image height
#' @param print.value Print pvalue in the plot? Default: TRUE
#' @param legend.position Legend position ("top", "right","left","bottom")
#' @param legend.title.position  Legend title position ("top", "right","left","bottom")
#' @param legend.ncols Number of columns of the legend
#' @param add.legend If true, legend is created. Otherwise names will
#' be added to the last point in the lines.
#' @param add.points If true, shows each death at the line of survival curves
#' @param dpi Figure quality
#' @importFrom GGally ggsurv
#' @importFrom survival survfit Surv
#' @importFrom scales percent
#' @importFrom ggthemes theme_base
#' @importFrom ggrepel geom_text_repel
#' @export
#' @return Survival plot
#' @examples
#' clin <- GDCquery_clinic("TCGA-LGG", type = "clinical", save.csv = FALSE)
#' TCGAanalyze_survival(clin, clusterCol="gender")
TCGAanalyze_survival <- function(data,
                                 clusterCol = NULL,
                                 legend = "Legend",
                                 labels = NULL,
                                 cutoff = 0,
                                 main = "Kaplan-Meier Overall Survival Curves",
                                 ylab = "Probability of survival",
                                 xlab = "Time since diagnosis (days)",
                                 filename = "survival.pdf",
                                 color = NULL,
                                 height = 8,
                                 width = 12,
                                 dpi = 300,
                                 legend.position = "inside",
                                 legend.title.position = "top",
                                 legend.ncols = 1,
                                 add.legend = TRUE,
                                 print.value = TRUE,
                                 add.points = TRUE
) {
    .e <- environment()

    if(!all(c("vital_status", "days_to_death","days_to_last_follow_up") %in% colnames(data)))
        stop("Columns vital_status, days_to_death and  days_to_last_follow_up should be in data frame")

    if(is.null(color)){
        color <- rainbow(length(unique(data[,clusterCol])))
    }

    group <- NULL
    if (is.null(clusterCol)) {
        stop("Please provide the clusterCol argument")
    } else if(length(unique(data[,clusterCol])) == 1) {
        stop( paste0("Sorry, but I'm expecting at least two groups\n",
                     "  Only this group found: ", unique(data[,clusterCol])))
    }
    notDead <- is.na(data$days_to_death)

    if (any(notDead == TRUE)) {
        data[notDead,"days_to_death"] <- data[notDead,"days_to_last_follow_up"]
    }
    # create a column to be used with survival package, info need
    # to be TRUE(DEAD)/FALSE (ALIVE)
    data$s <- grepl("dead",data$vital_status,ignore.case = TRUE)

    # Column with groups
    data$type <- as.factor(data[,clusterCol])

    # create the formula for survival analysis
    f.m <- formula(Surv(as.numeric(data$days_to_death),event=data$s) ~ data$type)
    fit <- survfit(f.m, data = data)

    # calculating p-value
    pvalue <- summary(coxph(
        Surv(as.numeric(data$days_to_death),event=data$s)
        ~ data$type))$logtest[3]


    surv <- ggsurv(fit, CI = "def", plot.cens = FALSE,
                   surv.col = "gg.def",
                   cens.col = "red", lty.est = 1,
                   lty.ci = 2, cens.shape = 3,
                   back.white = TRUE,
                   xlab = xlab, ylab = ylab, main = main)

    if (print.value){
        surv <- surv + annotate("text",x = -Inf,y = -Inf, hjust = -0.1,
                                vjust = -1.0, size = 6,
                                label = paste0("Log-Rank P-value = ",
                                               format(pvalue,
                                                      scientific = TRUE,
                                                      digits = 2)))
    }

    if (cutoff != 0) {
        surv <- surv + ggplot2::coord_cartesian(xlim = c(0, cutoff))
    }

    label.add.n <- function(x) {
        paste0(x, " (n = ",
               nrow(data[data[,clusterCol] == x,]), ")")
    }

    if(is.null(labels)){
        labels <- sapply(levels(data$type),label.add.n)
    }
    surv <- surv + scale_colour_manual(name = legend,
                                       labels = labels,
                                       values=color)
    if(add.points){
        surv <- surv + geom_point(aes(colour = group),
                                  shape = 3,size = 2)
    }
    surv <- surv + guides(linetype = FALSE) +
        scale_y_continuous(labels = scales::percent) +
        theme_base()

    if(add.legend == TRUE){
        if(legend.position == "inside"){
            surv <- surv +  theme(legend.justification=c(1,1),
                                  legend.background = element_rect(colour = "black"),
                                  legend.position=c(1,1))
        } else {
            surv <- surv +  theme(legend.position=legend.position)
        }
        surv <- surv +
            guides(color=guide_legend(override.aes=list(size=3)),
                   fill=guide_legend(ncol=legend.ncols,title.position = legend.title.position, title.hjust =0.5))

    }

    if(add.legend == FALSE){
        surv <- surv +  geom_text_repel(data=ddply(surv$data, .(group), function(x) x[nrow(x), ]),
                                        aes(label = group, color = factor(group)),
                                        segment.color = '#555555', segment.size = 0.0,
                                        size = 3, show.legend = FALSE) +
            theme(legend.position="none")
    }

    if(!is.null(filename)) {
        ggsave(surv, filename = filename, width = width, height = height, dpi = dpi)
    } else {
        return(surv)
    }
}
#' @title Mean methylation boxplot
#' @description
#'   Creates a mean methylation boxplot for groups (groupCol),
#'   subgroups will be highlited  as shapes if the subgroupCol was set.
#'
#'   Observation: Data is a summarizedExperiment.
#'
#' @param data SummarizedExperiment object obtained from TCGAPrepare
#' @param groupCol Columns in colData(data) that defines the groups. If no
#' columns defined a columns called "Patients" will be used
#' @param subgroupCol Columns in colData(data) that defines the subgroups.
#' @param shapes Shape vector of the subgroups. It must have the size of the levels
#' of the subgroups. Example: shapes = c(21,23) if for two levels
#' @param filename The name of the pdf that will be saved
#' @param subgroup.legend Name of the subgroup legend. DEFAULT: subgroupCol
#' @param group.legend Name of the group legend. DEFAULT: groupCol
#' @param color vector of colors to be used in graph
#' @param title main title in the plot
#' @param ylab y axis text in the plot
#' @param print.pvalue Print p-value for two groups
#' @param xlab x axis text in the plot
#' @param labels Labels of the groups
#' @param sort Sort boxplot by mean or median.
#' Possible values: mean.asc, mean.desc, median.asc, meadian.desc
#' @param plot.jitter Plot jitter? Default TRUE
#' @param jitter.size Plot jitter size? Default 3
#' @param height Plot height default:10
#' @param width Plot width default:10
#' @param dpi Pdf dpi default:600
#' @param order Order of the boxplots
#' @param axis.text.x.angle Angle of text in the x axis
#' @param y.limits Change lower/upper y-axis limit
#' @param legend.position Legend position ("top", "right","left","bottom")
#' @param legend.title.position  Legend title position ("top", "right","left","bottom")
#' @param legend.ncols Number of columns of the legend
#' @param add.axis.x.text Add text to x-axis? Default: FALSE
#' @import ggplot2 stats
#' @importFrom SummarizedExperiment colData rowRanges assay
#' @importFrom grDevices rainbow
# ' @importFrom gtools combinations
#' @importFrom plyr ddply . summarize
#' @importFrom knitr kable
#' @export
#' @return Save the pdf survival plot
#' @examples
#' nrows <- 200; ncols <- 21
#' counts <- matrix(runif(nrows * ncols, 0, 1), nrows)
#' rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                    IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                     strand=sample(c("+", "-"), 200, TRUE),
#'                     feature_id=sprintf("ID%03d", 1:200))
#'colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input","Other"), 7),
#'                     row.names=LETTERS[1:21],
#'                     group=rep(c("group1","group2","group3"),c(7,7,7)),
#'                     subgroup=rep(c("subgroup1","subgroup2","subgroup3"),7))
#'data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=counts),
#'          rowRanges=rowRanges,
#'          colData=colData)
#' TCGAvisualize_meanMethylation(data,groupCol  = "group")
#' # change lower/upper y-axis limit
#' TCGAvisualize_meanMethylation(data,groupCol  = "group", y.limits = c(0,1))
#' # change lower y-axis limit
#' TCGAvisualize_meanMethylation(data,groupCol  = "group", y.limits = 0)
#' TCGAvisualize_meanMethylation(data,groupCol  = "group", subgroupCol="subgroup")
#' TCGAvisualize_meanMethylation(data,groupCol  = "group")
#' TCGAvisualize_meanMethylation(data,groupCol  = "group",sort="mean.desc",filename="meandesc.pdf")
#' TCGAvisualize_meanMethylation(data,groupCol  = "group",sort="mean.asc",filename="meanasc.pdf")
#' TCGAvisualize_meanMethylation(data,groupCol  = "group",sort="median.asc",filename="medianasc.pdf")
#' TCGAvisualize_meanMethylation(data,groupCol  = "group",sort="median.desc",filename="mediandesc.pdf")
#' if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
TCGAvisualize_meanMethylation <- function(data,
                                          groupCol=NULL,
                                          subgroupCol=NULL,
                                          shapes = NULL,
                                          print.pvalue=FALSE,
                                          plot.jitter = TRUE,
                                          jitter.size=3,
                                          filename = "groupMeanMet.pdf",
                                          ylab = expression(
                                              paste("Mean DNA methylation (",
                                                    beta,"-values)")),
                                          xlab = NULL,
                                          title = "Mean DNA methylation",
                                          labels = NULL,
                                          group.legend = NULL,
                                          subgroup.legend = NULL,
                                          color = NULL,
                                          y.limits = NULL,
                                          sort,
                                          order,
                                          legend.position = "top",
                                          legend.title.position = "top",
                                          legend.ncols = 3,
                                          add.axis.x.text = FALSE,
                                          width=10,
                                          height=10,
                                          dpi=600,
                                          axis.text.x.angle = 90
) {
    .e <- environment()
    mean <- colMeans(assay(data),na.rm = TRUE)

    if (is.null(groupCol)){
        groups <- rep("Patient",length(mean))
    } else {
        if(!(groupCol %in% colnames(colData(data)))) stop("groupCol not found in the object")
        groups <- colData(data)[,groupCol]
    }

    if (is.null(subgroupCol)){
        subgroups <- NULL
    } else {
        if(!(subgroupCol %in% colnames(colData(data)))) stop("subgroupCol not found in the object")
        subgroups <- colData(data)[,subgroupCol]
    }

    if (!is.null(subgroupCol)){
        df <- data.frame(mean = mean, groups = groups, subgroups = subgroups)
    } else {
        df <- data.frame(mean = mean, groups = groups)
    }
    message("==================== DATA Summary ====================")
    data.summary <- ddply(df, .(groups), summarize,
                          Mean=mean(mean), Median=median(mean),
                          Max = max(mean),Min=min(mean))
    print(kable(data.summary))
    message("==================== END DATA Summary ====================")

    #comb2by2 <- combinations(length(levels(droplevels(df$groups))),
    #                  2,
    #                 levels(droplevels(df$groups)))
    groups <- levels(droplevels(df$groups))
    mat.pvalue <- matrix(ncol=length(groups),nrow=length(groups),
                         dimnames=list(groups,groups))
    if(length(groups) > 1){
        comb2by2 <- t(combn(levels(droplevels(df$groups)),2))

        for (i in 1:nrow(comb2by2)){
            try({
                aux <- t.test(mean ~ groups,
                              data = subset(df,subset=df$groups %in% comb2by2[i,]) )$p.value;
                mat.pvalue[comb2by2[i,1],comb2by2[i,2]] <- aux
                mat.pvalue[comb2by2[i,2],comb2by2[i,1]] <- aux
            },
            silent = TRUE
            )
        }
        message("==================== T test results ====================")
        print(kable(mat.pvalue))
        message("==================== END T test results ====================")

    }
    if(print.pvalue & length(levels(droplevels(df$groups))) == 2) {
        pvalue <- t.test(mean ~ groups, data = df)$p.value
    } else {
        print.pvalue <- FALSE
    }
    # Plot for methylation analysis Axis x: LGm clusters Axis y:
    # mean methylation
    label.add.n <- function(x) {
        paste0(x, " (n = ",
               nrow(subset(df,subset = (df$groups == x))), ")")
    }
    if (is.null(group.legend)) {
        group.legend <- groupCol
    }
    if (is.null(subgroup.legend)) {
        subgroup.legend <- subgroupCol
    }

    if(missing(sort)){
        if(missing(order)){
            x <- factor(df$groups)
        } else {
            x <- factor(df$groups,levels = order)
        }
    } else if(sort == "mean.asc") {
        x <- reorder(df$groups, df$mean, FUN="mean")
    } else  if(sort == "mean.desc") {
        x <- reorder(df$groups, -df$mean, FUN="mean")
    } else if(sort == "median.asc") {
        x <- reorder(df$groups, df$mean, FUN="median")
    } else if(sort == "median.desc") {
        x <- reorder(df$groups, -df$mean, FUN="median")
    }

    if (is.null(labels)) {
        labels <- levels(x)
        labels <-  sapply(labels,label.add.n)
    }

    if(is.null(color)){
        color <- rainbow(length(labels))
        color <- color[(match(levels(x),levels(factor(df$groups))))]
    }

    p <- ggplot(df, aes(x, df$mean),
                environment = .e) +
        geom_boxplot(aes(fill = x),
                     notchwidth = 0.25, outlier.shape = NA)

    if (plot.jitter){

        if (!is.null(subgroupCol)){

            p <- p + geom_jitter(aes(shape = subgroups,
                                     size =  subgroups),
                                 position = position_jitter(width = 0.1),
                                 size = jitter.size)
        } else {
            p <- p + geom_jitter(position = position_jitter(width = 0.1),
                                 size = jitter.size)
        }
    }
    if(add.axis.x.text){
        axis.text.x <- element_text(angle = axis.text.x.angle,
                                    vjust = 0.5,
                                    size = 16)
    } else {
        axis.text.x <-  element_blank()
    }
    p <- p + scale_fill_manual(values = color,labels = labels, name = group.legend)
    p <- p + scale_x_discrete(limits=levels(x))
    p <- p + ylab(ylab) + xlab(xlab) + labs(title = title) +
        labs(shape=subgroup.legend, color=group.legend) +
        theme_bw() +
        theme(axis.title.x = element_text(face = "bold", size = 20),
              axis.text.x = axis.text.x,
              axis.title.y = element_text(face = "bold",
                                          size = 20),
              axis.text.y = element_text(size = 16),
              plot.title = element_text(face = "bold", size = 16),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14),
              axis.text= element_text(size = 22),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line.x=element_line(colour = "black"),
              axis.line.y=element_line(colour = "black"),
              legend.position=legend.position,
              legend.key = element_rect(colour = 'white')) +
        guides(fill=guide_legend(ncol=legend.ncols,title.position = legend.title.position, title.hjust =0.5))

    if (!is.null(shapes)){
        p <- p + scale_shape_manual(values = shapes)
    }

    if (print.pvalue){
        p <- p + annotate("text",x = -Inf,y = -Inf, hjust = -0.1,
                          vjust = -1.0, size = 5,
                          label = paste0("P-value = ",
                                         format(pvalue,scientific = TRUE,
                                                digits = 2)))
    }
    if(!is.null(y.limits)){
        p <- p + expand_limits(x = 0, y = y.limits )
    }

    # saving box plot to analyse it
    if(!is.null(filename)){
        ggsave(p, filename = filename, width = width, height = height, dpi = dpi)
        message(paste("Plot saved in: ", file.path(getwd(),filename)))
    } else {
        return(p)
    }
}

#' @title Calculate pvalues
#' @description Calculate pvalues using wilcoxon test
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
#' @param cores Number of cores to be used
#' @param save Save histogram of pvalues
#' @return Data frame with cols p values/p values adjusted
#' @import graphics
#' @importFrom grDevices png dev.off pdf
#' @import stats
#' @importFrom parallel detectCores
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
#' @importFrom plyr adply
#' @importFrom stats wilcox.test
#' @importFrom doParallel registerDoParallel
#' @keywords internal
calculate.pvalues <- function(data,
                              groupCol = NULL,
                              group1 = NULL,
                              group2 = NULL,
                              paired = FALSE,
                              method = "BH",
                              exact = TRUE,
                              cores = 1, save = FALSE) {

    parallel <- FALSE
    if (cores > 1){
        if(is.windows()){
            if (cores > detectCores()) cores <- detectCores()
            registerDoParallel(cores)
            parallel = TRUE
        } else {
            if (cores > detectCores()) cores <- detectCores()
            registerDoParallel(cores)
            parallel = TRUE
        }
    }


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

    if(!is.factor(colData(data)[,groupCol])) {
        colData(data)[,groupCol] <- factor(
            colData(data)[,groupCol]
        )
    }

    p.value <- adply(assay(data),1,
                     function(x) {
                         aux <- data_frame(beta=x[c(idx1,idx2)],
                                           cluster=droplevels(
                                               colData(data)[c(idx1,idx2),
                                                             groupCol]))
                         wilcox.test(beta ~ cluster,
                                     data=aux, #exact = TRUE,
                                     paired = paired)$p.value
                     }, .progress = "text", .parallel = parallel
    )
    p.value <- p.value[,2]
    ## Plot a histogram
    if(save) {
        message("Saved histogram_pvalues.png...")
        png(filename = "histogram_pvalues.png")
        hist(p.value)
        dev.off()
    }

    ## Calculate the adjusted p-values by using Benjamini-Hochberg
    ## (BH) method
    p.value.adj <- p.adjust(p.value, method = method)

    ## Plot a histogram
    if(save) {
        message("Saved histogram_pvalues_adj.png")
        png(filename = "histogram_pvalues_adj.png")
        hist(p.value.adj)
        dev.off()
    }
    #Saving the values into the object
    group1.col <- gsub("[[:punct:]]| ", ".", group1)
    group2.col <- gsub("[[:punct:]]| ", ".", group2)
    colp <- paste("p.value",  group1.col,  group2.col, sep = ".")
    values(rowRanges(data))[,colp] <-  p.value
    coladj <- paste("p.value.adj", group1.col,  group2.col, sep = ".")
    values(rowRanges(data))[,coladj] <-  p.value.adj
    colp <- paste("p.value",  group2.col,  group1.col, sep = ".")
    values(rowRanges(data))[,colp] <-  p.value
    coladj <- paste("p.value.adj", group2.col,  group1.col, sep = ".")
    values(rowRanges(data))[,coladj] <-  p.value.adj

    return(data)
}

#' @title Creates a volcano plot for DNA methylation or expression
#' @description Creates a volcano plot from the
#' expression and methylation analysis.
#' @details
#'    Creates a volcano plot from the expression and methylation analysis.
#'    Please see the vignette for more information
#'    Observation: This function automatically is called by TCGAanalyse_DMR
#' @param x x-axis data
#' @param y y-axis data
#' @param y.cut p-values threshold. Default: 0.01
#' @param x.cut  x-axis threshold. Default: 0.0
#' @param filename Filename. Default: volcano.pdf, volcano.svg, volcano.png
#' @param legend Legend title
#' @param color vector of colors to be used in graph
#' @param title main title. If not specified it will be
#' "Volcano plot (group1 vs group2)
#' @param ylab y axis text
#' @param xlab x axis text
#' @param xlim x limits to cut image
#' @param ylim y limits to cut image
#' @param height Figure height
#' @param width Figure width
#' @param names Names to be ploted if significant.
#' Should be the same size of x and y
#' @param names.fill Names should be filled in a color box?  Default: TRUE
#' @param names.size Size of the names text
#' @param dpi Figure dpi
#' @param label vector of labels to be used in the figure.
#' Example: c("Not Significant","Hypermethylated in group1",
#' "Hypomethylated in group1"))#'
#' @param highlight List of genes/probes to be highlighted. It should be in the names argument.
#' @param highlight.color Color of the points highlighted
#' @param show.names What names will be showd? Possibilities: "both", "significant", "highlighted"
#' @export
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @return Saves the volcano plot in the current folder
#' @examples
#' x <- runif(200, -1, 1)
#' y <- runif(200, 0.01, 1)
#' TCGAVisualize_volcano(x,y)
#' TCGAVisualize_volcano(x,y,filename = NULL,y.cut = 10000000,x.cut=0.8,
#'                       names = rep("AAAA",length(x)), legend = "Status",
#'                       names.fill = FALSE)
#' TCGAVisualize_volcano(x,y,filename = NULL,y.cut = 10000000,x.cut=0.8,
#'                       names = as.character(1:length(x)), legend = "Status",
#'                       names.fill = TRUE, highlight = c("1","2"),show="both")
#' while (!(is.null(dev.list()["RStudioGD"]))){dev.off()}
TCGAVisualize_volcano <- function(x,y,
                                  filename = "volcano.pdf",
                                  ylab =  expression(paste(-Log[10],
                                                           " (FDR corrected -P values)")),
                                  xlab=NULL, title=NULL, legend=NULL,
                                  label=NULL, xlim=NULL, ylim=NULL,
                                  color = c("black", "red", "green"),
                                  names=NULL,
                                  names.fill= TRUE,
                                  show.names="significant",
                                  x.cut=0,
                                  y.cut=0.01,
                                  height=5,
                                  width=10,
                                  highlight=NULL,
                                  highlight.color = "orange",
                                  names.size = 4,
                                  dpi = 300){

    if(!is.null(names)) {
        if(all(grepl("\\|",names))){
            names <- strsplit(names,"\\|")
            names <- unlist(lapply(names,function(x) x[1]))
        }
    }
    .e <- environment()
    threshold <- rep("1",length(x))
    names(color) <- as.character(1:3)

    if(is.null(label)) {
        label = c("1" = "Not Significant",
                  "2" = "Up regulated",
                  "3" = "Down regulated")
    } else  {
        names(label) <- as.character(1:3)
    }
    # get significant data
    sig <-  y < y.cut
    sig[is.na(sig)] <- FALSE
    # hypermethylated/up regulated samples compared to old state
    up <- x  > x.cut
    up[is.na(up)] <- FALSE
    if (any(up & sig)) threshold[up & sig] <- "2"

    # hypomethylated/ down regulated samples compared to old state
    down <-  x < (-x.cut)
    down[is.na(down)] <- FALSE
    if (any(down & sig)) threshold[down & sig] <- "3"

    if(!is.null(highlight)){
        idx <- which(names %in% highlight)
        if(length(idx) >0 ){
            print(idx)
            threshold[which(names %in% highlight)]  <- "4"
            color <- c(color,highlight.color)
            names(color) <- as.character(1:4)
        }
    }
    df <- data.frame(x=x,y=y,threshold=threshold)

    # As last color should be the highlighthed, we need to order all the vectors
    if(!is.null(highlight)){
        order.idx <-  order(df$threshold)
        down <- down[order.idx]
        sig <- sig[order.idx]
        up <- up[order.idx]
        df <- df[order.idx, ]
        names <- names[order.idx]
    }
    # Plot a volcano plot
    p <- ggplot(data=df,
                aes(x = x , y = -1 * log10(y), colour = threshold ),
                environment = .e) +
        geom_point() +
        ggtitle(title) + ylab(ylab) + xlab(xlab) +
        geom_vline(aes(xintercept = -x.cut),
                   colour = "black",linetype = "dashed") +
        geom_vline(aes(xintercept = x.cut),
                   colour = "black", linetype = "dashed") +
        geom_hline(aes(yintercept = -1 * log10(y.cut)),
                   colour = "black", linetype = "dashed") +
        scale_color_manual(breaks = as.numeric(names(label)),
                           values = color,
                           labels = label,
                           name = legend) +
        theme_bw() + theme(panel.border = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.text = element_text(size = 10),
                           axis.line.x=element_line(colour = "black"),
                           axis.line.y=element_line(colour = "black"),
                           legend.position="top",
                           legend.key = element_rect(colour = 'white'))

    # Label points with the textxy function from the calibrate plot
    if(!is.null(names)){
        # With the names the user can highlight the significant genes, up and down
        # or the ones highlighted
        if(show.names == "significant"){
            idx <- (up & sig) | (down & sig)
            important <- c("2","3")
        } else if(show.names == "highlighted") {
            if(!is.null(highlight)){
                idx <- (names %in% highlight)
                important <- c("4")
            } else {
                message("Missing highlight argument")
                return(NULL)
            }
        } else if(show.names == "both"){
            if(!is.null(highlight)){
                idx <- (up & sig) | (down & sig) |  (names %in% highlight)
                important <- c("2","3","4")
            } else {
                message("Missing highlight argument")
                return(NULL)
            }
        } else {
            message("Wrong highlight argument")
            return(NULL)
        }

        if(any(threshold %in% important)){
            if(names.fill){
                p <- p + geom_label_repel(
                    data = subset(df, threshold %in% important),
                    aes(label = names[idx],fill=threshold),
                    size = names.size, show.legend = FALSE,
                    fontface = 'bold', color = 'white',
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")
                ) +   scale_fill_manual(values=color[as.numeric(important)])
            }  else {
                p <- p + geom_text_repel(
                    data = subset(df, threshold %in% important),
                    aes(label = names[idx]),
                    size = names.size, show.legend = FALSE,
                    fontface = 'bold', color = 'black',
                    point.padding = unit(0.3, "lines"),
                    box.padding = unit(0.5, 'lines')
                )
            }
        }
    }

    if(!is.null(filename)){
        ggsave(p, filename = filename, width = width, height = height, dpi = dpi)
    } else {
        return(p)
    }
}


#' @title Differentially methylated regions Analysis
#' @description
#'   This function will search for differentially methylated CpG sites,
#'   which are regarded as possible functional regions involved
#'   in gene transcriptional regulation.
#'
#'   In order to find these regions we use the beta-values (methylation values
#'   ranging from 0.0 to 1.0) to compare two groups.
#'
#'   Firstly, it calculates the difference between the mean methylation of each
#'   group for each probes. Secondly, it calculates the p-value using the
#'   wilcoxon test using the Benjamini-Hochberg adjustment method.
#'   The default parameters will require a minimum absolute beta values delta
#'   of 0.2 and a false discovery rate (FDR)-adjusted Wilcoxon rank-sum P-value
#'   of < 0.01 for the difference.
#'
#'   After these analysis, we save a volcano plot (x-axis:diff mean methylation,
#'   y-axis: significance) that will help the user identify the differentially
#'   methylated CpG sites and return the object with the calculus in the rowRanges.
#'
#'   If the calculus already exists in the object it will not recalculated.
#'   You should set overwrite parameter to TRUE to force it, or remove the
#'   collumns with the results from the object.
#'
#' @param data  SummarizedExperiment obtained from the TCGAPrepare
#' @param groupCol  Columns with the groups inside the SummarizedExperiment
#'  object. (This will be obtained by the function colData(data))
#' @param group1 In case our object has more than 2 groups, you should set
#' the name of the group
#' @param group2 In case our object has more than 2 groups, you should set
#' the name of the group
#' @param calculate.pvalues.probes In order to get the probes faster the user can select to calculate the pvalues
#' only for the probes with a difference in DNA methylation. The default is to calculate to all probes.
#' Possible values: "all", "differential". Default "all"
#' @param plot.filename Filename. Default: volcano.pdf, volcano.svg, volcano.png. If set to FALSE, there will be no plot.
#' @param legend Legend title
#' @param color vector of colors to be used in graph
#' @param title main title. If not specified it will be
#' "Volcano plot (group1 vs group2)
#' @param ylab y axis text
#' @param xlab x axis text
#' @param xlim x limits to cut image
#' @param ylim y limits to cut image
#' @param label vector of labels to be used in the figure.
#' Example: c("Not Significant","Hypermethylated in group1",
#' "Hypomethylated in group1"))
#' @param p.cut p values threshold. Default: 0.01
#' @param probe.names is probe.names
#' @param diffmean.cut diffmean threshold. Default: 0.2
#' @param adj.method Adjusted method for the p-value calculation
#' @param paired Wilcoxon paired parameter. Default: FALSE
#' @param overwrite Overwrite the pvalues and diffmean values if already in the object
#' for both groups? Default: FALSE
#' @param save Save object with results? Default: TRUE
#' @param save.directory Directory to save the files. Default: working directory
#' @param filename Name of the file to save the object.
#' @param cores Number of cores to be used in the non-parametric test
#' Default = groupCol.group1.group2.rda
#' @import ggplot2
#' @importFrom SummarizedExperiment colData rowRanges assay rowRanges<- values<- SummarizedExperiment metadata<-
#' @importFrom S4Vectors metadata
#' @importFrom dplyr data_frame
#' @importFrom methods as
#' @import readr
#' @import utils
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
#' SummarizedExperiment::colData(data)$group <- c(rep("group 1",ncol(data)/2),
#'                          rep("group 2",ncol(data)/2))
#' hypo.hyper <- TCGAanalyze_DMR(data, p.cut = 0.85,"group","group 1","group 2")
#' SummarizedExperiment::colData(data)$group2 <- c(rep("group_1",ncol(data)/2),
#'                          rep("group_2",ncol(data)/2))
#' hypo.hyper <- TCGAanalyze_DMR(data, p.cut = 0.85,"group2","group_1","group_2")
TCGAanalyze_DMR <- function(data,
                            groupCol=NULL,
                            group1=NULL,
                            group2=NULL,
                            calculate.pvalues.probes = "all",
                            plot.filename = "methylation_volcano.pdf",
                            ylab =  expression(paste(-Log[10],
                                                     " (FDR corrected -P values)")),
                            xlab =  expression(paste(
                                "DNA Methylation difference (",beta,"-values)")
                            ),
                            title = NULL,
                            legend = "Legend",
                            color = c("black",  "red", "darkgreen"),
                            label = NULL,
                            xlim = NULL,
                            ylim = NULL,
                            p.cut = 0.01,
                            probe.names = FALSE,
                            diffmean.cut = 0.2,
                            paired = FALSE,
                            adj.method="BH",
                            overwrite=FALSE,
                            cores = 1,
                            save=TRUE,
                            save.directory = ".",
                            filename=NULL) {
    .e <- environment()

    names(color) <- as.character(1:3)
    # Check if object is a summarized Experiment
    if(class(data)!= class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        stop(paste0("Sorry, but I'm expecting a Summarized Experiment object, but I got a: ", class(data)))
    }
    # Check if object has NAs
    if(any(rowSums(!is.na(assay(data))))== 0){
        stop(paste0("Sorry, but we found some probes with NA for all samples in your data, please either remove/or replace them"))
    }

    if (is.null(groupCol)) {
        message("Please, set the groupCol parameter")
        return(NULL)
    }
    if(!(groupCol %in% colnames(colData(data)))){
        stop(paste0("column ",groupCol, " not found in the object"))
    }

    if ( length(unique(colData(data)[,groupCol])) != 2 &&
         is.null(group1) && is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    } else if (length(unique(colData(data)[,groupCol])) == 2  && (
        is.null(group1) || is.null(group2)) ) {
        group1 <- unique(colData(data)[,groupCol])[1]
        group2 <- unique(colData(data)[,groupCol])[2]
    } else {
        message(paste0("Group1:", group1))
        message(paste0("Group2:", group2))
    }

    # Check if groups has at least one sample
    if(!any(colData(data)[,groupCol] == group1,na.rm = TRUE)){
        stop(paste0("Sorry, but ", group1, " has no samples" ))
    }
    if(!any(colData(data)[,groupCol] == group2,na.rm = TRUE)){
        stop(paste0("Sorry, but ", group2, " has no samples" ))
    }


    # defining title and label if not specified by the user
    if (is.null(title)) {
        title <- paste("Volcano plot", "(", group2, "vs", group1,")")
    }

    if (is.null(label)) {
        label <- c("Not Significant",
                   "Hypermethylated",
                   "Hypomethylated")
        label[2:3] <-  paste(label[2:3], "in", group2)
    }
    group1.col <- gsub("[[:punct:]]| ", ".", group1)
    group2.col <- gsub("[[:punct:]]| ", ".", group2)
    diffcol <- paste("diffmean", group1.col, group2.col,sep = ".")
    if (!(diffcol %in% colnames(values(data))) || overwrite) {
        data <- diffmean(data,groupCol, group1 = group1, group2 = group2, save = save)
        if (!(diffcol %in% colnames(values(rowRanges(data))))) {
            stop(paste0("Error! Not found ", diffcol))
        }
    }

    pcol <- paste("p.value.adj", group2.col, group1.col,sep = ".")
    if(!(pcol %in% colnames(values(data)))){
        pcol <- paste("p.value.adj", group1.col, group2.col, sep = ".")
    }
    if (!(pcol %in% colnames(values(data))) | overwrite) {
        if(calculate.pvalues.probes == "all"){
            data <- calculate.pvalues(data, groupCol, group1, group2,
                                      paired = paired,
                                      method = adj.method,
                                      cores = cores,
                                      save = save)
        } else  if(calculate.pvalues.probes == "differential"){
            message(paste0("Caculating p-values only for probes with a difference of mean methylation equal or higher than ", diffmean.cut))
            print(diffcol)
            print(colnames(values(data)))
            diff.probes <- abs(values(data)[,diffcol]) > diffmean.cut
            nb <- length(which(diff.probes == TRUE))
            if(nb == 0) {
                warning("No probes differenly methylated")
                return(NULL)
            }
            print(paste0("Number of probes differenly methylated: ",nb))
            data <- calculate.pvalues(data[diff.probes,], groupCol, group1, group2,
                                      paired = paired,
                                      method = adj.method,
                                      cores = cores,
                                      save = save)
        }

        # An error should not happen, if it happens (probably due to an incorret
        # user input) we will stop
        if (!(pcol %in% colnames(values(data))))  stop(paste0("Error! Not found ", pcol))
    }
    log <- paste0("TCGAanalyze_DMR.",gsub(" ", ".",group1),".",gsub(" ", ".",group2))
    assign(log,c("groupCol" = groupCol,
                 "group1" = group1.col,
                 "group2" = group2.col,
                 "plot.filename" = plot.filename,
                 "xlim" = xlim,
                 "ylim" = ylim,
                 "p.cut" = p.cut,
                 "diffmean.cut" = diffmean.cut,
                 "paired" = "paired",
                 "adj.method" = adj.method))
    metadata(data)[[log]] <- (eval(as.symbol(log)))
    statuscol <- paste("status",group1.col,group2.col,sep = ".")
    statuscol2 <- paste("status",group2.col,group1.col,sep = ".")
    values(data)[,statuscol] <-  "Not Significant"
    values(data)[,statuscol2] <-  "Not Significant"

    # get significant data
    sig <-  values(data)[,pcol] < p.cut
    sig[is.na(sig)] <- FALSE
    # hypermethylated samples compared to old state
    hyper <- values(data)[,diffcol]  > diffmean.cut
    hyper[is.na(hyper)] <- FALSE
    if (any(hyper & sig)) values(data)[hyper & sig,statuscol] <- "Hypermethylated"
    if (any(hyper & sig)) values(data)[hyper & sig,statuscol2] <- "Hypomethylated"

    # hypomethylated samples compared to old state
    hypo <-  values(data)[,diffcol] < (-diffmean.cut)
    hypo[is.na(hypo)] <- FALSE

    if (any(hypo & sig)) values(data)[hypo & sig,statuscol] <- "Hypomethylated"
    if (any(hypo & sig)) values(data)[hypo & sig,statuscol2] <- "Hypermethylated"

    # Plot a volcano plot
    names <- NULL
    if(probe.names) names <- values(data)$probeID

    if(plot.filename != FALSE) {
        TCGAVisualize_volcano(x = values(data)[,diffcol],
                              y = values(data)[,pcol],
                              filename = plot.filename,
                              ylab =  ylab,
                              xlab = xlab,
                              title = title,
                              legend= legend,
                              label = label,
                              names = names,
                              x.cut = diffmean.cut,
                              y.cut = p.cut)
    }
    if (save) {

        # saving results into a csv file
        csv <- paste0(paste("DMR_results",
                            gsub("_",".",groupCol),
                            group1.col,
                            group2.col,
                            "pcut",p.cut,"meancut",diffmean.cut,  sep = "_"),".csv")
        dir.create(save.directory,showWarnings = FALSE,recursive = TRUE)
        csv <- file.path(save.directory,csv)
        message(paste0("Saving the results also in a csv file: "), csv)
        df <- values(data)
        if (any(hyper & sig)) df[hyper & sig,statuscol] <- paste("Hypermethylated","in", group2)
        if (any(hyper & sig)) df[hyper & sig,statuscol2] <- paste("Hypomethylated","in", group1)
        if (any(hypo & sig)) df[hypo & sig,statuscol] <- paste("Hypomethylated","in", group2)
        if (any(hypo & sig)) df[hypo & sig,statuscol2] <- paste("Hypermethylated","in", group1)
        # get metadata not created by this function
        idx <- grep("mean|status|value",colnames(df),invert = TRUE)

        write_csv(as.data.frame(df[,
                                   c(colnames(df)[idx],
                                     paste("mean", group1.col,sep = "."),
                                     paste("mean", group2.col,sep = "."),
                                     paste("diffmean", group1.col, group2.col, sep = "."),
                                     paste("p.value", group1.col, group2.col, sep = "."),
                                     paste("p.value.adj", group1.col, group2.col, sep = "."),
                                     statuscol,
                                     paste("diffmean",group2.col,group1.col,sep = "."),
                                     paste("p.value",group2.col,group1.col,sep = "."),
                                     paste("p.value.adj",group2.col,group1.col,sep = "."),
                                     statuscol2)
                                   ]),path =  csv)
        if (is.null(filename)) {
            filename <- paste0(paste(
                gsub("_",".",groupCol),
                group1.col,
                group2.col,
                "pcut",p.cut,
                "meancut",diffmean.cut,
                sep = "_"),
                ".rda")
            filename <- file.path(save.directory, filename)
        }

        # saving results into R object
        save(data, file = filename)
    }
    return(data)
}

#' @title Create starburst plot
#'
#' @description
#'   Create Starburst plot for comparison of DNA methylation and gene expression.
#'    The log10 (FDR-corrected P value) is plotted for beta value for DNA
#'    methylation (x axis) and gene expression (y axis) for each gene.
#'
#'    The black dashed line shows the FDR-adjusted P value of 0.01.
#'
#'    You can set names to TRUE to get the names of the significant genes.
#'
#'    Candidate biologically significant genes will be circled in the plot.
#'
#'    Candidate biologically significant are the genes that respect the
#'    expression (logFC.cut), DNA methylation (diffmean.cut) and
#'    significance thresholds (exp.p.cut, met.p.cut)
#'
#' @details
#'    Input: data with gene expression/methylation expression
#'    Output: starburst plot
#'
#' @param met A SummarizedExperiment with methylation data obtained from the
#' TCGAPrepare or Data frame from DMR_results file. Expected colData columns: diffmean,  p.value.adj  and p.value
#' Execute volcanoPlot function in order to obtain these values for the object.
#' @param exp Object obtained by DEArnaSEQ function
#' @param filename The filename of the file (it can be pdf, svg, png, etc)
#' @param legend legend title
#' @param color vector of colors to be used in graph
#' @param label vector of labels to be used in graph
#' @param title main title
#' @param names Add the names of the significant genes? Default: FALSE
#' @param names.fill Names should be filled in a color box?  Default: TRUE
#' @param circle Circle pair gene/probe that respects diffmean.cut and logFC.cut
#' Default: TRUE
#' @param ylab y axis text
#' @param xlab x axis text
#' @param xlim x limits to cut image
#' @param ylim y limits to cut image
#' @param met.p.cut methylation p value cut-off
#' @param exp.p.cut expression p value cut-off
#' @param diffmean.cut If set, the probes with diffmean higher
#' than methylation cut-off will be
#'  highlighted in the plot. And the data frame return will be subseted.
#' @param logFC.cut If set, the probes with expression fold
#' change higher than methylation cut-off will be
#'  highlighted in the plot. And the data frame return will be subseted.
#' @param group1 The name of the group 1
#' Obs: Column p.value.adj.group1.group2 should exist
#' @param group2 The name of the group 2.
#' Obs: Column p.value.adj.group1.group2 should exist
#' @param return.plot If true only plot object will be returned (pdf will not be created)
#' @param height Figure height
#' @param width Figure width
#' @param dpi Figure dpi
#' @import ggplot2
#' @importFrom SummarizedExperiment rowRanges rowRanges<- values<-
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @export
#' @return Save a starburst plot
#' @examples
#' library(SummarizedExperiment)
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
#'                   logFC=runif(20000, -5, 5),
#'                   FDR=runif(20000, 0.01, 1))
#' rowRanges(met)$diffmean.g1.g2 <- c(runif(20000, -0.1, 0.1))
#' rowRanges(met)$diffmean.g2.g1 <- -1*(rowRanges(met)$diffmean.g1.g2)
#' rowRanges(met)$p.value.g1.g2 <- c(runif(20000, 0, 1))
#' rowRanges(met)$p.value.adj.g1.g2 <- c(runif(20000, 0, 1))
#' result <- TCGAvisualize_starburst(met,exp,
#'                                   exp.p.cut = 0.05, met.p.cut = 0.05,
#'                                   group1="g1",group2="g2",
#'                                   diffmean.cut=0.0,
#'                                   names=TRUE, circle = FALSE)
#' result <- TCGAvisualize_starburst(SummarizedExperiment::values(met),
#'                                   exp,
#'                                   exp.p.cut = 0.05, met.p.cut = 0.05,
#'                                   group1="g1",group2="g2",
#'                                   diffmean.cut=0.0,
#'                                   names=TRUE, circle = FALSE)
TCGAvisualize_starburst <- function(met,
                                    exp,
                                    group1=NULL,
                                    group2=NULL,
                                    exp.p.cut = 0.01,
                                    met.p.cut = 0.01,
                                    diffmean.cut = 0,
                                    logFC.cut = 0,
                                    names = FALSE,
                                    names.fill = TRUE,
                                    circle = TRUE,
                                    filename = "starburst.pdf",
                                    return.plot = FALSE,
                                    ylab = expression(atop("Gene Expression",
                                                           paste(Log[10],
                                                                 " (FDR corrected P values)"))),
                                    xlab = expression(atop("DNA Methylation",
                                                           paste(Log[10],
                                                                 " (FDR corrected P values)"))),
                                    title = "Starburst Plot",
                                    legend = "DNA Methylation/Expression Relation",
                                    color = NULL,
                                    label = c("Not Significant",
                                              "Up regulated & Hypo methylated",
                                              "Down regulated & Hypo methylated",
                                              "hypo methylated",
                                              "hyper methylated",
                                              "Up regulated",
                                              "Down regulated",
                                              "Up regulated & Hyper methylated",
                                              "Down regulated & Hyper methylated"),
                                    xlim = NULL, ylim = NULL,
                                    height=10,
                                    width=20,
                                    dpi=600
)
{
    .e <- environment()

    group1.col <- gsub("[[:punct:]]| ", ".",group1)
    group2.col <- gsub("[[:punct:]]| ", ".",group2)

    if(is.null(color)) color <- c("#000000", "#E69F00","#56B4E9", "#009E73",
                                  "red", "#0072B2","#D55E00", "#CC79A7",
                                  "purple")
    names(color) <- as.character(1:9)
    names(label) <- as.character(1:9)
    names.color <- color
    names(names.color) <- label

    if ( is.null(group1) || is.null(group2)) {
        message("Please, set the group1 and group2 parameters")
        return(NULL)
    }

    if (class(met) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
        met <- values(met)
    }

    # Preparing methylation
    pcol <- gsub("[[:punct:]]| ", ".",paste("p.value.adj",group1,group2,sep = "."))
    if(!(pcol %in%  colnames(met))){
        pcol <- gsub("[[:punct:]]| ", ".",paste("p.value.adj",group2,group1,sep = "."))
    }
    if(!(pcol %in%  colnames(met))){
        stop("Error! p-values adjusted not found. Please, run TCGAanalyze_DMR")
    }

    # Methylation matrix and expression matrix should have the same name column for merge
    idx <- grep("Gene_symbol",colnames(met),ignore.case = TRUE)
    colnames(met)[idx] <- "Gene_symbol"

    # Check if gene symbol columns exists
    if(!any(grepl("Gene_symbol",colnames(exp),ignore.case = FALSE))) {
        if("mRNA" %in% colnames(exp)) {
            if(all(grepl("\\|",exp$mRNA))) {
                exp$Gene_symbol <- unlist(lapply(strsplit(exp$mRNA,"\\|"),function(x) x[2]))
            } else {
                exp$Gene_symbol <- exp$mRNA
            }
        } else {
            aux <- rownames(exp)
            if(all(grepl("\\|",aux))) {
                exp$Gene_symbol <- unlist(lapply(strsplit(aux,"\\|"),function(x) x[2]))
            } else {
                exp$Gene_symbol <- aux
            }
        }
    } else {
        # Check if it has the same pattern
        idx <- grep("Gene_symbol",colnames(exp),ignore.case = TRUE)
        colnames(exp)[idx] <- "Gene_symbol"
    }


    volcano <- merge(met, exp, by = "Gene_symbol")
    volcano$ID <- paste(volcano$Gene_symbol,
                        volcano$probeID, sep = ".")

    # Preparing gene expression
    volcano$geFDR <- log10(volcano$FDR)
    volcano$geFDR2 <- volcano$geFDR
    volcano[volcano$logFC > 0, "geFDR2"] <-
        -1 * volcano[volcano$logFC > 0, "geFDR"]

    diffcol <- gsub("[[:punct:]]| ", ".",paste("diffmean",group1,group2,sep = "."))
    volcano$meFDR <- log10(volcano[,pcol])
    volcano$meFDR2 <- volcano$meFDR
    idx <- volcano[,diffcol] > 0
    idx[is.na(idx)] <- FALSE # handling NAs
    volcano[idx, "meFDR2"] <- -1 * volcano[idx, "meFDR"]

    label[2:9] <-  paste(label[2:9], "in", group2)

    # subseting by regulation (geFDR) and methylation level
    # (meFDR) down regulated up regulated lowerthr
    # |||||||||||||||| upperthr hypomethylated hipermethylated
    met.lowerthr <- log10(met.p.cut)
    met.upperthr <- (-met.lowerthr)

    exp.lowerthr <- log10(exp.p.cut)
    exp.upperthr <- (-exp.lowerthr)

    # Group 2:up regulated and hypomethylated
    a <- subset(volcano,
                volcano$geFDR2 > exp.upperthr &
                    volcano$meFDR2 < met.lowerthr)

    a.sig <- subset(a, abs(a[,diffcol]) > diffmean.cut &
                        abs(a$logFC) > logFC.cut)

    # Group 3: down regulated and hypomethylated
    b <- subset(volcano,
                volcano$geFDR2 < exp.lowerthr &
                    volcano$meFDR2 < met.lowerthr)

    b.sig <- subset(b, abs(b[,diffcol]) > diffmean.cut &
                        abs(b$logFC) > logFC.cut)

    # Group 4: hypomethylated
    c <- subset(volcano,
                volcano$geFDR2 > exp.lowerthr &
                    volcano$geFDR2 < exp.upperthr &
                    volcano$meFDR2 < met.lowerthr)

    # Group 5: hypermethylated
    d <- subset(volcano,
                volcano$geFDR2 > exp.lowerthr &
                    volcano$geFDR2 < exp.upperthr &
                    volcano$meFDR2 > met.upperthr)

    # Group 6: upregulated
    e <- subset(volcano,
                volcano$geFDR2 > exp.upperthr &
                    volcano$meFDR2 < met.upperthr &
                    volcano$meFDR2 > met.lowerthr)

    # Group 7: downregulated
    f <- subset(volcano,
                volcano$geFDR2 < exp.lowerthr &
                    volcano$meFDR2 < met.upperthr &
                    volcano$meFDR2 > met.lowerthr)

    # Group 8: upregulated and hypermethylated
    g <- subset(volcano,
                volcano$geFDR2 > exp.upperthr &
                    volcano$meFDR2 > met.upperthr)

    g.sig <- subset(g, abs(g[,diffcol]) > diffmean.cut &
                        abs(g$logFC) > logFC.cut)

    # Group 9: downregulated and hypermethylated
    h <- subset(volcano,
                volcano$geFDR2 < exp.lowerthr &
                    volcano$meFDR2 > met.upperthr)

    h.sig <- subset(h, abs(h[,diffcol]) > diffmean.cut &
                        abs(h$logFC) > logFC.cut)

    groups <- as.character(seq(2,9))

    # return methylation < 0, expressao >0
    volcano$starburst.status  <-  "Not Significant"
    volcano$shape <- "1"
    volcano$threshold.starburst <- "1"
    volcano$threshold.size <- "1"

    state <- c("Up regulated & Hypo methylated",
               "Down regulated & Hypo methylated",
               "hypo methylated",
               "hyper methylated",
               "Up regulated",
               "Down regulated",
               "Up regulated & Hyper methylated",
               "Down regulated & Hyper methylated")

    s <- list(a, b, c, d, e, f, g, h)
    for (i in seq_along(s)) {
        idx <- rownames(s[[i]])
        if (length(idx) > 0) {
            volcano[idx, "threshold.starburst"] <- groups[i]
            volcano[idx, "starburst.status"] <-  state[i]
        }
    }

    size <- rep(2,4)
    shape <-  as.character(rep(2,4))
    volcano_normal <- volcano
    significant <- NULL
    s <- list(a.sig,b.sig,g.sig,h.sig)
    for (i in seq_along(s)) {
        idx <- rownames(s[[i]])
        if (length(idx) > 0) {
            volcano[idx, "threshold.size"] <- size[i]
            volcano[idx, "shape"] <-  shape[i]
            significant <-  rbind(significant,volcano[idx,])
        }
    }

    ## starburst plot
    p <- ggplot(data = volcano_normal, environment = .e,
                aes(x = volcano_normal$meFDR2,
                    y = volcano_normal$geFDR2,
                    colour = volcano_normal$threshold.starburst)) +
        geom_point()
    #p <- p + scale_shape_discrete(
    #    labels = c("Candidate Biologically Significant"),
    #    name = "Biological importance")

    if(!is.null(significant) & circle){
        p <- p + geom_point( data = significant,
                             aes(x = significant$meFDR2,
                                 y = significant$geFDR2),
                             color = "black",
                             shape=1,
                             size = 8,
                             show.legend = FALSE)
    }

    if(names == TRUE & !is.null(significant)){
        message("Adding names to genes")
        if(names.fill){
            p <- p + geom_label_repel(
                data = significant,
                aes(x = significant$meFDR2, y =  significant$geFDR2,
                    label = significant$Gene_symbol, fill = as.factor(significant$starburst.status)),
                size = 4, show.legend = FALSE,
                fontface = 'bold', color = 'white',
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            ) + scale_fill_manual(values=names.color)
        }  else {
            p <- p + geom_text_repel(
                data = significant,
                aes(x = significant$meFDR2, y =  significant$geFDR2,
                    label = significant$Gene_symbol, fill = significant$starburst.status),
                size = 4, show.legend = FALSE,
                fontface = 'bold', color = 'black',
                point.padding = unit(0.3, "lines"),
                box.padding = unit(0.5, 'lines')
            )
        }
    }


    if (!is.null(xlim)) {
        p <- p + xlim(xlim)
    }
    if (!is.null(ylim)) {
        p <- p + ylim(ylim)
    }
    p <- p + ggtitle(title) + ylab(ylab) + xlab(xlab) + guides(size=FALSE)
    p <- p + scale_color_manual(values = color, labels = label, name = legend) +
        guides(col = guide_legend(nrow = 3))

    p <-  p + geom_hline(aes(yintercept = exp.lowerthr), colour = "black",
                         linetype = "dashed") +
        geom_hline(aes(yintercept = exp.upperthr), colour = "black",
                   linetype = "dashed") +
        geom_vline(aes(xintercept = met.lowerthr), colour = "black",
                   linetype = "dashed") +
        geom_vline(aes(xintercept = met.upperthr), colour = "black",
                   linetype = "dashed")  +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line.x=element_line(colour = "black"),
              axis.line.y=element_line(colour = "black"),
              legend.position="top",
              legend.key = element_rect(colour = 'white'),
              plot.title = element_text(face = "bold", size = 16),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14),
              axis.text= element_text(size = 14),
              axis.title.x = element_text(face = "bold", size = 14),
              axis.text.x = element_text(vjust = 0.5,
                                         size = 14),
              axis.title.y = element_text(face = "bold",
                                          size = 14),
              axis.text.y = element_text(size = 14))
    #p <- p + geom_point( data = significant,
    #                     aes(x = meFDR2,
    #                         y = geFDR2), shape=1, size = 10,show_guide = FALSE)
    #    p <- p + scale_shape_discrete(
    #        labels = c("Candidate Biologically Significant"),
    #        name = "Biological importance")
    if(!return.plot) ggsave(filename = filename, width = width, height = height, dpi = dpi)

    #statuscol <- paste("status", group1, group2, sep = ".")

    #volcano <- subset(volcano,select = c("Gene_Symbol",
    #                                     "probeID",statuscol,
    #                                     "starburst.status")
    #)
    volcano$shape <- NULL
    volcano$threshold.starburst <- NULL
    volcano$threshold.size <- NULL

    volcano <- subset(volcano, volcano$geFDR <= exp.lowerthr &
                          volcano$meFDR <= met.lowerthr)

    if (diffmean.cut != 0) {
        volcano <- subset(volcano, abs(volcano[,diffcol]) > diffmean.cut)
    }
    if (logFC.cut != 0){
        volcano <- subset(volcano, abs(volcano$logFC) >= logFC.cut)
    }

    if(return.plot) {
        return(list(plot=p,starburst=volcano))
    }

    return(volcano)
}
