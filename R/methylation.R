#' @title calculate diffmean methylation between two groups
#' @description
#'    calculate diffmean methylation between two groups
#' @param group1 Data frame probes vs patient
#' @param group2 Data frame probes vs patient
#' @import ggplot2
#' @export
#' @return dataframe with diffmean values
diffmean <- function(group1, group2) {
  g1 <- group1
  g1$mean.g1 <- apply(group1, 1, mean, na.rm = TRUE)
  g1$probeID <- rownames(group1)

  g2 <- group2
  g2$mean.g2 <- apply(group2, 1, mean, na.rm = TRUE)
  g2$probeID <- rownames(group2)

  g1.g2 <- merge(g1[, c("probeID", "mean.g1")],
                 g2[, c("probeID","mean.g2")],
                 by = "probeID")
  rownames(g1.g2) <- g1.g2[, 1]
  g1.g2$diffmean <- g1.g2$mean.g1 - g1.g2$mean.g2

  png(filename = "histogram_diffmean.png")
  hist(g1.g2$diffmean)
  dev.off()

  return(g1.g2)
}

#' @title creates survival analysis
#' @description Creates survival analysis
#'
#' @param met.md Data frame with the following columns: vital, os, cluster
#' @param legend Legend title of the figure
#' @param cutoff xlim
#' @param main main title
#' @param ylab y axis text
#' @param xlab x axis text
#' @param time.reference In order not to forget -> MONTHS DAY YEAR
#' @param filename output file name
#' @param default.plot You can use ggplot or plot
#' @param color not implemented yet
#' @importFrom GGally ggsurv
#' @importFrom survival survfit Surv
#' @export
#' @return Survival plot
survivalPlot <- function(met.md, legend = "Legend", cutoff = 0,
                         main = "Kaplan-Meier Overall Survival Curves",
                         ylab = "PROBABILITY OF SURVIVAL",
                         xlab = "TIME SINCE DIAGNOSIS",
                         time.reference = "DAYS",
                         filename = "survival.pdf",
                         color = c("green", "firebrick4", "orange3", "blue"),
                         default.plot = "ggplot") {
  with(met.md,{
    if (cutoff != 0) {
      # Axis-x cut-off
      aux <- subset(met.md, met.md$death_days_to > cutoff)
      # 'cut' info bigger than cutoff
      met.md[rownames(aux), "death_days_to"] <- cutoff
      # Pacient is alive (0:alive,1:dead)
      met.md[rownames(aux), "vital_status"] <- "Alive"
    }
    # create a column to be used with survival package, info need
    # to be TRUE(DEAD)/FALSE (ALIVE)
    met.md$s <- met.md$vital == "Dead"

    # Column with groups
    met.md$type <- as.factor(met.md$cluster)

    # create the formula for survival analysis
    f.m <- formula(Surv(as.numeric(met.md$death_days_to), s) ~ type)
    fit <- survfit(f.m, data = met.md)

    pdf(file = filename)
    if (default.plot == "plot") {
      plot(fit, lwd = 4, col = color, main = main,
           xlab = paste0(xlab,"(", time.reference, ")"),
           ylab = ylab, yscale = 100,
           bg = "black")
      box(col = "black", lwd = 3)
      par(xpd = TRUE)
      legend("right", legend = sapply(seq_along(fit$n), function(x) {
        paste0("group", x, " (n=", fit$n[x], ")")
      }), col = color, lwd = 3, title = legend, box.lwd = 3,
      bg = "white")
      dev.off()
    } else {
      ## Using ggplot
      surv <- ggsurv(fit, CI = "def", plot.cens = FALSE,
                     surv.col = "gg.def",
                     cens.col = "red", lty.est = 1,
                     lty.ci = 2, cens.shape = 3,
                     back.white = TRUE,
                     xlab = paste0(xlab, "(", time.reference,
                                   ")"), ylab = ylab, main = main)
      if (cutoff != 0) {
        surv <- surv + ggplot2::coord_cartesian(xlim = c(0, cutoff))
      }
      label.add.n <- function(x) {
        paste0(x, "(n=", nrow(met.md[met.md$cluster == x,]), ")")
      }
      surv <- surv + scale_colour_discrete(
        name = "legend",
        labels = sapply(levels(met.md$type), label.add.n))
      surv <- surv + geom_point(aes(colour = met.md$group), shape = 3,
                                size = 2)
      surv <- surv + guides(linetype = FALSE) +
        scale_y_continuous(formatter = 'percent')
      ggsave(surv, filename = filename)
    }
  })
}
#' @title organize TCGA Methylation MetaData
#'
#' @description
#'   Organize TCGA methylation metadata for the mean methylation analysis.
#'
#' @param wd Directory with the files
#' @examples
#' met.md <- organizeMethylationMetaDataFrame("data")
#' @export
#' @return \code{invisible (metadata)}
organizeMethylationMetaDataFrame <- function(wd = NULL) {
  oldwd <- getwd()
  setwd(wd)
  on.exit(setwd(oldwd))

  dirs <- list.dirs()
  files <- NULL
  for (i in seq_along(dirs)) {
    files <- c(files, paste0(dirs[i], "/", list.files(dirs[i])))
  }
  files <- files[grep(".*clinical_patient.*", files)]
  header <- read.table(files[1], nrows = 1, header = FALSE,
                       sep = "\t", stringsAsFactors = FALSE)
  metadata <- read.table(files[1], header = TRUE, sep = "\t",
                         skip = 2)
  colnames(metadata) <- unlist(header)
  invisible(metadata)
}

#' @title organize TCGA Methylation Data
#'
#' @description
#'   Organize TCGA methylation data for the mean methylation analysis.
#'
#' @details
#'    Input: a directory with DNA methylation data,
#'    Output: a data frame with DNA methylation values
#'    where rows are the probes names and columns are paciente ID
#'    Execution: read all files inside the directory and merge it by
#'    probes (Composite.Element.REF)
#' @examples
#' met <- organizeMethylationDataFrame(wd = "data")
#' @param wd Directory with the files
#' @return Methylation betavalues table
#' @export
organizeMethylationDataFrame <- function(wd = getwd()) {
  oldwd <- getwd()
  setwd(wd)
  on.exit(setwd(oldwd))

  dirs <- list.dirs()
  files <- NULL
  for (i in seq_along(dirs)) {
    files <- c(files, paste0(dirs[i], "/", list.files(dirs[i])))
  }
  files <- files[grep("TCGA.*txt", files)]

  for (i in seq_along(files)) {
    header <- readLines(files[i], n = 1)
    IDpatient <- substr(header, 19, 46)
    data <- read.table(files[i], header = TRUE, sep = "\t",
                       skip = 1)
    colnames(data)[2] <- IDpatient
    if (!exists("methylation")) {
      methylation <- data[, c(1, 3:5, 2)]
    } else {
      methylation <- merge(methylation, data[, c(1, 2)],
                           by = "Composite.Element.REF")
    }
  }

  # Use Composite.Element.REF as name, instead of column
  rownames(methylation) <- methylation$Composite.Element.REF
  # methylation$Composite.Element.REF <- NULL

  # remove X Y chromossomes
  message("Removing X Y chromossomes")
  methylation <- methylation[methylation$Chromosome != "X" &
                               methylation$Chromosome != "Y", ]
  # methylation$Chromosome <- NULL

  # remove NA lines
  message("Removing NA Lines")
  methylation <- na.omit(methylation)

  invisible(methylation)
}
#' @title Mean methylation boxplot
#' @description
#'   Mean methylation boxplot.
#'   Input: a dataframe with two columns:
#'    1 - values of mean methylation
#'    2 - Groups it belongs
#' @param data data frame first col mean methylation by patient, second groups
#' @param filename pdf filename
#' @param legend legend title
#' @param color vector of colors to be used in graph
#' @param title main title
#' @param ylab y axis text
#' @param xlab x axis text
#' @param sort Sort by mean methylation? False by default
#' @import ggplot2
#' @export
#' @return Save the survival plot
#'
#' @examples
#'# -----------------
#'# |      data       |
#'# -----------------
#'# | mean.met|cluster|
#'# -------------------
#'   avg <- runif(500,0,1)
#'   cluster <- c('Lgm1','Lgm2','Lgm3','Lgm4','Lgm5','Lgm6')
#'   cluster <- sample(cluster, 500,replace = TRUE)
#'   data <- data.frame(avg,cluster)
#'   metMeanBoxplot(data)
metMeanBoxplot <- function(data, sort = FALSE,
                           filename = "G-CIMP-mean.methylation.pdf",
                           ylab = "Mean DNA methylation",
                           xlab = "DNA Methylation Clusters",
                           title = "Mean DNA methylation by cluster",
                           legend = "Legend",
                           color = c("green", "red", "purple",
                                     "orange", "salmon", "grey")) {

  # Plot for methylation analysis Axis x: LGm clusters Axis y:
  # mean methylation
  if (sort) {
    p <- ggplot(data, aes(reorder(factor(data$cluster), data$avg),
                          data$avg))
    p <- p + geom_boxplot(aes(fill = reorder(factor(data$cluster),
                                             data$avg)),
                          notchwidth = 0.25) +
      geom_jitter(height = 0,position = position_jitter(width = 0.1),
                  size = 3) +
      scale_fill_manual(values = color,
                        labels = levels(reorder(factor(data$cluster),
                                                data$avg)),
                        name = legend)
  } else {
    p <- ggplot(data, aes(factor(data$cluster), data$avg))
    p <- p + geom_boxplot(aes(fill = factor(data$cluster)),
                          notchwidth = 0.25) +
      geom_jitter(height = 0,
                  position = position_jitter(width = 0.1),
                  size = 3) +
      scale_fill_manual(values = color,
                        labels = levels(factor(data$cluster)),
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

}

#' @title Calculate pvalues
#' @details
#'    Input: two matrix to be compared with wilcoxon test,
#'           a boolean value to do a paired or non-paired test
#'    Output: p-values (non-adj/adj) histograms, p-values (non-adj/adj)
#' @param values  Dataframe with values
#' @param idx1  Index of the values in group 1
#' @param idx2  Index of the values in group 2
#' @param paired  Do a paired wilcoxon test? Default: True
#' @param exact  Do a exact wilcoxon test? Default: True
#' @param cores  How many cores to be used Default: parallel::detectCores()
#' @return Data frame with cols p values/p values adjusted
#' @importFrom exactRankTests wilcox.exact
#' @importFrom parallel mclapply detectCores
#' @export
#' @return Data frame with two cols
#'         p-values/p-values adjusted
calculate.pvalues <- function(values, idx1, idx2, paired = TRUE,
                              exact = TRUE, cores = NULL) {
  if (is.null(cores)) {
    cores <- parallel::detectCores()
  }
  # Apply Wilcoxon test in order to calculate the p-values
  w.p.values <- unlist(mclapply(values, function(probe) {
    zz <- wilcox.exact(as.matrix(probe[idx1]), as.matrix(probe[idx2]),
                       exact = exact, paired = paired)
    z <- zz$p.value
    return(z)
  }, mc.cores = cores))

  ## Plot a histogram
  png(filename = "histogram_pvalues.png")
  hist(w.p.values)
  dev.off()

  ## Calculate the adjusted p-values by using Benjamini-Hochberg
  ## (BH) method
  w.p.values.adj <- p.adjust(w.p.values, method = "BH")

  png(filename = "histogram_pvalues_adj.png")

  ## Plot a histogram
  hist(w.p.values.adj)
  dev.off()

  return(data.frame(w.p.values, w.p.values.adj))
}

#' @title Volcano plot
#' @description
#'   Volcano plot - hypo/hyper methylation graphic.
#'   Input: a dataframe with columns:
#'    1 - p.value.adj - p value adjusted
#'    2 - diffmean - difference between old stated and new state of methylation
#' @param data data frame columns: diffmean, p.value.adj
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
#' @export
#' @return A dataframe with the Composite.Element.REF and
#'         the group it was classified
#'          group 1 = Not Significant
#'          group 2 = Hypermethylated
#'          group 3 = Hypomethylated
volcanoPlot <- function(data, filename = "volcano.pdf",
                        ylab = "- 1*log10 of the Significance",
                        xlab = "DNA Methylation",
                        title = "Volcano Plot TCGA GBM Tumors",
                        legend = "Legend",
                        color = c("1" = "black", "2" = "green",
                                  "3" = "red"),
                        label = c("1" = "Not Significant",
                                  "2" = "Hypermethylated",
                                  "3" = "Hypomethylated"),
                        xlim = NULL, ylim = NULL, p.cut = 0.05,
                        diffmean.cut = 0) {
  .e <- environment()
  data$threshold <- "1"

  # get significant data
  data.s <- subset(data, data$p.value.adj < p.cut)

  # hypermethylated samples compared to old state
  hyper <- subset(data.s, data.s$diffmean < (-diffmean.cut))
  data[rownames(hyper), "threshold"] <- "2"

  # hypomethylated samples compared to old state
  hypo <- subset(data.s, data$diffmean > diffmean.cut)
  data[rownames(hypo), "threshold"] <- "3"

  # Plot a volcano plot
  p <- ggplot(data = data, aes(x = diffmean,
                               y = -1 * log10(data$p.value.adj),
                               colour = data$threshold),
              environment = .e) +
    geom_point()
  if (!is.null(xlim)) {
    p <- p + xlim(xlim)
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }
  p <- p + labs(title = title) + ylab(ylab) + xlab(xlab)
  p <- p + geom_vline(aes(xintercept = -diffmean.cut), colour = "black",
                      linetype = "dashed") +
    geom_vline(aes(xintercept = diffmean.cut),
               colour = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = -1 * log10(p.cut)),
               colour = "black",
               linetype = "dashed")
  p <- p + scale_color_manual(breaks = c("1", "2", "3"),
                              values = color,
                              labels = label,
                              name = legend)
  # saving box plot to analyse it
  ggsave(p, filename = filename, width = 10, height = 5, dpi = 600)
  data <- subset(data, data$threshold == "2" | data$threshold == "3")
  return(data[, c("Composite.Element.REF", "threshold")])
}


# Get latest Genome Reference Consortium Human Build And save
# it as Genomic Ranges
#' @importFrom biomaRt useMart getBM
get.GRCh.bioMart <- function() {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  chrom <- c(1:22, "X", "Y")
  gene.location <- getBM(attributes = c("ensembl_gene_id",
                                        "ensembl_transcript_id",
                                        "chromosome_name",
                                        "start_position",
                                        "end_position", "strand",
                                        "external_gene_name",
                                        "external_transcript_name",
                                        "external_gene_source",
                                        "external_transcript_source_name",
                                        "hgnc_id", "entrezgene"),
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
#' # met Methylation data
#' # experssion Expression data
#' gene.met <- starbursAnalysis(met,expression)
#' starburstplot(gene.met)
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
  data.nearest.gene <- cbind(met,
                             gene.order.by.distance[, c("external_gene_name",
                                                        "entrezgene",
                                                        "distance")])

  # volcano$gene <- met$Gene_Symbol volcano$cgID <-
  # rownames(volcano)
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
    - 1 * volcano[volcano$diffmean > 0, "meFDR"]

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
#' @param diffmean.cut diffmean cut-off
#' @import ggplot2
#' @export
#' @return Save a starburst plot
#' @examples
#' # met Methylation data
#' # experssion Expression data
#' gene.met <- starbursAnalysis(met,expression)
#' starburstplot(gene.met)
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
                          xlim = NULL, ylim = NULL, p.cut = 0.05,
                          diffmean.cut = 0)
{
  .e <- environment()
  volcano.m <- data
  volcano.m$threshold.starburst <- "1"
  volcano.m$threshold.size <- "1"

  # subseting by regulation (geFDR) and methylation level
  # (meFDR) down regulated up regulated lowerthr
  # |||||||||||||||| upperthr hypomethylated hipermethylated
  lowerthr <- log10(0.05)  # - 1.30103
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
  f <- subset(volcano.m,
              volcano.m$geFDR2 > upperthr &
                volcano.m$meFDR2 > upperthr)

  # Group 9: downregulated and hypermethylated
  f <- subset(volcano.m,
              volcano.m$geFDR2 < lowerthr &
                volcano.m$meFDR2 > upperthr)

  size <- c("3", "3", "2", "2", "2", "2", "2","2")
  groups <- c("2", "3", "4", "5", "6", "7","8","9")
  s <- list(a, b, c, d, e, f)
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
