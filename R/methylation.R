#' @title calculate diffmean methylation between primary and recurent
#' @description
#'    calculate diffmean methylation between primary and recurent
#'    samples
#' @param data Data frame probes vs patient
#' @import ggplot2
#' @export
diffmean.prim.rec <- function(data){
  pat.prim <- "TCGA-[0-9a-z]{2}-[0-9a-z]{4}-01"
  pat.rec <- "TCGA-[0-9a-z]{2}-[0-9a-z]{4}-02"
  rec <- data[,grep(pat.rec,colnames(data))]
  prim <- data[,grep(pat.prim,colnames(data))]
  prim.rec <- diffmean(prim,rec)
  return(prim.rec)
}

diffmean <- function(group1,group2){
  g1 <- group1
  g1$mean.g1 <- apply(group1,1,mean,na.rm=TRUE)
  g1$probeID <- rownames(group1)

  g2 <- group2
  g2$mean.g2 <- apply(group2,1,mean,na.rm=TRUE)
  g2$probeID <- rownames(group2)

  g1.g2 <- merge(g1[,c("probeID","mean.g1")],g2[,c("probeID","mean.g2")],by="probeID")
  rownames(g1.g2) <- g1.g2[,1]
  g1.g2$diffmean <- g1.g2$mean.g1-g1.g2$mean.g2

  png(filename="histogram_diffmean.png")
  hist(g1.g2$diffmean)
  dev.off()

  return(g1.g2)
}

#' @title creates survival analysis
#' @description Creates survival analysis
#'
#' @param met.md Data frame with the following columns: vital, os, cluster
#'        legend Legend title of the figure
#'        cutoff xlim
#'        main main title
#'        ylab y axis text
#'        xlab x axis text
#'        time.reference In order not to forget -> MONTHS DAY YEAR
#'        filename output file name
#'        default.plot You can use ggplot or plot
#'        color not implemented yet
#' @import survival GGally
#' @export
survivalPlot <- function(met.md,
                         legend = "Legend",
                         cutoff=0,
                         main ="Kaplan-Meier Overall Survival Curves\nTCGA Samples\nLGG",
                         ylab = 'PROBABILITY OF SURVIVAL',
                         xlab = "TIME SINCE DIAGNOSIS",
                         time.reference = "DAYS",
                         filename = "survival.pdf",
                         color = c("green","firebrick4","orange3","blue"),
                         default.plot="ggplot"
){

  if(cutoff != 0){
    # Axis-x cut-off
    aux <- subset(met.md, death_days_to > cutoff)
    # 'cut' info bigger than cutoff
    met.md[rownames(aux),"death_days_to"] <- cutoff
    # Pacient is alive (0:alive,1:dead)
    met.md[rownames(aux), "vital_status"] <- "Alive"
  }
  # create a column to be used with survival package,
  # info need to be TRUE(DEAD)/FALSE (ALIVE)
  met.md$s <- met.md$vital == "Dead"

  # Column with groups
  met.md$type <- as.factor(met.md$cluster)

  # create the formula for survival analysis
  f.m = formula(Surv(as.numeric(death_days_to),s) ~ type)
  fit = survfit(f.m, data=met.md)

  pdf(file = filename)
  if(default.plot=="plot"){
    plot(fit,
         lwd=4,
         col=color,
         main=main,
         xlab=paste0(xlab,"(", time.reference,")"),
         ylab=ylab,
         yscale=100,
         bg="black"
    )
    box(col="black", lwd=3)
    par(xpd=TRUE)
    legend("right",
           legend=sapply(seq_along(fit$n),function(x){paste0("group",x," (n=", fit$n[x],")")}),
           col=color,
           lwd=3,
           title=legend,
           box.lwd=3,
           bg="white")
    dev.off()
  } else{
    ## Using ggplot
    surv <- ggsurv(fit,
                   CI         = 'def',
                   plot.cens  = FALSE,
                   surv.col   = 'gg.def',
                   cens.col   = 'red',
                   lty.est    = 1,
                   lty.ci     = 2,
                   cens.shape = 3,
                   back.white = TRUE,
                   xlab       = paste0(xlab,"(", time.reference,")"),
                   ylab       = ylab,
                   main       = main)
    if(cutoff != 0){
      surv <- surv + ggplot2::coord_cartesian(xlim = c(0, cutoff))
    }
    label.add.n <- function(x){paste0(x,"(n=",nrow(met.md[met.md$cluster==x,]),")")}
    surv <- surv + ggplot2::scale_colour_discrete(name = "legend",
                                                  labels = sapply(levels(met.md$type),label.add.n)
    )
    surv <- surv + geom_point(aes(colour=group),shape=3,size=2)
    surv <- surv + ggplot2::guides(linetype = FALSE) +
      scale_y_continuous(labels = percent)
    ggsave(surv, file=filename)
  }
}
#' @title organize TCGA Methylation MetaData
#'
#' @description
#'   Organize TCGA methylation metadata for the mean methylation analysis.
#'
#' @param dir Directory with the files
#' @import stringr
#' @export
#' @return \code{invisible (metadata)}
organizeMethylationMetaDataFrame <- function(wd = NULL){
  oldwd <- getwd()
  setwd(wd)
  on.exit(setwd(oldwd))

  dirs <- list.dirs()
  files <- NULL
  for(i in seq_along(dirs)){
    files <- c(files,paste0(dirs[i],"/",list.files(dirs[i])))
  }
  files <- files[grep(".*clinical_patient.*",files)]
  header <- read.table(files[1], nrows = 1, header = FALSE, sep ="\t", stringsAsFactors = FALSE)
  metadata <- read.table(files[1], header=T, sep="\t", skip=2)
  colnames( metadata ) <- unlist(header)
  invisible (metadata)
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
#'
#' @author
#' Thathi Malta  \email{tathimalta@@gmail.com}
#'
#' Maintainer: Tiago Chedraoui Silva \email{tiagochst@@gmail.com}
#'
#' @param wd Directory with the files
#' @return Methylation table
#' @import stringr
#' @export
organizeMethylationDataFrame <- function(wd = getwd()){
  oldwd <- getwd()
  setwd(wd)
  on.exit(setwd(oldwd))

  dirs <- list.dirs()
  files <- NULL
  for(i in seq_along(dirs)){
    files <- c(files,paste0(dirs[i],"/",list.files(dirs[i])))
  }
  files <- files[grep("TCGA.*txt",files)]

  for(i in seq_along(files)) {  #enquanto files[i] não for na, faz o que tá dizendo
    header <- readLines(files[i], n=1)
    IDpatient <- substr(header, 19, 46)
    data <- read.table(files[i], header=T, sep="\t", skip=1)
    colnames(data)[2] <- IDpatient
    if(!exists("methylation")) {
      methylation <- data[, c(1,3:5,2)]
    }
    else{
      methylation <- merge(methylation, data[,c(1,2)], by="Composite.Element.REF")
    }
  }

  # Use Composite.Element.REF as name, instead of column
  rownames(methylation) <- methylation$Composite.Element.REF
  #methylation$Composite.Element.REF <- NULL

  # remove X Y chromossomes
  message("Removing X Y chromossomes")
  methylation <- methylation[methylation$Chromosome !="X" & methylation$Chromosome!="Y",]
  #methylation$Chromosome <- NULL

  # remove NA lines
  message("Removing NA Lines")
  methylation <- na.omit(methylation)

  invisible (methylation)
}
#' @title Mean methylation boxplot
#' @description
#'   Mean methylation boxplot.
#'   Input: a dataframe with two columns:
#'    1 - values of mean methylation
#'    2 - Groups it belongs
#' @param data data frame first column mean methylation by patient, second groups
#' @param filename pdf filename
#' @param legend legend title
#' @param color vector of colors to be used in graph
#' @param title main title
#' @param ylab y axis text
#' @param xlab x axis text
#' @param sort Sort by mean methylation? False by default
#' @import ggplot2
#' @export
#'
#' @examples
#'# -----------------
#'# |      data       |
#'# -----------------
#'# | mean.met|cluster|
#'# -------------------
#' \dontrun{
#'   met.mean.boxplot(data)
#'}
#' \dontrun{
#'   mean <- runif(500,0,1)
#'   cluster <- c("Lgm1","Lgm2","Lgm3","Lgm4","Lgm5","Lgm6")
#'   cluster.vec <- sample(cluster, 500,replace = T)
#'   data <- data.frame(mean,cluster.vec)
#'   met.mean.boxplot(data)
#'}
met.mean.boxplot <- function (data,
                              sort = F,
                              filename="G-CIMP-mean.methylation.pdf",
                              ylab   = "Mean DNA methylation",
                              xlab   = "DNA Methylation Clusters",
                              title  = "Mean DNA methylation by cluster - G-CIMP samples",
                              legend ="Legend",
                              color = c("green","red","purple","orange","salmon","grey")
){

  # Plot for methylation analysis
  #   Axis x: LGm clusters
  #   Axis y:  mean methylation
  if(sort){
    p <- ggplot(data, aes(reorder(factor(cluster),avg), avg))
    p <- p +
      geom_boxplot(aes(fill = reorder(factor(cluster),avg)), notchwidth = 0.25) +
      geom_jitter(height = 0, position = position_jitter(width = .1), size=3) +
      scale_fill_manual(
        values = color,
        labels = levels(reorder(factor(box.plot$cluster),box.plot$avg)),
        name = legend
      )
  }
  else{
    p <- ggplot(data, aes(factor(cluster), avg))
    p <- p +
      geom_boxplot(aes(fill = factor(cluster)), notchwidth = 0.25) +
      geom_jitter(height = 0, position = position_jitter(width = .1), size=3) +
      scale_fill_manual(
        values = color,
        labels = levels(factor(data$cluster)),
        name = legend
      )
  }
  p <- p + ylab(ylab) + xlab(xlab) + labs(title = title) +
    theme(axis.title.x = element_text (face  = "bold", size = 20),
          axis.text.x  = element_text (angle = 90, vjust = 0.5, size = 16),
          axis.title.y = element_text (face  = "bold", size = 20),
          axis.text.y  = element_text (size  = 16),
          plot.title   = element_text (face  = "bold", size = 16)
    )

  # saving box plot to analyse it
  ggsave(p, filename=filename,
         width = 10, heigh = 10, dpi = 600)

}

#' @title Calculate pvalues
#' @details
#'    Input: two matrix to be compared with wilcoxon test,
#'           a boolean value to do a paired or non-paired test
#'    Output: p-values (non-adj/adj) histograms, p-values (non-adj/adj)
#' @param matrix1  Matrix to be compared with matrix 2
#' @param matrix2  Matrix to be compared with matrix 1
#' @param paired  Do a paired wilcoxon test? Default: True
#' @return Methylation table
#' @import exactRankTests parallel
#' @export
calculate.pvalues <- function (values,idx1,idx2,paired=TRUE,exact=TRUE,mc.cores=parallel::detectCores()){
  # Apply Wilcoxon test in order to calculate the p-values
  #print(idx1)
  #print(idx2)
  w.p.values <- unlist(mclapply(values,function(probe) {
    zz <- wilcox.exact(as.matrix(probe[idx1]),as.matrix(probe[idx2]), exact=exact,paired=paired)
    z <- zz$p.value
    return(z)
  }, mc.cores=mc.cores))
  #print(w.p.values)
  ##Plot a histogram
  png(filename="histogram_pvalues.png")
  hist(w.p.values)
  dev.off()
  ##Calculate the adjusted p-values by using the Benjamini-Hochberg (BH) method
  w.p.values.adj <- p.adjust(w.p.values,method="BH")

  png(filename="histogram_pvalues_adj.png")
  ##Plot a histogram
  hist(w.p.values.adj)
  dev.off()

  return(data.frame(w.p.values,w.p.values.adj))
}
calculate.pvalues2 <- function (values,idx1,idx2,paired=TRUE){
  # Apply Wilcoxon test in order to calculate the p-values
  w.p.values <- unlist(mclapply(values,function(probe) {
    x <- coin::wilcoxsign_test(as.matrix(probe[idx1]) ~ as.matrix(probe[idx2]),
                               zero.method="Wilcoxon",
                               dist="exact")
    z <- coin::pvalue(x)
    return(z)
  }, mc.cores=detectCores()))
  print(w.p.values)
  ##Plot a histogram
  png(filename="histogram_pvalues.png")
  hist(w.p.values)
  dev.off()
  ##Calculate the adjusted p-values by using the Benjamini-Hochberg (BH) method
  w.p.values.adj <- p.adjust(w.p.values,method="BH")
  print(length(w.p.values.adj))

  png(filename="histogram_pvalues_adj.png")
  ##Plot a histogram
  hist(w.p.values.adj)
  dev.off()

  return(data.frame(w.p.values,w.p.values.adj))
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
#' @import ggplot2
#' @export
volcano.plot <- function(data,
                         filename="volcano.pdf",
                         ylab   = "-1*log10 of the Significance",
                         xlab   = "DNA Methylation",
                         title  = "Volcano Plot TCGA GBM Tumors",
                         legend ="Legend",
                         color = c("1"="black","2"="green","3"="red"),
                         label = c("1"="Not Significant",
                                   "2"="Hypermethylated",
                                   "3"="Hypomethylated"),
                         xlim=NULL,
                         ylim=NULL,
                         p.cut=0.05,
                         diffmean.cut = 0
){
  .e <- environment()
  data$threshold <- "1"

  # get significant data
  data.s  <- subset(data,p.value.adj < p.cut)

  # hypermethylated samples compared to old state
  hyper <- subset(data.s,diffmean < (-diffmean.cut))
  data[rownames(hyper),"threshold"] <- "2"

  # hypomethylated samples compared to old state
  hypo <- subset(data.s,diffmean > diffmean.cut)
  data[rownames(hypo),"threshold"] <- "3"

  # Plot a volcano plot
  p <- ggplot(data=data,aes(x=diffmean,y=-1*log10(p.value.adj),colour=threshold),environment = .e) + geom_point()
  if(!is.null(xlim)) {p <- p + xlim(xlim)}
  if(!is.null(ylim)) {p <- p + ylim(ylim)}
  p <- p + labs(title=title) + ylab(ylab) + xlab(xlab)
  p <- p + geom_vline(aes(xintercept= -diffmean.cut), colour="black", linetype = "dashed") +
    geom_vline(aes(xintercept=  diffmean.cut), colour="black", linetype = "dashed") +
    geom_hline(aes(yintercept=  -1*log10(p.cut)),colour="black", linetype = "dashed")
  p <- p + scale_color_manual(breaks=c("1","2","3"),
                              values = color,
                              labels = label,
                              name   = legend)
  # saving box plot to analyse it
  ggsave(p, filename=filename, width = 10, heigh = 5, dpi = 600)
  data <- subset(data,threshold == "2" | threshold == "3" )
  return(data[,c("Composite.Element.REF","threshold")])
}


#  Get latest Genome Reference Consortium Human Build
#  And save it as Genomic Ranges
#' @import biomaRt
get.GRCh.bioMart <- function(){
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  chrom <- c(1:22,"X","Y")
  gene.location <-getBM(attributes= c("ensembl_gene_id","ensembl_transcript_id","chromosome_name",
                                      "start_position","end_position", "strand","external_gene_name","external_transcript_name",
                                      "external_gene_source","external_transcript_source_name","hgnc_id","entrezgene"),
                        filters=c("chromosome_name"),
                        values=list(chrom), mart=ensembl)

  gene.location <- gene.location[!is.na (gene.location$entrezgene),]
  #the duplicates are different transcripts, not different coordinates
  gene.location <- gene.location[!duplicated(gene.location$entrezgene),]

  save(gene.location,file = "inst/extdata/GRCh.rda")
}

#' @title Preparing methylation data and expression data
#'
#' @description
#'  Preparing methylation data and expression data
#'
#' @param met methylation data
#' @param volcano data
#' @param expression expression data
#' @return Methylation table
#' @import GenomicRanges
#' @export
starbursanalysis <- function(met,expression){
  ####fix methylation gene names before merging.
  ### map gene ID to genomic coordinates
  gene.GR <- GRanges(seqnames = paste0("chr",gene.location$chromosome_name),
                     ranges   = IRanges(start = gene.location$start_position,
                                        end   = gene.location$end_position),
                     strand   = gene.location$strand,
                     symbol   = gene.location$external_gene_name,
                     EntrezID = gene.location$entrezgene
  )

  probe.info <- GRanges(seqnames = paste0("chr",met$Chromosome),
                        ranges = IRanges(start = met$Genomic_Coordinate,
                                         end = met$Genomic_Coordinate),
                        probeID = met$Composite.Element.REF)

  distance <- as.data.frame(distanceToNearest(probe.info,gene.GR)) #closest gene to each 450k probe ##your data
  rownames(gene.location) <- NULL
  gene.order.by.distance <- gene.location[distance$subjectHits,]
  gene.order.by.distance$distance <- as.matrix(distance$distance)
  data.nearest.gene <- cbind(met, gene.order.by.distance[,c("external_gene_name","entrezgene","distance")])

  #volcano$gene <- met$Gene_Symbol
  #volcano$cgID <- rownames(volcano)
  volcano <- merge(met, expression, by.x = "Gene_Symbol", by.y = "GeneSymbol", incomparables = NA)
  volcano$ID <- paste(volcano$Gene_Symbol, volcano$probeID, volcano$cgID, sep = ".")

  # Preparing gene expression
  volcano$geFDR <- log10(volcano$FDR)
  volcano$geFDR2 <- volcano$geFDR
  volcano[volcano$dm > 0, "geFDR2"] <- -1 * volcano[volcano$dm > 0, "geFDR"]

  # Preparing methylation
  volcano$meFDR <- log10(volcano$p.value.adj)
  volcano$meFDR2 <- volcano$meFDR
  volcano[volcano$diffmean > 0, "meFDR2"] <- -1 * volcano[volcano$diffmean > 0, "meFDR"]

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
#' @param title main title
#' @param ylab y axis text
#' @param xlab x axis text
#' @param xlim x limits to cut image
#' @param ylim y limits to cut image
#' @param p.cut p value cut-off
#' @param diffmean.cut diffmean cut-off
#' @import ggplot2
#' @import GenomicRanges
#' @export
starburstplot <- function(data,
                          filename="volcano.pdf",
                          ylab   = "Gene Expression\nlog10 of the adjusted Significance (FDR)",
                          xlab   = "DNA Methylation\nlog10 of the adjusted Significance (FDR)",
                          title  = "Starburst Plot",
                          legend = "Methylation/Expression Relation",
                          color = c("1"="black","2"="purple","3"="darkgreen","4"="blue","5"="darkred","6"="red","7"="green"),
                          label = c("1"="Not Significant","2"="Up & Hypo","3"="Down & Hypo","4"= "hypo","5"="hyper","6"="Up","7"="Down"),
                          xlim=NULL,
                          ylim=NULL,
                          p.cut=0.05,
                          diffmean.cut = 0
){
  .e <- environment()
  volcano.m <- data
  volcano.m$threshold.starburst <- "1"
  volcano.m$threshold.size <- "1"

  # subseting by regulation (geFDR) and methylation level (meFDR)
  #   down regulated                up regulated
  #     lowerThr  ||||||||||||||||    upperThr
  #   hypomethylated               hipermethylated
  lowerThr <- log10(0.05) # -1.30103
  upperThr <- (-lowerThr) # +1.30103

  # Group 2:up regulated and hypomethylated
  a <- subset(volcano.m,geFDR2 > upperThr & meFDR2 < lowerThr)

  # Group 3: down regulated and hypomethylated
  b <- subset(volcano.m, geFDR2 < lowerThr & meFDR2 < lowerThr)

  # Group 4: hypomethylated
  c <- subset(volcano.m, geFDR2 > lowerThr & geFDR2 < upperThr & meFDR2 < lowerThr)

  # Group 5: hypermethylated
  d <- subset(volcano.m,geFDR2 > lowerThr & geFDR2 < upperThr & meFDR2 > upperThr)

  # Group 6: upregulated
  e <- subset(volcano.m,geFDR2 > upperThr & meFDR2 < upperThr & meFDR2 > lowerThr)

  # Group 7: downregulated
  f <- subset(volcano.m,geFDR2 < lowerThr & meFDR2 < upperThr & meFDR2 > lowerThr)

  size <- c("3","3","2","2","2","2")
  groups <- c("2","3","4","5","6","7")
  s <- list(a,b,c,d,e,f)
  for(i in seq_along(s)){
    idx <- rownames(s[[i]])
    if(length(idx)>0){
      volcano.m[idx,"threshold.starburst"] <- groups[i]
      volcano.m[idx,"threshold.size"] <-  size [i]
    }
  }

  ##starburst plot
  p <- ggplot(data=volcano.m, environment = .e,
              aes(x=meFDR2, y=geFDR2, colour=threshold.starburst, size = threshold.size)
  ) + geom_point()
  if(!is.null(xlim)) {p <- p + xlim(xlim)}
  if(!is.null(ylim)) {p <- p + ylim(ylim)}
  p <- p + labs(title) + ylab(ylab) + xlab(xlab)
  p <- p + scale_color_manual(
    values=color,
    labels=label,
    name=legend
  )
  p + geom_hline( aes( yintercept = lowerThr))+
    geom_hline( aes( yintercept = upperThr)) +
    geom_vline( aes( xintercept = lowerThr)) +
    geom_vline( aes( xintercept = upperThr))

  ggsave(file = "starbust.gcimp.pdf", width=14, height = 10)

  # return methylation < 0, expressao >0
}

#' @import matlab gplots
heatmap.plus.sm <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
                             distfun = dist, hclustfun = hclust, reorderfun = function(d,
                                                                                       w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv,
                                                                                                                                                  "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,
                             margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 +
                               1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                             labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,
                             verbose = getOption("verbose"), breaks, key = TRUE, ...)
{
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (is.null(Rowv))
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv))
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram"))
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv)
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr)))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram"))
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv)
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc)))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow))
    if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol))
    if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
    sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
    sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0,
            4)
  if (!missing(ColSideColors)) {
    if (!is.matrix(ColSideColors))
      stop("'ColSideColors' must be a matrix")
    if (!is.character(ColSideColors) || dim(ColSideColors)[1] !=
          nc)
      stop("'ColSideColors' dim()[2] must be of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1], 0.6, lhei[2])
  }
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors))
      stop("'RowSideColors' must be a matrix")
    if (!is.character(RowSideColors) || dim(RowSideColors)[1] !=
          nr)
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1], 0.2, lwid[2])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei,
        "; lmat=\n")
    print(lmat)
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    rsc = RowSideColors[rowInd, ]
    rsc.colors = matrix()
    rsc.names = names(table(rsc))
    rsc.i = 1
    for (rsc.name in rsc.names) {
      rsc.colors[rsc.i] = rsc.name
      rsc[rsc == rsc.name] = rsc.i
      rsc.i = rsc.i + 1
    }
    rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
    image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
    if (length(colnames(RowSideColors)) > 0) {
      axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors),
           las = 2, tick = FALSE)
    }
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    csc = ColSideColors[colInd, ]
    csc.colors = matrix()
    csc.names = names(table(csc))
    csc.i = 1
    for (csc.name in csc.names) {
      csc.colors[csc.i] = csc.name
      csc[csc == csc.name] = csc.i
      csc.i = csc.i + 1
    }
    csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
    image(csc, col = as.vector(csc.colors), axes = FALSE)
    if (length(colnames(ColSideColors)) > 0) {
      axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors),
           las = 2, tick = FALSE)
    }
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
  }
  if (revC) {
    iy <- nr:1
    ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  par(mar = c(margins[1], 0, 0, 0))
  if (doRdend)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
  if (doCdend)
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main))
    frame()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
                                                                doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
  ##################KEY###############
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x))
      tmpbreaks[length(tmpbreaks)] <- max(abs(x))
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, "Value", line = 2)
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }

    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)

}
