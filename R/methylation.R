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
  prim$mean.primary <- apply(prim,1,mean,na.rm=TRUE)
  rec$mean.recurrence <- apply(rec,1,mean,na.rm=TRUE)
  prim$probeID <- rownames(prim)
  rec$probeID <- rownames(rec)
  prim.rec <- merge(prim[,c("probeID","mean.primary")],rec[,c("probeID","mean.recurrence")],by="probeID")
  rownames(prim.rec) <- prim.rec[,1]
  prim.rec$diffmean <- prim.rec$mean.primary-prim.rec$mean.recurrence

  png(filename="histogram_diffmean.png")
  hist(prim.rec$diffmean)
  dev.off()

  return(prim.rec)
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
  files <- files[grep(".*clinical.*",files)]
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
      methylation <- data[, c(1,4,2)]
    }
    else{
      methylation <- merge(methylation, data[,c(1,2)], by="Composite.Element.REF")
    }
  }

  # Use Composite.Element.REF as name, instead of column
  rownames(methylation) <- methylation$Composite.Element.REF
  methylation$Composite.Element.REF <- NULL

  # remove X Y chromossomes
  methylation <- methylation[methylation$Chromosome !="X" & methylation$Chromosome!="Y",]
  methylation$Chromosome <- NULL

  # remove NA lines
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
calculate.pvalues <- function (values,idx1,idx2,paired=TRUE){
  # Apply Wilcoxon test in order to calculate the p-values
  w.p.values <- unlist(mclapply(values,function(probe) {
    zz <- wilcox.exact(as.matrix(probe[idx1]),as.matrix(probe[idx2]), exact=TRUE,paired=paired)
    z <- zz$p.value
    return(z)
  }, mc.cores=detectCores()))
  ##Plot a histogram
  print(w.p.values)
  png(filename="histogram_pvalues.png")
  hist(w.p.values)
  dev.off()
  print(length(w.p.values))
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
                         color = c("green","red","purple","orange","salmon","grey"),
                         label = c("Not Significant",
                                   "Hypermethylated in Recurrent GBM",
                                   "Hypomethylated in Recurrent GBM"),
                         xlim=NULL,
                         ylim=NULL,
                         p.cut=0.05
){

  data$threshold <- "1"

  # get significant data
  data.s  <- subset(data,p.value.adj < p.cut)

  # hypermethylated samples compared to old state
  hyper <- subset(data.s,diffmean<0)
  data[rownames(hyper),"threshold"] <- "2"

  # hypomethylated samples compared to old state
  hypo <- subset(data.s,diffmean>0)
  data[rownames(hypo),"threshold"] <- "3"

  # Plot a volcano plot
  p <- ggplot(data=data,aes(x=diffmean,y=-1*log10(p.value.adj),colour=threshold)) + geom_point()
  if(!is.null(xlim)) {p <- p + xlim(xlim)}
  if(!is.null(ylim)) {p <- p + ylim(ylim)}
  p <- p + labs(title=title)  + ylab(ylab) + xlab(xlab)
  p <- p + scale_color_manual(breaks=c("1","2","3"),
                              values = color,
                              labels = label,
                              name   = legend)
  # saving box plot to analyse it
  ggsave(p, filename=filename, width = 10, heigh = 5, dpi = 600)
}

##Starburst
# Note: Received two files from Michele ('DEGs_agi_gcimp.txt' and 'DEGs_affy_gcimp.txt').
# He performed a supervised analysis using gene expression platforms for the same samples.

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
#' @import GenomicRanges
#' @export
starburstAnalysis <- function(wd = NULL){

  gc.affy <- read.delim(file = "/dados/ResearchProjects/tathi//LGG.GBM_Tathi/DEGs_affy_gcimp.txt", sep = " ")
  gc.affy$probeID <- rownames(gc.affy)
  gc.agi <- read.delim(file = "/dados/ResearchProjects/tathi//LGG.GBM_Tathi/DEGs_agi_gcimp.txt", sep = " ")
  gc.agi$probeID <- rownames(gc.agi)

  ####fix methylation gene names before merging.
  ### map gene ID to genomic coordinates
  gene.location <- read.delim("/dados/ResearchProjects/thais/TCGA/Normals/mRNA_Brain/Galaxy1-[Homo_sapiens_genes_(GRCh37.p13)].tabular")
  gene.location <- gene.location[gene.location$Chromosome.Name %in% c(1:22,"X","Y"),]
  gene.location <- gene.location[!is.na (gene.location$EntrezGene.ID),]
  gene.location <- gene.location[!duplicated(gene.location$EntrezGene.ID),] #the duplicates are different transcripts, not different coordinates

  gene.GR <- GRanges(seqnames = paste0("chr",gene.location$Chromosome.Name),
                     ranges   = IRanges(start = gene.location$Gene.Start..bp.,
                                        end   = gene.location$Gene.End..bp.),
                     strand   = gene.location$Strand,
                     symbol   = gene.location$Associated.Gene.Name ,
                     EntrezID = gene.location$EntrezGene.ID
  )
  probe.info <- GRanges(seqnames = paste0("chr",LGG.GBM.250914$Chromosome), ranges = IRanges(start = LGG.GBM.250914$Genomic_Coor, end = LGG.GBM.250914$Genomic_Coor), probeID = LGG.GBM.250914$Composite.El)
  distance <- as.data.frame(distanceToNearest(probe.info ,gene.GR)) #closest gene to each 450k probe ##your data
  rownames(gene.location) <- NULL
  gene.order.by.distance <- gene.location[distance$subjectHits,]
  gene.order.by.distance$distance <- as.matrix(distance$distance)
  GBM.LGG.27.450k.nearest.gene <- cbind(LGG.GBM.250914, gene.order.by.distance[,c("Associated.Gene.Name ","EntrezGene.ID","distance")])

  all(rownames(volcano) == LGG.GBM.250914$Composite.Element.REF)
  volcano$gene <- LGG.GBM.250914$Gene_Symbol
  volcano$cgID <- rownames(volcano)
  volcano.m <- merge(volcano, gc.affy, by.x = "gene", by.y = "GeneSymbol", incomparables = NA)
  volcano.m$ID <- paste(volcano.m$gene, volcano.m$probeID, volcano.m$cgID, sep = ".")
  volcano.m$geFDR <- log10(volcano.m$FDR)
  volcano.m$geFDR2 <- volcano.m$geFDR
  volcano.m$meFDR <- log10(volcano.m$p.value.adj.LGm1xLGm2.3)
  volcano.m[volcano.m$dm > 0, "geFDR2"] <- -1 * volcano.m[volcano.m$dm > 0, "geFDR"]
  volcano.m$meFDR2 <- volcano.m$meFDR
  volcano.m[volcano.m$diffmean.LGm1_LGm2.3 > 0, "meFDR2"] <- -1 * volcano.m[volcano.m$diffmean.LGm1_LGm2.3 > 0, "meFDR"]
  volcano.m$threshold.starburst <- "1"
  volcano.m$threshold.size <- "1"

  # subseting by regulation (geFDR) and methylation level (meFDR)
  #   down regulated                up regulated
  #     lowerThr  ||||||||||||||||    upperThr
  #   hypomethylated               hipermethylated
  lowerThr <- log10(0.05) # -1.30103
  upperThr <- (-lowerThr) # +1.30103

  # up regulated and hypomethylated in lgm1
  a <- subset(volcano.m,geFDR2 > upperThr & meFDR2 < lowerThr)

  # down regulated and hypomethylated in lgm1
  b <- subset(volcano.m, geFDR2 < lowerThr & meFDR2 < lowerThr)

  # hypomethylated in lgm1
  c <- subset(volcano.m, geFDR2 > lowerThr & geFDR2 < upperThr & meFDR2 < lowerThr)

  # hypermethylated in lgm1
  d <- subset(volcano.m,geFDR2 > lowerThr & geFDR2 < upperThr & meFDR2 > upperThr)

  # upregulated in lmg1
  e <- subset(volcano.m,geFDR2 > upperThr &meFDR2 < upperThr & meFDR2 > lowerThr)

  # downregulated in lgm1
  f <- subset(volcano.m,geFDR2 < lowerThr & meFDR2 < upperThr & meFDR2 > lowerThr)

  tsSize      <- c("3","3","2","2","2","2")
  tsStarburst <- c("2","3","4","5","6","7")
  s <- c(a,b,c,d,e,f)
  mapply(function(x,y,z){
    volcano.m[rownames(x),"threshold.starburst"] <- y
    volcano.m[rownames(x),"threshold.size"] <- z
  },
  x = s,
  y = tsStarburst,
  z = tsSize,
  )

  ##starburst plot
  p <- ggplot(data=volcano.m,
              aes(x=meFDR2, y=geFDR2, colour=threshold.starburst, size = threshold.size)
  ) +
    geom_point() +
    #scale_color_manual(values = c("black", "red", "green")) +
    xlim(c(-2.5,2.5)) + ylim(c(-3.5,3.5)) +
    xlab("DNA Methylation\nlog10 of the adjusted Significance (FDR)") +
    ylab("Gene Expression\nlog10 of the adjusted Significance (FDR)") +
    labs(title = "Starburst Plot\nLGm1-G-CIMP vs. LGm2+LGm3-G-CIMP")
  p <- p + scale_color_manual(
    breaks=c( "1","2","3","4", "5", "6","7"), # color scale (for points)
    values=c("black", "purple", "darkgreen","blue","darkred","red","green"),
    labels=c("Not Significant","Up & Hypo","Down & Hypo", "hypo", "hyper", "Up", "Down"),
    name="Relation to LGm1"
  )
  p <- p + scale_size_manual(name="Relation to LGm1",
                             values=c(1,3,4),
                             labels=c("Not Significant",
                                      "Significant",
                                      "Significantly Epigenetically Activate")
  )
  p + geom_hline( aes( yintercept = lowerThr)) +
    geom_hline( aes( yintercept = upperThr)) +
    geom_vline( aes( xintercept = lowerThr)) +
    geom_vline( aes( xintercept = upperThr))

  ggsave(file = "starbust.gcimp.pdf", width=14, height = 10)
}
