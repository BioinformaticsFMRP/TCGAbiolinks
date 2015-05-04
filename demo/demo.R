library(downloader)
library(survival)
library(GGally)
library(ggplot2)
library(coin)
library(parallel)
#------------------ Downloading tcga data ----------------------------
dir.create("data")
query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",level="3")
TCGADownload(query[1:2,],path="data")

query <- TCGAQuery(tumor = "gbm", platform = "bio", level="2")
TCGADownload(query,path="data")

#-------------------------- Preparing data and metadata
all.met <- organizeMethylationDataFrame("data")
met <- all.met[,5:ncol(all.met)]
met.md <- organizeMethylationMetaDataFrame("data")
samples <- colnames(met)
idx <- is.element(met.md$bcr_patient_barcode,strtrim(samples,12))
met.md <- met.md[idx,]

#---- random split of pacients into groups -----
met.md$cluster <- c(rep("group1",nrow(met.md)/4),rep("group2",nrow(met.md)/4),
                    rep("group3",nrow(met.md)/4),rep("group4",nrow(met.md)-3*(floor(nrow(met.md)/4))))
rownames(met.md) <- met.md$bcr_patient_barcode

#---------------------- survival
survivalPlot(met.md)

#----------------------mean methylation
met.mean <- data.frame(apply(met,2,mean,na.rm=TRUE))
colnames(met.mean) <- "avg"
met.mean$patient <-  strtrim(rownames(met.mean),12)
aux <- merge(met.mean,met.md, by.x="patient",by.y="bcr_patient_barcode")
met.mean.boxplot(aux)

#----------------------- calculate.pvalues (just an example)
met.t <-  data.frame(t(met))
#usuário   sistema decorrido
#229.528    13.091   108.977
pvalues <- calculate.pvalues(met.t,c(1,3,5,7,9,11),c(2,4,6,8,10,12))
system.time(pvalues <- calculate.pvalues(met.t,c(1,3,5),c(2,4,6)))
system.time(pvalues2 <- calculate.pvalues2(met.t,c(1,3,5),c(2,4,6)))

#------------- volcano plot
met$p.value <- pvalues[,1]
met$p.value.adj <-pvalues[,2]
prim.rec <- diffmean.prim.rec(met)
volcano <- merge(met,prim.rec[,c(1,4)],by=0)
volcano.plot(volcano,p.cut=0.05)
