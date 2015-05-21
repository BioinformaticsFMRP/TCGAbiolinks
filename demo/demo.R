#------------------ Downloading tcga data ----------------------------
dir.create("dataDemo")
query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",level = "3")
TCGADownload(query,path = "data")

query <- TCGAQuery(tumor = "gbm", platform = "bio", level = "2")
TCGADownload(query,path = "data")

#-------------------------- Preparing data and metadata
# met: data frame with probe info (probe id, ge) and beta values
# beta: data frame with only beta values (probe vs patient)
# probe: data frame with only probe info(probe vs patient)
# Conclusion: met = probe +d beta
# met.md = patient metadata
met <- organizeMethylationDataFrame("data")
probe  <- met[,1:4]
beta   <- met[,5:ncol(met)]

met.md <- organizeMethylationMetaDataFrame("data")
samples <- colnames(beta)
idx <- is.element(met.md$bcr_patient_barcode,strtrim(samples,12))
met.md <- met.md[idx,]

# random split of pacients into groups
met.md$cluster <- c(rep("group1",nrow(met.md)/4),
                    rep("group2",nrow(met.md)/4),
                    rep("group3",nrow(met.md)/4),
                    rep("group4",nrow(met.md)-3*(floor(nrow(met.md)/4))))
rownames(met.md) <- met.md$bcr_patient_barcode

#---------------------- survival
survivalPlot(met.md)

#----------------------mean methylation
met.mean <- data.frame(apply(beta,2,mean,na.rm=TRUE))
colnames(met.mean) <- "avg"
met.mean$patient <-  strtrim(rownames(met.mean),12)
aux <- merge(met.mean,met.md, by.x="patient",by.y="bcr_patient_barcode")
met.mean.boxplot(aux)

#----------------------- calculate.pvalues (just an example)
beta.t <-  data.frame(t(beta))
# TODO verify class of the object
pvalues <- calculate.pvalues(beta.t,c(1,3,5,7,9,11),c(2,4,6,8,10,12))

#------------- volcano plot
probe$p.value <- pvalues[,1]
probe$p.value.adj <- pvalues[,2]
prim.rec <- diffmean.prim.rec(beta)
probe <- cbind(probe,prim.rec)
hypo.hyper <- volcano.plot(probe,p.cut = 0.619)

#------------------- heatmap
#get bvalues, for the probes hypo/hyper
hypo.hyper.probe <- beta[rownames(hypo.hyper),]
hypo.hyper.probe <- hypo.hyper.probe[1:10,]
heatmap.plus.sm(as.matrix(hypo.hyper.probe),
                Colv = NA,
                Rowv = NA,
                scale = "none", col = jet.colors(75))
dev.off()

heatmap.2 (as.matrix(hypo.hyper.probe))
dev.off()

#----- starburst plot
gc.affy <- read.delim(file = "DEGs_affy_gcimp.txt", sep = " ")
gc.affy$probeID <- rownames(gc.affy)
gene.met <- starbursanalysis(probe,gc.affy)
starburstplot(gene.met)
