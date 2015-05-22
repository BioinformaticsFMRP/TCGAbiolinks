#-------------------------- Preparing data and metadata
# met: data frame with probe info (probe id, ge) and beta values
# beta: data frame with only beta values (probe vs patient)
# probe: data frame with only probe info(probe vs patient)
# Conclusion: met = probe +d beta
# met.md = patient metadata
sample <- c("TCGA-06-0147-01A-01D-0218-05", "06-0145-01A-06D-0218-05")
query <- TCGAQuery(tumor = "GBM",samples = sample, level = 3)
TCGADownload(query,path = "exampleData",samples = sample, quiet = TRUE)
met <- TCGAPrepare(query, dir="exampleData")
probe <- met[,1:4]
beta  <- met[,5:ncol(met)]

query <- TCGAQuery(tumor = "GBM",platform = "bio",level = 2)
TCGADownload(query,path = "exampleMetaData",type="clinical_patient",
             quiet = TRUE)
met.md <- TCGAPrepare(query,"exampleMetaData")
met.md <- met.md[-(1:2),]
samples <- colnames(beta)
idx <- is.element(met.md$bcr_patient_barcode,strtrim(samples,12))
met.md <- met.md[idx,]

# random split of pacients into groups
met.md$cluster <- c("group1","group2")
rownames(met.md) <- met.md$bcr_patient_barcode

#---------------------- survival
survivalPlot(met.md)

#----------------------mean methylation
met.mean <- data.frame(apply(beta, 2, mean, na.rm = TRUE))
colnames(met.mean) <- "avg"
met.mean$patient <-  strtrim(rownames(met.mean),12)
aux <- merge(met.mean,met.md,
             by.x = "patient",
             by.y = "bcr_patient_barcode")
metMeanBoxplot(aux)

#----------------------- calculate.pvalues (just an example)
beta.t <-  data.frame(t(beta))
# TODO verify class of the object
pvalues <- calculate.pvalues(beta.t,c(1,3,5,7,9,11),c(2,4,6,8,10,12))

#------------- volcano plot
probe$p.value <- pvalues[,1]
probe$p.value.adj <- pvalues[,2]
prim.rec <- diffmean(beta,c(1,3,5,7,9,11),c(2,4,6,8,10,12))
probe <- cbind(probe,prim.rec)
hypo.hyper <- volcanoPlot(probe,p.cut = 0.619)
