samples <- c("TCGA-06-0125-01A-01D-A45W-05","TCGA-06-0152-01A-02D-A45W-05",
             "TCGA-06-0171-01A-02D-A45W-05","TCGA-06-0190-01A-01D-A45W-05",
             "TCGA-06-0210-01A-01D-A45W-05","TCGA-06-0211-01A-01D-A45W-05",
             "TCGA-06-0221-01A-01D-A45W-05","TCGA-14-0736-01A-01D-A45W-05",
             "TCGA-19-0957-01C-01D-A45W-05","TCGA-19-1389-01A-01D-A45W-05",
             "TCGA-14-1402-01A-01D-A45W-05","TCGA-19-4065-01A-01D-2004-05",
             "TCGA-06-0125-02A-11D-2004-05","TCGA-06-0152-02A-01D-2004-05",
             "TCGA-06-0171-02A-11D-2004-05","TCGA-06-0190-02A-01D-2004-05",
             "TCGA-06-0210-02A-01D-2004-05","TCGA-06-0211-02A-02D-2004-05",
             "TCGA-06-0221-02A-11D-2004-05","TCGA-14-0736-02A-01D-2004-05",
             "TCGA-19-0957-02A-11D-2004-05","TCGA-19-1389-02A-21D-2004-05",
             "TCGA-14-1402-02A-01D-2004-05","TCGA-19-4065-02A-11D-2004-05")
query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",
                   level = "3", samples = samples)
TCGADownload(query,path = "dataDemo2")
query <- TCGAQuery(tumor = "gbm", platform = "bio", level = "2")
TCGADownload(query, path = "dataDemo2")

#-------------------------- Preparing data and metadata
# met: data frame with probe info (probe id, ge) and beta values
# beta: data frame with only beta values (probe vs patient)
# probe: data frame with only probe info(probe vs patient)
# Conclusion: met = probe +d beta
# met.md = patient metadata
met <- organizeMethylationDataFrame("dataDemo2")
probe  <- met[,1:4]
beta   <- met[,5:ncol(met)]

met.md <- organizeMethylationMetaDataFrame("dataDemo2")
samples <- colnames(beta)
idx <- is.element(met.md$bcr_patient_barcode,strtrim(samples,12))
met.md <- met.md[idx,]

pat.prim <- "TCGA-[0-9a-z]{2}-[0-9a-z]{4}-01"
pat.rec <- "TCGA-[0-9a-z]{2}-[0-9a-z]{4}-02"
rec <- grep(pat.rec,colnames(beta))
prim <- grep(pat.prim,colnames(beta))
beta.t <-  data.frame(t(beta))
pvalues <- calculate.pvalues(beta.t,prim,rec)

probe$p.value <- pvalues[,1]
probe$p.value.adj <- pvalues[,2]
prim.rec <- diffmean(beta[1:2,], beta[3:4,])
probe <- cbind(probe,prim.rec)
hypo.hyper <- volcanoPlot(probe,p.cut = 0.05,diffmean.cut = 0.2)
