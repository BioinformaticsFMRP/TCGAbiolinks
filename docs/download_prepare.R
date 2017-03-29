## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results = 'hide', message=FALSE, warning=FALSE----------------------
query <- GDCquery(project = "TCGA-GBM",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           file.type  = "normalized_results",
                           experimental.strategy = "RNA-Seq",
                           barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
                           legacy = TRUE)
GDCdownload(query, method = "api", chunks.per.download = 10)
data <- GDCprepare(query)

## ----message=FALSE, warning=FALSE----------------------------------------
# Gene expression aligned against hg19.
datatable(as.data.frame(colData(data)), 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)
# Only first 100 to make render faster
datatable(assay(data)[1:100,], 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = TRUE)

rowRanges(data)

## ----results = 'hide', message=FALSE, warning=FALSE----------------------
# Gene expression aligned against hg38
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
tryCatch(GDCdownload(query, method = "client"),error = function(e)GDCdownload(query))
data <- GDCprepare(query)

## ----message=FALSE, warning=FALSE----------------------------------------
datatable(as.data.frame(colData(data)), 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)

datatable(assay(data)[1:100,], 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = TRUE)

rowRanges(data)

