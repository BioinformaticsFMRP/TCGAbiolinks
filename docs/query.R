## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ---- eval = TRUE, echo = FALSE------------------------------------------
datatable(TCGAbiolinks:::getGDCprojects(),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
          rownames = FALSE,
          caption = "List of projects")

## ---- eval = TRUE, echo = FALSE------------------------------------------
datatable(TCGAbiolinks:::getBarcodeDefinition(),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
          rownames = FALSE,
          caption = "List sample types")

## ------------------------------------------------------------------------
datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454"),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

## ------------------------------------------------------------------------
datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=1817673686"),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 40), 
          rownames = FALSE)

## ----message=FALSE, warning=FALSE----------------------------------------
query <- GDCquery(project = c("TCGA-GBM", "TCGA-LGG"),
                  data.category = "DNA Methylation",
                  legacy = FALSE,
                  platform = c("Illumina Human Methylation 450"),
                  sample.type = "Recurrent Solid Tumor"
)
datatable(getResults(query), 
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)

## ----message=FALSE, warning=FALSE----------------------------------------
query.met <- GDCquery(project = "TCGA-COAD",
                  data.category = "DNA Methylation",
                  legacy = FALSE,
                  platform = c("Illumina Human Methylation 450"))
query.exp <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(substr(getResults(query.met, cols = "cases"), 1, 12),
                             substr(getResults(query.exp, cols = "cases"), 1, 12))

# Only seelct the first 5 patients
query.met <- GDCquery(project = "TCGA-COAD",
                  data.category = "DNA Methylation",
                  legacy = FALSE,
                  platform = c("Illumina Human Methylation 450"),
                  barcode = common.patients[1:5])
query.exp <- GDCquery(project = "TCGA-COAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ",
                  barcode = common.patients[1:5])
datatable(getResults(query.met, cols = c("data_type","cases")),
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)
datatable(getResults(query.exp, cols = c("data_type","cases")), 
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)



## ----message=FALSE, warning=FALSE----------------------------------------
query <- GDCquery(project = c("TCGA-BRCA"),
                  data.category = "Raw Sequencing Data",  
                  sample.type = "Primary solid Tumor")
# Only first 100 to make render faster
datatable(getResults(query, rows = 1:100,cols = c("file_name","cases")), 
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)

## ----message=FALSE, warning=FALSE----------------------------------------
query <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"),
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = c("Illumina Human Methylation 450", "Illumina Human Methylation 27"))
datatable(getResults(query, rows = 1:100), 
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)

## ----message=FALSE, warning=FALSE----------------------------------------
# Gene expression aligned against hg19.
query.exp.hg19 <- GDCquery(project = "TCGA-GBM",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
                  legacy = TRUE)
datatable(getResults(query.exp.hg19), 
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)

