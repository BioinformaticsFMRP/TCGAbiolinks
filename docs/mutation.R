## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
acc.maf <- GDCquery_Maf("ACC", pipelines = "muse")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
# Only first 50 to make render faster
datatable(acc.maf[1:50,],
              filter = 'top',
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
# Check maf availables
datatable(select(getResults(query.maf.hg19),-contains("cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
          rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
query.maf.hg19 <- GDCquery(project = "TCGA-CHOL", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.maf",
                           legacy = TRUE)
GDCdownload(query.maf.hg19)
maf <- GDCprepare(query.maf.hg19)

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
# Only first 50 to make render faster
datatable(maf[1:50,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

