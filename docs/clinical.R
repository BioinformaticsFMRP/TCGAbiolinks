## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results='hide', echo=TRUE, message=FALSE, warning=FALSE-------------
clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
datatable(clinical, filter = 'top', 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
          rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
datatable(clinical, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
datatable(clinical.drug, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
datatable(clinical.radiation, options = list(scrollX = TRUE,  keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
datatable(clinical.admin, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Tissue slide image files
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Tissue slide image",
                  legacy = TRUE,
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")) 

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
query %>% getResults %>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Pathology report
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Pathology report",
                  legacy = TRUE,
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))  

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
query %>% getResults %>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Tissue slide image
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Tissue slide image",
                  legacy = TRUE,
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")) 

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
query %>% getResults %>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Clinical Supplement
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Clinical Supplement",
                  legacy = TRUE,
                  barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")) 

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
query %>% getResults %>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Clinical data
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Clinical data",
                  legacy = TRUE,
                  file.type = "txt")  

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
query %>% getResults %>% select(-matches("cases"))%>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
GDCdownload(query)
clinical.biotab <- GDCprepare(query)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
names(clinical.biotab)
datatable(clinical.biotab$clinical_radiation_coad)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Get XML files and parse them
clin.query <- GDCquery(project = "TCGA-READ", data.category = "Clinical", barcode = "TCGA-F5-6702")
GDCdownload(clin.query)
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")

# Get indexed data
clinical.index <- GDCquery_clinic("TCGA-READ")

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
select(clinical.patient,vital_status,days_to_death,days_to_last_followup) %>% datatable
select(clinical.patient.followup, vital_status,days_to_death,days_to_last_followup) %>% datatable
# Vital status should be the same in the follow up table 
filter(clinical.index,submitter_id == "TCGA-F5-6702") %>% select(vital_status,days_to_death,days_to_last_follow_up) %>% datatable

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE-----------
# Get XML files and parse them
recurrent.samples <- GDCquery(project = "TCGA-LIHC",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts",
                             sample.type = 	"Recurrent Solid Tumor")$results[[1]] %>% select(cases)
recurrent.patients <- unique(substr(recurrent.samples$cases,1,12))
clin.query <- GDCquery(project = "TCGA-LIHC", data.category = "Clinical", barcode = recurrent.patients)
GDCdownload(clin.query)
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient") 

## ----echo = TRUE, message = FALSE, warning = FALSE-----------------------
# Get indexed data
GDCquery_clinic("TCGA-LIHC") %>% filter(submitter_id %in% recurrent.patients) %>% 
    select(progression_or_recurrence,days_to_recurrence,tumor_grade) %>% datatable

# XML data
clinical.patient %>% select(bcr_patient_barcode,neoplasm_histologic_grade) %>% datatable


