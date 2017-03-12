## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

## ---- echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE----------------
devtools::load_all(".")

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ---- eval = FALSE-------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("TCGAbiolinks")

## ---- eval = FALSE-------------------------------------------------------
#  #---------------------------------------------------------------
#  #  For available entries and combinations please se table below
#  #---------------------------------------------------------------
#  
#  # Gene expression aligned against Hg38
#  query <- GDCquery(project = "TARGET-AML",
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - Counts")
#  
#  # All DNA methylation data for TCGA-GBM and TCGA-GBM
#  query.met <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"),
#                        legacy = TRUE,
#                        data.category = "DNA methylation",
#                        platform = c("Illumina Human Methylation 450", "Illumina Human Methylation 27"))
#  
#  # Using sample type to get only Primary solid Tumor samples and Solid Tissue Normal
#  query.mirna <- GDCquery(project = "TCGA-ACC",
#                          data.category = "Transcriptome Profiling",
#                          data.type = "miRNA Expression Quantification",
#                          sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
#  
#  # Example Using legacy to accessing hg19 and filtering by barcode
#  query <- GDCquery(project = "TCGA-GBM",
#                    data.category = "DNA methylation",
#                    platform = "Illumina Human Methylation 27",
#                    legacy = TRUE,
#                    barcode = c("TCGA-02-0047-01A-01D-0186-05","TCGA-06-2559-01A-01D-0788-05"))
#  
#  # Gene expression aligned against hg19.
#  query.exp.hg19 <- GDCquery(project = "TCGA-GBM",
#                    data.category = "Gene expression",
#                    data.type = "Gene expression quantification",
#                    platform = "Illumina HiSeq",
#                    file.type  = "normalized_results",
#                    experimental.strategy = "RNA-Seq",
#                    barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
#                    legacy = TRUE)
#  
#  # Searching idat file for DNA methylation
#  query <- GDCquery(project = "TCGA-OV",
#                    data.category = "Raw microarray data",
#                    data.type = "Raw intensities",
#                    experimental.strategy = "Methylation array",
#                    legacy = TRUE,
#                    file.type = ".idat",
#                    platform = "Illumina Human Methylation 450")
#  

## ---- eval = TRUE, echo = FALSE------------------------------------------
knitr::kable(getGDCprojects(), digits = 2,
             caption = "List of projects",row.names = FALSE)

## ---- eval = TRUE, echo = FALSE------------------------------------------
knitr::kable(getBarcodeDefinition(), digits = 2,
             caption = "List sample types",row.names = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  acc.muse.maf <- GDCquery_Maf("ACC", pipelines = "muse")
#  acc.varscan2.maf <- GDCquery_Maf("ACC", pipelines = "varscan2")
#  acc.somaticsniper.maf <- GDCquery_Maf("ACC", pipelines = "somaticsniper")
#  acc.mutect.maf <- GDCquery_Maf("ACC", pipelines = "mutect")

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
#  mut <- GDCquery_Maf("ACC", pipelines = "mutect2")
#  clin <- GDCquery_clinic("TCGA-ACC","clinical")
#  clin <- clin[,c("bcr_patient_barcode","disease","gender","tumor_stage","race","vital_status")]
#  TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:20],
#                          filename = "oncoprint.pdf",
#                          annotation = clin,
#                          color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
#                          rows.font.size=10,
#                          width = 10,
#                          heatmap.legend.side = "right",
#                          dist.col = 0,
#                          label.font.size = 10)

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
library(png)
library(grid)
img <- readPNG("oncoprint.png")
grid.raster(img)

## ---- message=FALSE, warning=FALSE---------------------------------------
query.maf.hg19 <- GDCquery(project = "TCGA-COAD", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)
# Check maf availables
knitr::kable(getResults(query.maf.hg19)[,c("created_datetime","file_name")]) 

## ---- eval = FALSE-------------------------------------------------------
#  query.maf.hg19 <- GDCquery(project = "TCGA-COAD",
#                             data.category = "Simple nucleotide variation",
#                             data.type = "Simple somatic mutation",
#                             access = "open",
#                             file.type = "gsc_COAD_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf",
#                             legacy = TRUE)
#  GDCdownload(query.maf.hg19)
#  coad.mutect.maf <- GDCprepare(query.maf.hg19)

## ---- fig.height=6, message=FALSE, warning=FALSE, include=FALSE----------
clin.query <- GDCquery(project = "TCGA-BLCA", data.category = "Clinical", barcode = "TCGA-FD-A5C0")
 json  <- tryCatch(GDCdownload(clin.query), 
                   error = function(e) GDCdownload(clin.query, method = "client"))
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
clinical.index <- GDCquery_clinic("TCGA-BLCA")

## ---- eval=FALSE, fig.height=6, message=FALSE, warning=FALSE, include=TRUE----
#  clin.query <- GDCquery(project = "TCGA-BLCA", data.category = "Clinical", barcode = "TCGA-FD-A5C0")
#   json  <- tryCatch(GDCdownload(clin.query),
#                     error = function(e) GDCdownload(clin.query, method = "client"))
#  clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
#  clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
#  clinical.index <- GDCquery_clinic("TCGA-BLCA")

## ---- echo=TRUE, message=FALSE, warning=FALSE, fig.height=6--------------
# Example of the second difference:
clinical.patient[,c("vital_status","days_to_death","days_to_last_followup")]
clinical.patient.followup[,c("vital_status","days_to_death","days_to_last_followup")]
# indexed data is equivalent to follow ups information
clinical.index[clinical.index$submitter_id=="TCGA-FD-A5C0",
               c("vital_status","days_to_death","days_to_last_follow_up")]

## ---- eval = FALSE-------------------------------------------------------
#  clin <- GDCquery_clinic("TCGA-ACC", type = "clinical", save.csv = TRUE)
#  clin <- GDCquery_clinic("TCGA-ACC", type = "biospecimen", save.csv = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  query <- GDCquery(project = "TCGA-COAD",
#                    data.category = "Clinical",
#                    barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
#  GDCdownload(query)
#  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
#  clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
#  clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
#  clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")
#  
#  query <- GDCquery(project = "TCGA-COAD",
#                    data.category = "Biospecimen",
#                    barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
#  GDCdownload(query)
#  clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")
#  clinical.sample <- GDCprepare_clinic(query, clinical.info = "sample")
#  clinical.slide <- GDCprepare_clinic(query, clinical.info = "slide")
#  clinical.portion <- GDCprepare_clinic(query, clinical.info = "portion")

## ---- eval = FALSE-------------------------------------------------------
#  # This code will get all clinical indexed data from TCGA
#  library(TCGAbiolinks)
#  library(data.table)
#  clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>%
#              regexPipes::grep("TCGA",value=T) %>%
#              sort %>%
#              plyr::alply(1,GDCquery_clinic, .progress = "text") %>%
#              rbindlist
#  readr::write_csv(clinical,path = paste0("all_clin_indexed.csv"))
#  
#  # This code will get all clinical XML data from TCGA
#  getclinical <- function(proj){
#      message(proj)
#      while(1){
#          result = tryCatch({
#              query <- GDCquery(project = proj, data.category = "Clinical")
#              GDCdownload(query)
#              clinical <- GDCprepare_clinic(query, clinical.info = "patient")
#              for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
#                  message(i)
#                  aux <- GDCprepare_clinic(query, clinical.info = i)
#                  if(is.null(aux)) next
#                  # add suffix manually if it already exists
#                  replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
#                  colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
#                  if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
#              }
#              readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
#              return(clinical)
#          }, error = function(e) {
#              message(paste0("Error clinical: ", proj))
#          })
#      }
#  }
#  clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>%
#      regexPipes::grep("TCGA",value=T) %>% sort %>%
#      plyr::alply(1,getclinical, .progress = "text") %>%
#      rbindlist(fill = TRUE) %>% setDF %>% subset(!duplicated(clinical))
#  readr::write_csv(clinical,path = "all_clin_XML.csv")
#  # result: https://drive.google.com/open?id=0B0-8N2fjttG-WWxSVE5MSGpva1U
#  # Obs: this table has multiple lines for each patient, as the patient might have several followups, drug treatments,
#  # new tumor events etc...

## ---- eval = TRUE--------------------------------------------------------
bar <- c("TCGA-G9-6378-02A-11R-1789-07", "TCGA-CH-5767-04A-11R-1789-07",  
         "TCGA-G9-6332-60A-11R-1789-07", "TCGA-G9-6336-01A-11R-1789-07",
         "TCGA-G9-6336-11A-11R-1789-07", "TCGA-G9-7336-11A-11R-1789-07",
         "TCGA-G9-7336-04A-11R-1789-07", "TCGA-G9-7336-14A-11R-1789-07",
         "TCGA-G9-7036-04A-11R-1789-07", "TCGA-G9-7036-02A-11R-1789-07",
         "TCGA-G9-7036-11A-11R-1789-07", "TCGA-G9-7036-03A-11R-1789-07",
         "TCGA-G9-7036-10A-11R-1789-07", "TCGA-BH-A1ES-10A-11R-1789-07",
         "TCGA-BH-A1F0-10A-11R-1789-07", "TCGA-BH-A0BZ-02A-11R-1789-07",
         "TCGA-B6-A0WY-04A-11R-1789-07", "TCGA-BH-A1FG-04A-11R-1789-08",
         "TCGA-D8-A1JS-04A-11R-2089-08", "TCGA-AN-A0FN-11A-11R-8789-08",
         "TCGA-AR-A2LQ-12A-11R-8799-08", "TCGA-AR-A2LH-03A-11R-1789-07",
         "TCGA-BH-A1F8-04A-11R-5789-07", "TCGA-AR-A24T-04A-55R-1789-07",
         "TCGA-AO-A0J5-05A-11R-1789-07", "TCGA-BH-A0B4-11A-12R-1789-07",
         "TCGA-B6-A1KN-60A-13R-1789-07", "TCGA-AO-A0J5-01A-11R-1789-07",
         "TCGA-AO-A0J5-01A-11R-1789-07", "TCGA-G9-6336-11A-11R-1789-07",
         "TCGA-G9-6380-11A-11R-1789-07", "TCGA-G9-6380-01A-11R-1789-07",
         "TCGA-G9-6340-01A-11R-1789-07", "TCGA-G9-6340-11A-11R-1789-07")

S <- TCGAquery_SampleTypes(bar,"TP")
S2 <- TCGAquery_SampleTypes(bar,"NB")

# Retrieve multiple tissue types  NOT FROM THE SAME PATIENTS
SS <- TCGAquery_SampleTypes(bar,c("TP","NB"))

# Retrieve multiple tissue types  FROM THE SAME PATIENTS
SSS <- TCGAquery_MatchedCoupledSampleTypes(bar,c("NT","TP"))

## ---- eval = FALSE-------------------------------------------------------
#  # Check with subtypes from TCGAprepare and update examples
#  GBM_path_subtypes <- TCGAquery_subtype(tumor = "gbm")
#  
#  LGG_path_subtypes <- TCGAquery_subtype(tumor = "lgg")

## ---- eval = TRUE, echo = FALSE------------------------------------------
datatable(lgg.gbm.subtype[1:10,],
          caption = "Table with LGG molecular subtypes from TCGAquery_subtype",
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  query <- GDCquery(project = "TCGA-ACC", data.category = "Copy Number Variation",
#                    data.type = "Copy Number Segment",
#                    barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
#  GDCdownload(query)
#  data <- GDCprepare(query)
#  
#  #--------------------------------------
#  # Gene expression
#  #--------------------------------------
#  # Aligned against Hg38
#  # mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
#  query.exp.hg38 <- GDCquery(project = "TCGA-GBM",
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - FPKM-UQ",
#                    barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
#  GDCdownload(query.exp.hg38)
#  expdat <- GDCprepare(query = query.exp.hg38,
#                       save = TRUE,
#                       save.filename = "exp.rda")
#  
#  # Aligned against Hg19
#  query.exp.hg19 <- GDCquery(project = "TCGA-GBM",
#                    data.category = "Gene expression",
#                    data.type = "Gene expression quantification",
#                    platform = "Illumina HiSeq",
#                    file.type  = "normalized_results",
#                    experimental.strategy = "RNA-Seq",
#                    barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
#                    legacy = TRUE)
#  GDCdownload(query.exp.hg19)
#  data <- GDCprepare(query.exp.hg19)
#  
#  
#  #--------------------------------------
#  # DNA methylation data
#  #--------------------------------------
#  # DNA methylation aligned to hg38
#  query_met.hg38 <- GDCquery(project= "TCGA-LGG",
#                             data.category = "DNA Methylation",
#                             platform = "Illumina Human Methylation 450",
#                             barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"))
#  GDCdownload(query_met.hg38)
#  data.hg38 <- GDCprepare(query_met.hg38)
#  
#  # DNA methylation aligned to hg19
#  query_meth.hg19 <- GDCquery(project= "TCGA-LGG",
#                              data.category = "DNA methylation",
#                              platform = "Illumina Human Methylation 450",
#                              barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"),
#                              legacy = TRUE)
#  GDCdownload(query_meth.hg19)
#  data.hg19 <- GDCprepare(query_meth.hg19)
#  
#  # A function to download only 20 samples
#  legacyPipeline <- function(project, data.category, platform, n = 20){
#      query <- GDCquery(project = project,
#                        data.category = data.category,
#                        platform = platform,
#                        legacy = TRUE)
#      cases <- getResults(query, rows = 1:n, cols="cases") # Get two samples from the search
#      query <- GDCquery(project = project,
#                        data.category = data.category,
#                        platform = platform,
#                        legacy = TRUE,
#                        barcode = cases)
#      GDCdownload(query,chunks.per.download = 5)
#      data <- GDCprepare(query)
#      return(data)
#  }
#  # DNA methylation
#  data <- legacyPipeline("TCGA-GBM","DNA methylation","Illumina Human Methylation 27")
#  data <- legacyPipeline("TCGA-GBM","DNA methylation","Illumina Human Methylation 450")
#  data <- legacyPipeline("TCGA-GBM","DNA methylation","Illumina DNA Methylation OMA003 CPI")
#  data <- legacyPipeline("TCGA-GBM","DNA methylation","Illumina DNA Methylation OMA002 CPI")
#  
#  #-------------------------------------------------------
#  # Example to  idat files from TCGA projects
#  #-------------------------------------------------------
#  projects <- TCGAbiolinks:::getGDCprojects()$project_id
#  projects <- projects[grepl('^TCGA',projects,perl=T)]
#  match.file.cases.all <- NULL
#  for(proj in projects){
#      print(proj)
#      query <- GDCquery(project = proj,
#                        data.category = "Raw microarray data",
#                        data.type = "Raw intensities",
#                        experimental.strategy = "Methylation array",
#                        legacy = TRUE,
#                        file.type = ".idat",
#                        platform = "Illumina Human Methylation 450")
#      match.file.cases <- getResults(query,cols=c("cases","file_name"))
#      match.file.cases$project <- proj
#      match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
#      tryCatch(GDCdownload(query, method = "api",chunks.per.download = 20),
#               error = function(e) GDCdownload(query, method = "client"))
#  }
#  # This will create a map between idat file name, cases (barcode) and project
#  readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
#  # code to move all files to local folder
#  for(file in dir(".",pattern = ".idat", recursive = T)){
#      TCGAbiolinks::move(file,basename(file))
#  }

## ---- eval = FALSE-------------------------------------------------------
#  library(TCGAbiolinks)
#  # Downloading and prepare
#  query <- GDCquery(project = "TARGET-AML",
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - FPKM-UQ")
#  GDCdownload(query)
#  data <- GDCprepare(query)
#  
#  # Downloading and prepare using legacy
#  query <- GDCquery(project = "TCGA-GBM",
#                    data.category = "Protein expression",
#                    legacy = TRUE,
#                    barcode = c("TCGA-OX-A56R-01A-21-A44T-20","TCGA-08-0357-01A-21-1898-20"))
#  GDCdownload(query)
#  data <- GDCprepare(query, save = TRUE,
#                     save.filename = "gbmProteinExpression.rda",
#                     remove.files.prepared = TRUE)
#  
#  # Downloading and prepare using legacy
#  query <- GDCquery(project = "TCGA-GBM",
#                    data.category = "DNA methylation",
#                    platform = "Illumina Human Methylation 27",legacy = TRUE,
#                    barcode = c("TCGA-02-0047-01A-01D-0186-05","TCGA-06-2559-01A-01D-0788-05"))
#  GDCdownload(query)
#  data <- GDCprepare(query, add.gistic2.mut = c("PTEN","FOXJ1"))
#  
#  # To view gistic and mutation information please access the samples information matrix in the summarized Experiment object
#  library(SummarizedExperiment)
#  samples.information <- colData(data)
#  

## ---- eval = FALSE-------------------------------------------------------
#  # You can define a list of samples to query and download providing relative TCGA barcodes.
#  listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
#                   "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
#                   "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
#                   "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
#                   "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")
#  
#  # Query platform Illumina HiSeq with a list of barcode
#  query <- GDCquery(project = "TCGA-BRCA",
#                    data.category = "Gene expression",
#                    data.type = "Gene expression quantification",
#                    experimental.strategy = "RNA-Seq",
#                    platform = "Illumina HiSeq",
#                    file.type = "results",
#                    barcode = listSamples,
#                    legacy = TRUE)
#  
#  # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
#  GDCdownload(query)
#  
#  # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
#  # rsem.genes.results as values
#  BRCARnaseqSE <- GDCprepare(query)
#  
#  BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
#  
#  # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
#  BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)

## ---- eval = TRUE, echo = FALSE,size = 8---------------------------------
library(TCGAbiolinks)
dataGE <- dataBRCA[sample(rownames(dataBRCA),10),sample(colnames(dataBRCA),7)]

knitr::kable(dataGE[1:10,2:3], digits = 2, 
             caption = "Example of a matrix of gene expression (10 genes in rows and 2 samples in columns)",
             row.names = TRUE)

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
library(png)
library(grid)
img <- readPNG("PreprocessingOutput.png")
grid.raster(img)

## ---- eval = FALSE-------------------------------------------------------
#  # Downstream analysis using gene expression data
#  # TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
#  # save(dataBRCA, geneInfo , file = "dataGeneExpression.rda")
#  library(TCGAbiolinks)
#  
#  # normalization of genes
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)
#  
#  # quantile filter of genes
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    method = "quantile",
#                                    qnt.cut =  0.25)
#  
#  # selection of normal samples "NT"
#  samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
#                                     typesample = c("NT"))
#  
#  # selection of tumor samples "TP"
#  samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
#                                     typesample = c("TP"))
#  
#  # Diff.expr.analysis (DEA)
#  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
#                              mat2 = dataFilt[,samplesTP],
#                              Cond1type = "Normal",
#                              Cond2type = "Tumor",
#                              fdr.cut = 0.01 ,
#                              logFC.cut = 1,
#                              method = "glmLRT")
#  
#  # DEGs table with expression values in normal and tumor samples
#  dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
#                                            dataFilt[,samplesTP],dataFilt[,samplesNT])
#  

## ---- eval = TRUE, echo = FALSE------------------------------------------
library(TCGAbiolinks)
dataDEGsFiltLevel$FDR <- format(dataDEGsFiltLevel$FDR, scientific = TRUE)
knitr::kable(dataDEGsFiltLevel[1:10,], digits = 2,
             caption = "Table of DEGs after DEA", row.names = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  library(TCGAbiolinks)
#  # Enrichment Analysis EA
#  # Gene Ontology (GO) and Pathway enrichment by DEGs list
#  Genelist <- rownames(dataDEGsFiltLevel)
#  
#  system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
#  
#  # Enrichment Analysis EA (TCGAVisualize)
#  # Gene Ontology (GO) and Pathway enrichment barPlot
#  
#  TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
#                          GOBPTab = ansEA$ResBP,
#                          GOCCTab = ansEA$ResCC,
#                          GOMFTab = ansEA$ResMF,
#                          PathTab = ansEA$ResPat,
#                          nRGTab = Genelist,
#                          nBar = 10)
#  

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
library(png)
library(grid)
img <- readPNG("EAplot.png")
grid.raster(img)

## ---- eval = FALSE-------------------------------------------------------
#  clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
#  TCGAanalyze_survival(clin.gbm,
#                       "gender",
#                       main = "TCGA Set\n GBM",height = 10, width=10)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case2_surv.png")
grid.raster(img)

## ---- eval = FALSE-------------------------------------------------------
#  library(TCGAbiolinks)
#  # Survival Analysis SA
#  
#  clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
#  dataBRCAcomplete <- log2(BRCA_rnaseqv2)
#  
#  tokenStop<- 1
#  
#  tabSurvKMcomplete <- NULL
#  
#  for( i in 1: round(nrow(dataBRCAcomplete)/100)){
#      message( paste( i, "of ", round(nrow(dataBRCAcomplete)/100)))
#      tokenStart <- tokenStop
#      tokenStop <-100*i
#      tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
#                                        dataBRCAcomplete,
#                                        Genelist = rownames(dataBRCAcomplete)[tokenStart:tokenStop],
#                                        Survresult = F,
#                                        ThreshTop=0.67,
#                                        ThreshDown=0.33)
#  
#      tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
#  }
#  
#  tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
#  tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]
#  
#  tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
#      rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
#      ]

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
tabSurvKMcompleteDEGs$pvalue <- format(tabSurvKMcompleteDEGs$pvalue, scientific = TRUE)
knitr::kable(tabSurvKMcompleteDEGs[1:5,1:4], 
             digits = 2,
             caption = "Table KM-survival genes after SA",
             row.names = TRUE)
knitr::kable(tabSurvKMcompleteDEGs[1:5,5:7], 
             digits = 2,
             row.names = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  data <- TCGAanalyze_DMR(data, groupCol = "methylation_subtype",
#                          group1 = "CIMP.H",
#                          group2="CIMP.L",
#                          p.cut = 10^-5,
#                          diffmean.cut = 0.25,
#                          legend = "State",
#                          plot.filename = "coad_CIMPHvsCIMPL_metvolcano.png")

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5met.png")
grid.raster(img)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case2_Heatmap.png")
grid.raster(img)

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
#  mut <- GDCquery_Maf(tumor = "ACC", pipelines = "muse")
#  clin <- GDCquery_clinic("TCGA-ACC","clinical")
#  clin <- clin[,c("bcr_patient_barcode","disease","gender","tumor_stage","race","vital_status")]
#  TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:20],
#                          filename = "oncoprint.pdf",
#                          annotation = clin,
#                          color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
#                          rows.font.size=10,
#                          heatmap.legend.side = "right",
#                          dist.col = 0,
#                          width = 5,
#                          label.font.size = 10)

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
library(png)
library(grid)
img <- readPNG("oncoprint.png")
grid.raster(img)

## ---- eval = FALSE-------------------------------------------------------
#  # normalization of genes
#  dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#  
#  # quantile filter of genes
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    method = "quantile",
#                                    qnt.cut =  0.25)
#  
#  # selection of normal samples "NT"
#  group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#  # selection of normal samples "TP"
#  group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#  
#  # Principal Component Analysis plot for ntop selected DEGs
#   pca <- TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel, ntopgenes = 200, group1, group2)

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
library(png)
library(grid)
img <- readPNG("PCAtop200DEGs.png")
grid.raster(img)

## ---- eval = FALSE-------------------------------------------------------
#  library(TCGAbiolinks)
#  # Survival Analysis SA
#  
#  clinical_patient_Cancer <- TCGAquery_clinic("brca","clinical_patient")
#  dataBRCAcomplete <- log2(BRCA_rnaseqv2)
#  
#  tokenStop<- 1
#  
#  tabSurvKMcomplete <- NULL
#  
#  for( i in 1: round(nrow(dataBRCAcomplete)/100)){
#      message( paste( i, "of ", round(nrow(dataBRCAcomplete)/100)))
#      tokenStart <- tokenStop
#      tokenStop <-100*i
#      tabSurvKM <- TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
#                                          dataBRCAcomplete,
#                                          Genelist = rownames(dataBRCAcomplete)[tokenStart:tokenStop],
#                                          Survresult = F,ThreshTop=0.67,ThreshDown=0.33)
#      tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
#  }
#  
#  tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
#  tabSurvKMcomplete <- tabSurvKMcomplete[!duplicated(tabSurvKMcomplete$mRNA),]
#  rownames(tabSurvKMcomplete) <-tabSurvKMcomplete$mRNA
#  tabSurvKMcomplete <- tabSurvKMcomplete[,-1]
#  tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]
#  
#  tabSurvKMcompleteDEGs <- tabSurvKMcomplete[rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,]
#  
#  tflist <- EAGenes[EAGenes$Family == "transcription regulator","Gene"]
#  tabSurvKMcomplete_onlyTF <- tabSurvKMcomplete[rownames(tabSurvKMcomplete) %in% tflist,]
#  
#  TabCoxNet <- TCGAvisualize_SurvivalCoxNET(clinical_patient_Cancer,dataBRCAcomplete,
#                                            Genelist = rownames(tabSurvKMcompleteDEGs),
#                                            scoreConfidence = 700,
#                                            titlePlot = "TCGAvisualize_SurvivalCoxNET Example")

## ---- fig.width=6, fig.height=4, echo=FALSE, fig.align="center"----------
library(png)
library(grid)
img <- readPNG("SurvivalCoxNETOutput.png")
grid.raster(img)

## ----include=FALSE,echo=FALSE, fig.height=5, message=FALSE, warning=FALSE----
data <- tryCatch({
    query <- GDCquery(project = "TCGA-GBM",
                  data.category = "DNA methylation",
                  platform = "Illumina Human Methylation 27",
                  legacy = TRUE, 
                  barcode = c("TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
                              "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
                              "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
                              "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"))
            GDCdownload(query, method = "api", chunks.per.download = 2)
            GDCdownload(query, method = "api")
            data <- GDCprepare(query)
            data
    }, error = function(e) {
         nrows <- 200; ncols <- 21
         counts <- matrix(runif(nrows * ncols, 0, 1), nrows)
         rowRanges <- GenomicRanges::GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                            IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                            strand=sample(c("+", "-"), 200, TRUE),
                            feature_id=sprintf("ID%03d", 1:200))
        colData <- S4Vectors::DataFrame(shortLetterCode=rep(c("NT", "TP","TR"), 7),
                         row.names=LETTERS[1:21],
                         subtype_Pan.Glioma.DNA.Methylation.Cluster=rep(c("group1","group2","group3"),c(7,7,7)),
                         vital_status=rep(c("DEAD","ALIVE","DEAD"),7))
        data <- SummarizedExperiment::SummarizedExperiment(
                    assays=S4Vectors::SimpleList(counts=counts),
                    rowRanges=rowRanges,
                    colData=colData)
        data

        }
    )

## ---- eval=FALSE, echo=TRUE, fig.height=5, message=FALSE, warning=FALSE----
#  query <- GDCquery(project = "TCGA-GBM",
#                    data.category = "DNA methylation",
#                    platform = "Illumina Human Methylation 27",
#                    legacy = TRUE,
#                    barcode = c("TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
#                                "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
#                                "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
#                                "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"))
#  GDCdownload(query, method = "api")
#  data <- GDCprepare(query)

## ---- echo=TRUE, fig.height=4, fig.width=3, out.width = 3, out.heigh=5, message=FALSE, warning=FALSE----
# "shortLetterCode" is a column in the SummarizedExperiment::colData(data) matrix
TCGAvisualize_meanMethylation(data, groupCol = "shortLetterCode",filename = NULL)

## ---- echo=TRUE, fig.height=4, fig.width=7, out.width = 7, out.heigh=5, message=FALSE, warning=FALSE----
# setting y limits: lower 0, upper 1
TCGAvisualize_meanMethylation(data,groupCol = "shortLetterCode", 
                              filename = NULL, y.limits = c(0,1))
# setting y limits: lower 0
TCGAvisualize_meanMethylation(data,groupCol = "shortLetterCode", 
                              filename = NULL, y.limits = 0)

# Changing shapes of jitters to show subgroups
TCGAvisualize_meanMethylation(data,
                              groupCol = "subtype_Pan.Glioma.DNA.Methylation.Cluster", 
                              subgroupCol ="vital_status", filename = NULL)

# Sorting bars by descending mean methylation
TCGAvisualize_meanMethylation(data,
                              groupCol = "subtype_Pan.Glioma.DNA.Methylation.Cluster",
                              sort="mean.desc",
                              filename=NULL)
# Sorting bars by asc mean methylation
TCGAvisualize_meanMethylation(data,
                              groupCol = "subtype_Pan.Glioma.DNA.Methylation.Cluster",
                              sort = "mean.asc",
                              filename=NULL)

TCGAvisualize_meanMethylation(data,
                              groupCol = "vital_status",
                              sort = "mean.asc",
                              filename=NULL)

## ---- eval = FALSE-------------------------------------------------------
#  starburst <- TCGAvisualize_starburst(coad.SummarizeExperiment,
#                                       different.experssion.analysis.data,
#                                       group1 = "CIMP.H",
#                                       group2 = "CIMP.L",
#                                       met.p.cut = 10^-5,
#                                       exp.p.cut=10^-5,
#                                       names = TRUE)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5star.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  CancerProject <- "TCGA-BRCA"
#  DataDirectory <- paste0("../GDC/",gsub("-","_",CancerProject))
#  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")
#  
#  query <- GDCquery(project = CancerProject,
#                    data.category = "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "HTSeq - Counts")
#  
#  samplesDown <- getResults(query,cols=c("cases"))
#  
#  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                    typesample = "TP")
#  
#  dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                    typesample = "NT")
#  dataSmTP_short <- dataSmTP[1:10]
#  dataSmNT_short <- dataSmNT[1:10]
#  
#  queryDown <- GDCquery(project = CancerProject,
#                        data.category = "Transcriptome Profiling",
#                        data.type = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - Counts",
#                        barcode = c(dataSmTP_short, dataSmNT_short))
#  
#  GDCdownload(query = queryDown,
#              directory = DataDirectory)
#  
#  dataPrep <- GDCprepare(query = queryDown,
#                         save = TRUE,
#                         directory =  DataDirectory,
#                         save.filename = FileNameData)
#  
#  dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep,
#                                        cor.cut = 0.6,
#                                        datatype = "HTSeq - Counts")
#  
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                        geneInfo = geneInfoHT,
#                                        method = "gcContent")
#  
#  boxplot(dataPrep, outline = FALSE)
#  
#  boxplot(dataNorm, outline = FALSE)
#  
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    method = "quantile",
#                                    qnt.cut =  0.25)
#  
#  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmTP_short],
#                              mat2 = dataFilt[,dataSmNT_short],
#                              Cond1type = "Normal",
#                              Cond2type = "Tumor",
#                              fdr.cut = 0.01 ,
#                              logFC.cut = 1,
#                              method = "glmLRT")
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  require(TCGAbiolinks)
#  
#  CancerProject <- "TCGA-BRCA"
#  DataDirectory <- paste0("../GDC/",gsub("-","_",CancerProject))
#  FileNameData <- paste0(DataDirectory, "_","miRNA_gene_quantification",".rda")
#  
#  query.miR <- GDCquery(project = CancerProject,
#                    data.category = "Gene expression",
#                    data.type = "miRNA gene quantification",
#                    file.type = "hg19.mirna",
#                    legacy = TRUE)
#  
#  samplesDown.miR <- getResults(query.miR,cols=c("cases"))
#  
#  dataSmTP.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR,
#                                    typesample = "TP")
#  
#  dataSmNT.miR <- TCGAquery_SampleTypes(barcode = samplesDown.miR,
#                                    typesample = "NT")
#  dataSmTP_short.miR <- dataSmTP.miR[1:10]
#  dataSmNT_short.miR <- dataSmNT.miR[1:10]
#  
#  queryDown.miR <- GDCquery(project = CancerProject,
#                        data.category = "Gene expression",
#                        data.type = "miRNA gene quantification",
#                        file.type = "hg19.mirna",
#                        legacy = TRUE,
#                        barcode = c(dataSmTP_short.miR, dataSmNT_short.miR))
#  
#  GDCdownload(query = queryDown.miR,
#              directory = DataDirectory)
#  
#  dataAssy.miR <- GDCprepare(query = queryDown.miR,
#                             save = TRUE,
#                             save.filename = FileNameData,
#                             summarizedExperiment = TRUE,
#                             directory =DataDirectory )
#  rownames(dataAssy.miR) <- dataAssy.miR$miRNA_ID
#  
#  # using read_count's data
#  read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))]
#  dataAssy.miR <- dataAssy.miR[,read_countData]
#  colnames(dataAssy.miR) <- gsub("read_count_","", colnames(dataAssy.miR))
#  
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataAssy.miR,
#                                    method = "quantile",
#                                    qnt.cut =  0.25)
#  
#  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT_short.miR],
#                              mat2 = dataFilt[,dataSmTP_short.miR],
#                              Cond1type = "Normal",
#                              Cond2type = "Tumor",
#                              fdr.cut = 0.01 ,
#                              logFC.cut = 1,
#                              method = "glmLRT")
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  library(SummarizedExperiment)
#  library(TCGAbiolinks)
#  query.exp <- GDCquery(project = "TCGA-BRCA",
#                        legacy = TRUE,
#                        data.category = "Gene expression",
#                        data.type = "Gene expression quantification",
#                        platform = "Illumina HiSeq",
#                        file.type = "results",
#                        experimental.strategy = "RNA-Seq",
#                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
#  GDCdownload(query.exp)
#  brca.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "brcaExp.rda")
#  
#  # get subtype information
#  dataSubt <- TCGAquery_subtype(tumor = "BRCA")
#  
#  # get clinical data
#  dataClin <- GDCquery_clinic(project = "TCGA-BRCA","clinical")
#  
#  # Which samples are primary solid tumor
#  dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP")
#  # which samples are solid tissue normal
#  dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  dataPrep <- TCGAanalyze_Preprocessing(object = brca.exp, cor.cut = 0.6)
#  
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                        geneInfo = geneInfo,
#                                        method = "gcContent")
#  
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    method = "quantile",
#                                    qnt.cut =  0.25)
#  
#  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
#                              mat2 = dataFilt[,dataSmTP],
#                              Cond1type = "Normal",
#                              Cond2type = "Tumor",
#                              fdr.cut = 0.01 ,
#                              logFC.cut = 1,
#                              method = "glmLRT")
#  
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",
#                                  RegulonList = rownames(dataDEGs))
#  
#  TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
#                          GOBPTab = ansEA$ResBP,
#                          GOCCTab = ansEA$ResCC,
#                          GOMFTab = ansEA$ResMF,
#                          PathTab = ansEA$ResPat,
#                          nRGTab = rownames(dataDEGs),
#                          nBar = 20)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case1_EA.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#  group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#  
#  dataSurv <- TCGAanalyze_SurvivalKM(clinical_patient = dataClin,
#                                     dataGE = dataFilt,
#                                     Genelist = rownames(dataDEGs),
#                                     Survresult = FALSE,
#                                     ThreshTop = 0.67,
#                                     ThreshDown = 0.33,
#                                     p.cut = 0.05, group1, group2)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  require(dnet)  # to change
#  org.Hs.string <- dRDataLoader(RData = "org.Hs.string")
#  
#  TabCoxNet <- TCGAvisualize_SurvivalCoxNET(dataClin,
#                                            dataFilt,
#                                            Genelist = rownames(dataSurv),
#                                            scoreConfidence = 700,
#                                            org.Hs.string = org.Hs.string,
#                                            titlePlot = "Case Study n.1 dnet")

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case1_dnet.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  
#  query.exp <- GDCquery(project = "TCGA-LGG",
#                        legacy = TRUE,
#                        data.category = "Gene expression",
#                        data.type = "Gene expression quantification",
#                        platform = "Illumina HiSeq",
#                        file.type = "results",
#                        experimental.strategy = "RNA-Seq",
#                        sample.type = "Primary solid Tumor")
#  GDCdownload(query.exp)
#  lgg.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "lggExp.rda")
#  
#  # get subtype information
#  dataSubt <- TCGAquery_subtype(tumor = "LGG")
#  
#  # get indexed clinical data
#  dataClin <- GDCquery_clinic(project = "TCGA-LGG", "Clinical")
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  # expression data with molecular subtypes
#  lgg.exp <- subset(lgg.exp, select = colData(lgg.exp)$patient %in% dataSubt$patient)
#  
#  dataPrep <- TCGAanalyze_Preprocessing(object = lgg.exp,cor.cut = 0.6)
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                        geneInfo = geneInfo,
#                                        method = "gcContent")
#  
#  datFilt1 <- TCGAanalyze_Filtering(tabDF = dataNorm,method = "varFilter")
#  datFilt2 <- TCGAanalyze_Filtering(tabDF = datFilt1,method = "filter1")
#  datFilt <- TCGAanalyze_Filtering(tabDF = datFilt2,method = "filter2")
#  
#  data_Hc1 <- TCGAanalyze_Clustering(tabDF = datFilt,
#                                     method = "hclust",
#                                     methodHC = "ward.D2")
#  data_Hc2 <- TCGAanalyze_Clustering(tabDF = datFilt,
#                                     method = "consensus",
#                                     methodHC = "ward.D2")

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  #------  Add cluster information
#  cluster <- data.frame("groupsHC" = data_Hc2[[4]]$consensusClass)
#  cluster$groupsHC <- paste0("EC",cluster$groupsHC)
#  cluster$patient <-  substr(colData(lgg.exp)$patient,1,12)
#  
#  # Add information about gropus from consensus Cluster in clinical data
#  dataClin <- merge(dataClin,cluster, by.x="bcr_patient_barcode", by.y="patient")
#  
#  # Merge subtype and clinical data
#  clin_subt <- merge(dataClin,dataSubt, by.x="bcr_patient_barcode", by.y="patient")
#  clin_subt_all <- merge(dataClin,dataSubt,
#                         by.x="bcr_patient_barcode", by.y="patient", all.x = TRUE)
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  #----------- VISUALIZE --------------------
#  # plotting survival for groups EC1, EC2, EC3, EC4
#  TCGAanalyze_survival(data = clin_subt_all,
#                       clusterCol = "groupsHC",
#                       main = "TCGA kaplan meier survival plot from consensus cluster",
#                       legend = "RNA Group",
#                       color = c("black","red","blue","green3"),
#                       filename = "case2_surv.png")

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case2_surv.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  TCGAvisualize_Heatmap(t(datFilt),
#                        col.metadata =  clin_subt[,c("bcr_patient_barcode",
#                                                     "groupsHC",
#                                                     "Histology",
#                                                     "IDH.codel.subtype")],
#                        col.colors =  list(
#                            groupsHC = c("EC1"="black",
#                                         "EC2"="red",
#                                         "EC3"="blue",
#                                         "EC4"="green3"),
#                            Histology=c("astrocytoma"="navy",
#                                        "oligoastrocytoma"="green3",
#                                        "oligodendroglioma"="red"),
#                            IDH.codel.subtype = c("IDHmut-codel"="tomato",
#                                                  "IDHmut-non-codel"="navy",
#                                                  "IDHwt"="gold","NA"="white")),
#                        sortCol = "groupsHC",
#                        type = "expression", # sets default color
#                        scale = "row", # use z-scores for better visualization
#                        title = "Heatmap from concensus cluster",
#                        filename = "case2_Heatmap.pdf",
#                        cluster_rows = TRUE)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case2_Heatmap.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "muse")
#  # Selecting gene
#  mRNAsel <- "ATRX"
#  LGGselected <- LGGmut[LGGmut$Hugo_Symbol == mRNAsel,]
#  
#  dataMut <- LGGselected[!duplicated(LGGselected$Tumor_Sample_Barcode),]
#  dataMut$Tumor_Sample_Barcode <- substr(dataMut$Tumor_Sample_Barcode,1,12)
#  
#  # Adding the Expression Cluster classification found before
#  dataMut <- merge(dataMut, cluster, by.y="patient", by.x="Tumor_Sample_Barcode")
#  dataMut <- dataMut[dataMut$Variant_Classification!=0,]

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  dir.create("case3")
#  setwd("case3")
#  #-----------------------------------
#  # STEP 1: Search, download, prepare |
#  #-----------------------------------
#  # 1.1 - DNA methylation
#  # ----------------------------------
#  query.met <- GDCquery(project = "TCGA-ACC",
#                        legacy = TRUE,
#                        data.category = "DNA methylation",
#                        platform = "Illumina Human Methylation 450")
#  GDCdownload(query.met)
#  
#  acc.met <- GDCprepare(query = query.met,
#                        save = TRUE,
#                        save.filename = "accDNAmet.rda",
#                        summarizedExperiment = TRUE)
#  
#  #-----------------------------------
#  # 1.2 - RNA expression
#  # ----------------------------------
#  query.exp <- GDCquery(project = "TCGA-ACC",
#                        legacy = TRUE,
#                        data.category = "Gene expression",
#                        data.type = "Gene expression quantification",
#                        platform = "Illumina HiSeq",
#                        file.type = "results")
#  GDCdownload(query.exp)
#  acc.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "accExp.rda")

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  #--------------------------------------------
#  # STEP 2: Analysis
#  #--------------------------------------------
#  # 2.1 - Mean methylation of samples
#  # -------------------------------------------
#  TCGAvisualize_meanMethylation(acc.met,
#                                groupCol = "subtype_MethyLevel",
#                                subgroupCol = "subtype_Histology",
#                                group.legend  = "Groups",
#                                subgroup.legend = "Histology",
#                                filename = "acc_mean.png")

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("acc_mean.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  # na.omit
#  acc.met <- subset(acc.met,subset = (rowSums(is.na(assay(acc.met))) == 0))
#  
#  # Volcano plot
#  acc.met <- TCGAanalyze_DMR(acc.met, groupCol = "subtype_MethyLevel",
#                             group1 = "CIMP-high",
#                             group2="CIMP-low",
#                             p.cut = 10^-5,
#                             diffmean.cut = 0.25,
#                             legend = "State",
#                             plot.filename = "CIMP-highvsCIMP-low_metvolcano.png")

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("CIMP-highvsCIMP-low_metvolcano.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  #-------------------------------------------------
#  # 2.3 - DEA - Expression analysis - volcano plot
#  # ------------------------------------------------
#  acc.exp.aux <- subset(acc.exp,
#                        select = colData(acc.exp)$subtype_MethyLevel %in% c("CIMP-high","CIMP-low"))
#  
#  idx <- colData(acc.exp.aux)$subtype_MethyLevel %in% c("CIMP-high")
#  idx2 <- colData(acc.exp.aux)$subtype_MethyLevel %in% c("CIMP-low")
#  
#  dataPrep <- TCGAanalyze_Preprocessing(object = acc.exp.aux, cor.cut = 0.6)
#  
#  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                        geneInfo = geneInfo,
#                                        method = "gcContent")
#  
#  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    qnt.cut = 0.25,
#                                    method='quantile')
#  
#  dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,idx],
#                              mat2 = dataFilt[,idx2],
#                              Cond1type = "CIMP-high",
#                              Cond2type = "CIMP-low",
#                              method = "glmLRT")
#  
#  TCGAVisualize_volcano(dataDEGs$logFC,dataDEGs$FDR,
#                        filename = "Case3_volcanoexp.png",
#                        x.cut = 3,
#                        y.cut = 10^-5,
#                        names = rownames(dataDEGs),
#                        color = c("black","red","darkgreen"),
#                        names.size = 2,
#                        xlab = " Gene expression fold change (Log2)",
#                        legend = "State",
#                        title = "Volcano plot (CIMP-high vs CIMP-low)",
#                        width = 10)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5exp.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  
#  #------------------------------------------
#  # 2.4 - Starburst plot
#  # -----------------------------------------
#  # If true the argument names of the geenes in circle
#  # (biologically significant genes, has a change in gene
#  # expression and DNA methylation and respects all the thresholds)
#  # will be shown
#  # these genes are returned by the function see starburst object after the function is executed
#  starburst <- TCGAvisualize_starburst(acc.met, dataDEGs,"CIMP-high","CIMP-low",
#                                       filename = "starburst.png",
#                                       met.p.cut = 10^-5,
#                                       exp.p.cut = 10^-5,
#                                       diffmean.cut = 0.25,
#                                       logFC.cut = 3,
#                                       names = FALSE,
#                                       height=10,
#                                       width=15,
#                                       dpi=300)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("figure5star.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  library(TCGAbiolinks)
#  library(SummarizedExperiment)
#  library(ELMER)
#  library(parallel)
#  dir.create("case4")
#  setwd("case4")
#  #-----------------------------------
#  # STEP 1: Search, download, prepare |
#  #-----------------------------------
#  # 1.1 - DNA methylation
#  # ----------------------------------
#  query.met <- GDCquery(project = "TCGA-KIRC",
#                        legacy = TRUE,
#                        data.category = "DNA methylation",
#                        platform = "Illumina Human Methylation 450")
#  GDCdownload(query.met)
#  kirc.met <- GDCprepare(query = query.met,
#                         save = TRUE,
#                         save.filename = "kircDNAmet.rda",
#                         summarizedExperiment = TRUE)
#  
#  
#  kirc.met <- TCGAprepare_elmer(kirc.met,
#                                platform = "HumanMethylation450",
#                                save = TRUE,
#                                met.na.cut = 0.2)
#  

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  # Step 1.2 download expression data
#  #-----------------------------------
#  # 1.2 - RNA expression
#  # ----------------------------------
#  query.exp <- GDCquery(project = "TCGA-KIRC",
#                        legacy = TRUE,
#                        data.category = "Gene expression",
#                        data.type = "Gene expression quantification",
#                        platform = "Illumina HiSeq",
#                        file.type = "normalized_results")
#  GDCdownload(query.exp)
#  kirc.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "kirkExp.rda")
#  
#  kirc.exp <- TCGAprepare_elmer(kirc.exp,
#                                save = TRUE,
#                                platform = "IlluminaHiSeq_RNASeqV2")

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  #-----------------------------------
#  # STEP 2: ELMER integration         |
#  #-----------------------------------
#  # Step 2.1 prepare mee object       |
#  # -----------------------------------
#  library(ELMER)
#  library(parallel)
#  
#  geneAnnot <- txs()
#  geneAnnot$GENEID <- paste0("ID",geneAnnot$GENEID)
#  geneInfo <- promoters(geneAnnot,upstream = 0, downstream = 0)
#  probe <- get.feature.probe()
#  mee <- fetch.mee(meth = kirc.met, exp = kirc.exp, TCGA = TRUE,
#                   probeInfo = probe, geneInfo = geneInfo)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  direction <- c("hyper","hypo")
#  
#  for (j in direction){
#      print(j)
#      dir.out <- paste0("kirc/",j)
#      dir.create(dir.out, recursive = TRUE)
#      #--------------------------------------
#      # STEP 3: Analysis                     |
#      #--------------------------------------
#      # Step 3.1: Get diff methylated probes |
#      #--------------------------------------
#      Sig.probes <- get.diff.meth(mee, cores=detectCores(),
#                                  dir.out =dir.out,
#                                  diff.dir=j,
#                                  pvalue = 0.01)
#  
#      #-------------------------------------------------------------
#      # Step 3.2: Identify significant probe-gene pairs            |
#      #-------------------------------------------------------------
#      # Collect nearby 20 genes for Sig.probes
#      nearGenes <- GetNearGenes(TRange=getProbeInfo(mee, probe=Sig.probes$probe),
#                                cores=detectCores(),
#                                geneAnnot=getGeneInfo(mee))
#  
#      pair <- get.pair(mee=mee,
#                       probes=Sig.probes$probe,
#                       nearGenes=nearGenes,
#                       permu.dir=paste0(dir.out,"/permu"),
#                       dir.out=dir.out,
#                       cores=detectCores(),
#                       label= j,
#                       permu.size=10000,
#                       Pe = 0.001)
#  
#      Sig.probes.paired <- fetch.pair(pair=pair,
#                                      probeInfo = getProbeInfo(mee),
#                                      geneInfo = getGeneInfo(mee))
#      Sig.probes.paired <-read.csv(paste0(dir.out,"/getPair.",j,".pairs.significant.csv"),
#                                   stringsAsFactors=FALSE)[,1]
#  
#  
#      #-------------------------------------------------------------
#      # Step 3.3: Motif enrichment analysis on the selected probes |
#      #-------------------------------------------------------------
#      if(length(Sig.probes.paired) > 0 ){
#          #-------------------------------------------------------------
#          # Step 3.3: Motif enrichment analysis on the selected probes |
#          #-------------------------------------------------------------
#          enriched.motif <- get.enriched.motif(probes=Sig.probes.paired,
#                                               dir.out=dir.out, label=j,
#                                               background.probes = probe$name)
#          if(length(enriched.motif) > 0){
#              #-------------------------------------------------------------
#              # Step 3.4: Identifying regulatory TFs                        |
#              #-------------------------------------------------------------
#              print("get.TFs")
#  
#              TF <- get.TFs(mee = mee,
#                            enriched.motif = enriched.motif,
#                            dir.out = dir.out,
#                            cores = detectCores(), label = j)
#              save(TF, enriched.motif, Sig.probes.paired,
#                   pair, nearGenes, Sig.probes,
#                   file=paste0(dir.out,"/ELMER_results_",j,".rda"))
#          }
#      }
#  }

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  scatter.plot(mee,byProbe=list(probe=c("cg00328720"),geneNum=20), category="TN", lm_line = TRUE)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("case4_elmer.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  scatter.plot(mee,byTF = list(TF = c("ZNF677","PEG3"),
#                            probe = enriched.motif[["UA6"]]), category = "TN",
#                            save = TRUE, lm_line = TRUE, dir.out = "kirc")

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("elmer1.png")
grid.raster(img)

## ----eval=FALSE,echo=TRUE,message=FALSE,warning=FALSE--------------------
#  pair <- fetch.pair(pair="./kirc/getPair.hypo.pairs.significant.csv",
#                     probeInfo = mee@probeInfo, geneInfo = mee@geneInfo)
#  schematic.plot(pair=pair, byProbe="cg15862394",save=FALSE)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("elmer2.png")
grid.raster(img)

## ---- fig.width=6, fig.height=4, echo = FALSE, fig.align="center",hide=TRUE, message=FALSE,warning=FALSE----
library(png)
library(grid)
img <- readPNG("elmer3.png")
grid.raster(img)

## ----sessionInfo---------------------------------------------------------
sessionInfo()

