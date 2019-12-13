context("Download AND PREPARE")

test_that("GDCdownload API method is working ", {
    cases <-  c("TCGA-OR-A5JX-01A-11R-A29S-07",
                "TCGA-OR-A5KY-01A-11R-A29S-07",
                "TCGA-PK-A5HA-01A-11R-A29S-07")
    acc.gbm <- GDCquery(project =  c("TCGA-ACC"),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM-UQ",
                        barcode = substr(cases,1,12))
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex",summarizedExperiment = FALSE)
    expect_true(all(substr(colnames(obj)[-1],1,12) == substr(cases,1,12)))

    obj <- GDCprepare(acc.gbm,  directory = "ex",summarizedExperiment = TRUE)
    expect_true(all(substr(colData(obj) %>% rownames(),1,12) == substr(cases,1,12)))
    expect_true(all(obj$barcode == cases))

    query <- GDCquery(
        project = "CPTAC-3",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "HTSeq - Counts",
        barcode = c("CPT0010260013","CPT0000870008","CPT0105190006","CPT0077490006")
    )
    GDCdownload(query)
    data <- GDCprepare(query)
    expect_true(all(query$results[[1]]$sample.submitter_id == colnames(data)))
    expect_true(all(query$results[[1]]$sample.submitter_id == data$sample_submitter_id))
})




test_that("getBarcodeInfo works", {
    cols <- c("gender","project_id","days_to_last_follow_up","alcohol_history")
    x <- getBarcodeInfo(c("TCGA-OR-A5LR-01A", "TCGA-OR-A5LJ-01A"))
    expect_true(all(cols %in% colnames(x)))

    cols <- c("gender","project_id")
    x <- getBarcodeInfo(c("TARGET-20-PARUDL-03A"))
    expect_true(all(cols %in% colnames(x)))

    samples <- c("HCM-CSHL-0063-C18-85A",
                 "HCM-CSHL-0065-C20-06A",
                 "HCM-CSHL-0065-C20-85A",
                 "HCM-CSHL-0063-C18-01A")
    x <- colDataPrepare(samples)
    expect_true(all(rownames(x) == samples))
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","gender"] == "male")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","tumor_grade"] == "G2")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","ajcc_pathologic_stage"] == "Stage IVA")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","sample_type"] == "Metastatic")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0063-C18-85A","sample_type"] == "Next Generation Cancer Model")

    x <- getBarcodeInfo(c("HCM-CSHL-0063-C18-85A",
                          "HCM-CSHL-0065-C20-06A",
                          "HCM-CSHL-0065-C20-85A",
                          "TARGET-20-PARUDL-03A",
                          "TCGA-OR-A5LR-01A",
                          "HCM-CSHL-0063-C18-01A"))
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","gender"] == "male")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","tumor_grade"] == "G2")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","ajcc_pathologic_stage"] == "Stage IVA")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","sample_type"] == "Metastatic")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0063-C18-85A","sample_type"] == "Next Generation Cancer Model")



})

test_that("colDataPrepare handle replicates", {
    barcodes <- c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01")
    x <- colDataPrepare(barcodes)
    expect_true(nrow(x) == 2)
    expect_true(all(rownames(x) == c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01")))
    expect_true(all(x$barcode == c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01")))
})

test_that("GDCprepare accepts more than one project", {
    cases <-  c("TCGA-OR-A5JX-01A", "TCGA-OR-A5J3-01A",
                "TCGA-06-0680-11A","TCGA-14-0871-01A")
    expect_true(all(c("TCGA-ACC","TCGA-GBM") %in% colDataPrepare(cases)$project_id))
    acc.gbm <- GDCquery(project =  c("TCGA-ACC","TCGA-GBM"),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM-UQ",
                        barcode = substr(cases,1,12))
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex")
    expect_true(all(c("TCGA-ACC","TCGA-GBM") %in% SummarizedExperiment::colData(obj)$project_id))
})

test_that("Non TCGA data is processed", {
    proj <- "MMRF-COMMPASS"
    query <- GDCquery(
        project = proj,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
    )
    query <- GDCquery(
        project = proj,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        barcode = getResults(query)$cases[1:4]
    )
    GDCdownload(query)
    data <- GDCprepare(query)
})

test_that("GISTIC2 data is being correclty prepare", {
    query <- GDCquery(project = "TCGA-COAD",
                      data.category = "Copy Number Variation",
                      data.type = "Gene Level Copy Number Scores",
                      access = "open")
    GDCdownload(query,directory = "ex")
    data <- GDCprepare(query,directory = "ex")

    files <- dir("ex",pattern = "focal_score",recursive = TRUE,full.names = TRUE)
    raw.data <- readr::read_tsv(files)
    idx <- match(c( "79a12e57-0154-4de3-a6a4-80b6323b7cb3",
                    "cfd4127e-cd08-4f8c-b5b2-e440b452e044"),
                 colnames(raw.data))
    expect_true(all(substr(colnames(data)[idx],1,12) == c("TCGA-A6-5664","TCGA-AY-A71X")))
    unlink("ex",recursive = TRUE,force = TRUE)
})

test_that("IDAT files is processed", {
proj <- "TCGA-LUAD"
query <- GDCquery(project = proj,
                  data.category = "Raw microarray data",
                  data.type = "Raw intensities",
                  experimental.strategy = "Methylation array",
                  legacy = TRUE,
                  file.type = ".idat",
                  barcode = "TCGA-55-7724",
                  platform = "Illumina Human Methylation 450")

    tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
         error = function(e) GDCdownload(query, method = "client"))
    betas <- GDCprepare(query)
    expect_true(nrow(betas) == 485577)
    expect_true(ncol(betas) == 1)
})

test_that("Prepare Samples without clinical data", {
    # x <-  GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
    # x[is.na(x$diagnosis_id),]
    x <- colDataPrepare(c("TCGA-80-5608-01A","TCGA-17-Z053-01A","TCGA-78-7158-01A"))
   expect_true(nrow(x) == 3)
})
