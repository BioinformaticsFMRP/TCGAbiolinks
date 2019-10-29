context("Download AND PREPARE")

test_that("GDCdownload API method for one files is working ", {
    cases <-  c("TCGA-OR-A5JX-01A-11D-A29H-01")
    acc.gbm <- GDCquery(project =  c("TCGA-ACC"),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM-UQ",
                        barcode = substr(cases,1,12))
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex",summarizedExperiment = FALSE)
    expect_true( substr(colnames(obj)[2],1,12) == substr(cases,1,12))
})

test_that("getBarcodeInfo works", {
    cols <- c("gender","project_id","days_to_last_follow_up","alcohol_history")
    x <- getBarcodeInfo(c("TCGA-OR-A5LR-01A", "TCGA-OR-A5LJ-01A"))
    expect_true(all(cols %in% colnames(x)))

    cols <- c("gender","project_id")
    x <- getBarcodeInfo(c("TARGET-20-PARUDL-03A"))
    expect_true(all(cols %in% colnames(x)))

    x <- colDataPrepare(c("HCM-CSHL-0063-C18-85A",
                          "HCM-CSHL-0065-C20-06A",
                          "HCM-CSHL-0065-C20-85A",
                          "HCM-CSHL-0063-C18-01A"))
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
