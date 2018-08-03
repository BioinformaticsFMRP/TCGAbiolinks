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
    cols <- c("gender","project_id","days_to_last_follow_up","alcohol_history","cigarettes_per_day")
    x <- getBarcodeInfo(c("TCGA-OR-A5LR", "TCGA-OR-A5LJ"))
    expect_true(all(cols %in% colnames(x)))

    cols <- c("gender","project_id")
    x <- getBarcodeInfo(c("TARGET-20-PARUDL"))
    expect_true(all(cols %in% colnames(x)))

})

test_that("GDCprepare accepts more than one project", {
    cases <-  c("TCGA-OR-A5JX-01A-11D-A29H-01", "TCGA-OR-A5J3-01A-11D-A29H-01",
                "TCGA-02-0010-10A-01D-0182-01","TCGA-14-0871-01A-01D-0384-01")
    acc.gbm<- GDCquery(project =  c("TCGA-ACC","TCGA-GBM"),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM-UQ",
                        barcode = substr(cases,1,12))
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex")
    expect_true(all(c("TCGA-ACC","TCGA-GBM") %in% SummarizedExperiment::colData(obj)$project_id))
})
