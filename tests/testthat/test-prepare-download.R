context("Download AND PREPARE")

test_that("GDCdownload API method for two files is working ", {
    sink("/dev/null");
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy number variation",
                      legacy = TRUE,
                      file.type = "hg19.seg",
                      barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01", "TCGA-OR-A5LJ-10A-01D-A29K-01"))
    # data will be saved in  GDCdata/TCGA-ACC/legacy/Copy_number_variation/Copy_number_segmentation
    GDCdownload(query, method = "api")
    files <- file.path("GDCdata/TCGA-ACC/legacy/",
                       gsub(" ","_",query$results[[1]]$data_category),
                       gsub(" ","_",query$results[[1]]$data_type),
                       query$results[[1]]$file_id,
                       query$results[[1]]$file_name)
    expect_true(all(file.exists(files)))
    unlink("GDCdata",recursive = TRUE, force = TRUE)
})
test_that("GDCdownload API method for one files is working ", {
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy number variation",
                      legacy = TRUE,
                      file.type = "hg19.seg",
                      barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
    # data will be saved in  GDCdata/TCGA-ACC/legacy/Copy_number_variation/Copy_number_segmentation
    #GDCdownload(query, method = "api", directory = "example_data_dir")
    GDCdownload(query, method = "api", directory = "ex")
    files <- file.path("ex/TCGA-ACC/legacy/",
                       gsub(" ","_",query$results[[1]]$data_category),
                       gsub(" ","_",query$results[[1]]$data_type),
                       query$results[[1]]$file_id,
                       query$results[[1]]$file_name)
    expect_true(all(file.exists(files)))
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
    acc.gbm <- GDCquery(project =  c("TCGA-ACC","TCGA-GBM"),
                        data.category = "Copy Number Variation",
                        data.type = "Copy Number Segment",
                        barcode = cases)
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex")
    expect_true(all(cases %in% obj$Sample))
    acc.gbm<- GDCquery(project =  c("TCGA-ACC","TCGA-GBM"),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts",
                        barcode = substr(cases,1,12))
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex")
    expect_true(all(c("TCGA-ACC","TCGA-GBM") %in% SummarizedExperiment::colData(obj)$project_id))
})
test_that("Accecpts more than one platform", {
    cases <- c("TCGA-27-1831-01A-01D-0788-05","TCGA-S9-A6WP-01A-12D-A34D-05")
    query.met <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"),
                          legacy = TRUE,
                          data.category = "DNA methylation",
                          barcode = cases,
                          platform = c("Illumina Human Methylation 450", "Illumina Human Methylation 27"))
    GDCdownload(query.met, method = "api", directory = "ex")
    obj <- GDCprepare(query.met,  directory = "ex")
    expect_true(all(c("TCGA-LGG","TCGA-GBM") %in% SummarizedExperiment::colData(obj)$project_id))
    unlink("ex",recursive = TRUE, force = TRUE)
    unlink("MANIFEST.txt",recursive = TRUE, force = TRUE)
    unlink("Homo_sapiens_gene*",recursive = TRUE, force = TRUE)
})
