context("Download AND PREPARE")



test_that("GDCdownload API method is working ", {
    skip_on_bioc()
    skip_if_offline()

    cases <-  c("TCGA-OR-A5JX-01A-11R-A29S-07",
                "TCGA-OR-A5KY-01A-11R-A29S-07",
                "TCGA-PK-A5HA-01A-11R-A29S-07"
    )
    acc.gbm <- GDCquery(
        project =  c("TCGA-ACC"),
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        barcode = substr(cases,1,12)
    )
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
    skip_on_bioc()
    skip_if_offline()

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

    x <- getBarcodeInfo(
        c("HCM-CSHL-0063-C18-85A",
          "HCM-CSHL-0065-C20-06A",
          "HCM-CSHL-0065-C20-85A",
          "TARGET-20-PARUDL-03A",
          "TCGA-OR-A5LR-01A",
          "HCM-CSHL-0063-C18-01A")
    )
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","gender"] == "male")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","tumor_grade"] == "G2")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","ajcc_pathologic_stage"] == "Stage IVA")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0065-C20-06A","sample_type"] == "Metastatic")
    expect_true(x[x$sample_submitter_id == "HCM-CSHL-0063-C18-85A","sample_type"] == "Next Generation Cancer Model")

})

test_that("colDataPrepare handle replicates", {
    skip_on_bioc()
    skip_if_offline()
    barcodes <- c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01")
    x <- colDataPrepare(barcodes)
    expect_true(nrow(x) == 2)
    expect_true(all(rownames(x) == c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01")))
    expect_true(all(x$barcode == c("TCGA-06-0156-01A-02R-1849-01","TCGA-06-0156-01A-03R-1849-01")))
})

test_that("GDCprepare accepts more than one project", {
    skip_on_bioc()
    skip_if_offline()
    cases <-  c("TCGA-OR-A5JX-01A", "TCGA-OR-A5J3-01A",
                "TCGA-06-0680-11A","TCGA-14-0871-01A")
    expect_true(all(c("TCGA-ACC","TCGA-GBM") %in% colDataPrepare(cases)$project_id))
    acc.gbm <- GDCquery(project =  c("TCGA-ACC","TCGA-GBM"),
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        barcode = substr(cases,1,12))
    GDCdownload(acc.gbm, method = "api", directory = "ex")
    obj <- GDCprepare(acc.gbm,  directory = "ex")
    expect_true(all(c("TCGA-ACC","TCGA-GBM") %in% SummarizedExperiment::colData(obj)$project_id))
})

test_that("Non TCGA data is processed", {
    skip_on_bioc()
    skip_if_offline()
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
    #GDCdownload(query)
    #data <- GDCprepare(query)
})

test_that("GISTIC2 data is being correclty prepare", {
    skip_on_bioc()
    skip_if_offline()

    query <- GDCquery(
        project = "TCGA-COAD",
        data.category = "Copy Number Variation",
        data.type = "Gene Level Copy Number",
        access = "open",
        barcode = c("TCGA-AA-3522-10A","TCGA-A6-2672-10A")
    )
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
    skip_on_bioc()
    skip_if_offline()

    proj <- "TCGA-LUAD"
    query <- GDCquery(
        project = proj,
        data.category = "Raw microarray data",
        data.type = "Raw intensities",
        experimental.strategy = "Methylation array",
        legacy = TRUE,
        file.type = ".idat",
        barcode = "TCGA-55-7724",
        platform = "Illumina Human Methylation 450"
    )
    #tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
    #         error = function(e) GDCdownload(query, method = "client"))
    #betas <- GDCprepare(query)
    #expect_true(nrow(betas) == 485577)
    #expect_true(ncol(betas) == 1)
})

test_that("Prepare samples without clinical data", {
    skip_on_bioc()
    skip_if_offline()

    # x <-  GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
    # x[is.na(x$diagnosis_id),]
    x <- colDataPrepare(c("TCGA-80-5608-01A","TCGA-17-Z053-01A","TCGA-78-7158-01A"))
    expect_true(nrow(x) == 3)
})

test_that("Prepare multiple samples from the same patient", {
    skip_on_bioc()
    skip_if_offline()

    # https://portal.gdc.cancer.gov/cases/d7d3de82-802d-4664-8e42-d40408b129b0?bioId=548a300f-a7eb-4dc0-b9bc-5a643ef03d5d
    x <- colDataPrepare(c("BA2691R","BA2577R","BA2748R","BA2577D"))
    expect_true(nrow(x) == 4)
    expect_equal(x["BA2748R","sample_type"],"Primary Blood Derived Cancer - Bone Marrow")
    expect_equal(x["BA2577D","sample_type"],"Recurrent Blood Derived Cancer - Bone Marrow")
    expect_true("age_at_diagnosis" %in% colnames(x))
})

test_that("Preparing HT_HG-U133A as SE works", {
    skip_on_bioc()
    skip_if_offline()

    query <- GDCquery(
        project = "TCGA-GBM",
        legacy = TRUE,
        data.category = "Gene expression",
        data.type = "Gene expression quantification",
        platform = c("HT_HG-U133A")
    )
    query$results[[1]] <- query$results[[1]][1:2,]
    GDCdownload(query, method = "api", files.per.chunk = 100)
    se <- GDCprepare(query, summarizedExperiment = TRUE)

    expect_true(is(se,"SummarizedExperiment"))
})


test_that("Preparing RRPA files with number of proteins works", {
    skip_on_bioc()
    skip_if_offline()


    query_rppa <- GDCquery(
        project = c("TCGA-COAD"),
        data.category = "Proteome Profiling",
        experimental.strategy = "Reverse Phase Protein Array",
        platform = "RPPA",
        barcode = c("TCGA-CM-6165-01A","TCGA-DM-A28M-01A"),
        data.type = "Protein Expression Quantification"
    )

    GDCdownload(query_rppa)

    expect_message(object = {
        data_rppa <- GDCprepare(query_rppa)
    },regexp = "Some files differ in the number of proteins, we will introduce NA for the missing values")

    expect_true(is(data_rppa,"data.frame"))
})

