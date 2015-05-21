context("Download")

test_that("It can download a file", {
    sample <- "TCGA-33-6737-01"
    suppressMessages(
        query <- TCGAQuery(tumor = "lusc", samples = sample,
                           platform = "IlluminaHiSeq_RNASeqV2", level = "3")
    )
    suppressMessages(
        out <- TCGADownload(query, path = "dataTest",
                            type = "rsem.genes.results",
                            samples = sample, quiet = TRUE)
    )
    folder <- gsub(".tar.gz","",basename(query$deployLocation))
    file <- file.path("dataTest",folder,
                      paste0("unc.edu.38dbab79-7aeb-4153-",
                             "b3ce-8feeb36b48e5.1096372.rsem.genes.results"))
    expect_true(file.exists(file))
    unlink("dataTest", recursive = TRUE)
    unlink("mages", recursive = TRUE)
})

test_that("It can download a folder", {
    query <- TCGAQuery(tumor = "lusc", level = "mage-tab",
                       platform = "IlluminaHiSeq_RNASeqV2")

    suppressMessages(
        out <- TCGADownload(query, path = "dataTest", quiet = TRUE)
    )
    folder <- basename(query$deployLocation)
    file <- file.path("dataTest",folder)
    expect_true(file.exists(file))
    unlink("dataTest", recursive = TRUE)
    unlink("mages", recursive = TRUE)
})
