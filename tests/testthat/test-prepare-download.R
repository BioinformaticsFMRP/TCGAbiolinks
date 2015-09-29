context("Download AND PREPARE")

test_that("It can download a file", {
    sample <- c("TCGA-33-6737-01","TCGA-33-6737-11")
        query <- TCGAquery(tumor = "lusc", samples = sample,
                           platform = "IlluminaHiSeq_RNASeqV2", level = "3")
        out <- TCGAdownload(query, path = "dataTest",
                            type = "rsem.genes.results",
                            samples = sample)
    folder <- gsub(".tar.gz","",basename(query$deployLocation))
    file <- file.path("dataTest",folder,
                      paste0("unc.edu.38dbab79-7aeb-4153-",
                             "b3ce-8feeb36b48e5.1096372.rsem.genes.results"))
    expect_true(file.exists(file))
})

test_that("It can download a folder", {
    query <- TCGAquery(tumor = "lusc", level = "mage-tab",
                       platform = "IlluminaHiSeq_RNASeqV2")

    suppressMessages(
        out <- TCGAdownload(query, path = "dataTest")
    )
    folder <- basename(query$deployLocation)
    file <- file.path("dataTest",folder)
    expect_true(file.exists(file))
})


