context("Download AND PREPARE")

test_that("It can download a file", {
    sample <- c("TCGA-33-6737-01","TCGA-33-6737-11")
    query <- TCGAquery(tumor = "lusc", samples = sample,
                       platform = "IlluminaHiSeq_RNASeqV2", level = "3")
    capture.output(out <- TCGAdownload(query, path = "dataTest",
                                       type = "rsem.genes.results",
                                       samples = sample))
    folder <- gsub(".tar.gz","",basename(query$deployLocation))
    file <- file.path("dataTest",folder,
                      paste0("unc.edu.38dbab79-7aeb-4153-",
                             "b3ce-8feeb36b48e5.1096372.rsem.genes.results"))
    expect_true(file.exists(file))
})

test_that("It can download a folder", {
    query <- TCGAquery(tumor = "lusc", level = "mage-tab",
                       platform = "IlluminaHiSeq_RNASeqV2")

    capture.output(
        out <- TCGAdownload(query, path = "dataTest")
    )
    folder <- basename(query$deployLocation)
    file <- file.path("dataTest",folder)
    expect_true(file.exists(file))
})

test_that("It can download and prepare using types", {
    listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07")

    type <- c("junction_quantification","rsem.genes.results",
              "rsem.isoforms.results", "rsem.genes.normalized_results",
              "rsem.isoforms.normalized_results", "bt.exon_quantification")

    for(i in type){
        message(i)

        # Query platform IlluminaHiSeq_RNASeqV2 with a list of barcode
        query <- TCGAquery(tumor = "brca", samples = listSamples,
                           platform = "IlluminaHiSeq_RNASeqV2", level = "3")
        capture.output(
            # Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
            TCGAdownload(query, path = "dataBrca", type = i, samples = listSamples)
        )
        # Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
        # rsem.genes.results as values
        capture.output(
            df <- TCGAprepare(query,"dataBrca", type = i, summarizedExperiment = TRUE)
        )
        expect_true(any(grepl(i,dir("dataBrca",full.names = TRUE,recursive = T))))
        expect_false(is.null(class(df)))
    }
})

