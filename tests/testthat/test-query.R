context("Query")

test_that("Results of query by center return only this center", {
    res <- TCGAquery(center = "broad.mit.edu")
    expect_equal(length(unique(res$Center)),1)
    expect_equal(unique(res$Center), "broad.mit.edu")
})

test_that("Results of query by platform return only this platform", {
    res <- TCGAquery(platform = "pathology_reports")
    expect_equal(length(unique(res$Platform)), 1)
    expect_equal(unique(res$Platform), "pathology_reports")
})

test_that("Results of query by tumor return only this tumor", {
    res <- TCGAquery( tumor = "gbm")
    expect_identical(unique(res$Disease), "GBM")
    expect_equal(length(unique(res$Disease)),1)
    res <- TCGAquery( tumor = "GBM")
    expect_identical(unique(res$Disease), "GBM")
})

test_that("Results of query by barcode return this barcode", {
    sample <- "TCGA-06-0125-01A-01D-A45W-05"
    res <- TCGAquery( sample = sample )
    expect_true(all(grepl(sample,res$barcode)))
})

test_that("Results of query by barcode return this barcode", {
    samples <- c("TCGA-06-0125-01A-01D-A45W-05","TCGA-06-0152-01A-02D-A45W-05")
    res <- TCGAquery( sample = samples )
    expect_true(all(apply(sapply(samples,
                                 function(x) grepl(x,res$barcode)),1,any)))
})

test_that("Results of query by center/platform/disease return only them", {
    res <- TCGAquery(tumor = "gbm", platform = "bio",
                     center = "nationwidechildrens.org")
    expect_identical(unique(res$Disease), "GBM")
    expect_identical(unique(res$Platform), "bio")
    expect_equal(length(unique(res$Disease)), 1)
    expect_equal(length(unique(res$Platform)), 1)
    expect_identical(unique(res$Center),"nationwidechildrens.org")
    expect_equal(length(unique(res$Center)), 1)
})

test_that("Results of query by center/platform/disease does not exist
          return empty data frame", {
              res <- TCGAquery(tumor = "gbm", platform = "bio",
                               center = "jhu-usc.edu")
              expect_equal(nrow(res), 0)
          })

test_that("If platform argument is wrong, result is NULL ", {
    capture.output(res <- TCGAquery( platform = "omg"))
    expect_null(res)
})

test_that("If disease argument is wrong, result is NULL ", {
    capture.output(res <- TCGAquery( tumor = "omg"))
    expect_null(res)
})
test_that("If level argument is wrong, result is NULL ", {
    suppressMessages(res <- TCGAquery( level = 4))
    expect_null(res)
})

test_that("If center argument is wrong, result is NULL ", {
    capture.output(res <- TCGAquery( center = "omg"))
    expect_null(res)
})

test_that("If the version specified it is updated", {
    res <- TCGAquery(tumor = c("gbm"),
                     platform = c("HumanMethylation27"),
                     level = 3,
                     samples = "TCGA-14-4157-01A-01D-1228-05",
                     version = list(c("HumanMethylation27","GBM",5)))
    expect_true( grepl("3.9.5.0",res$name))
})


test_that("If no result exists return is NULL ", {
    capture.output(res <- TCGAquery(tumor = "gbm",platform = "bio",level=3))
    expect_null(res)
})

test_that("TCGAquery_maf works", {
    capture.output(mutation <- TCGAquery_maf(tumor = "lgg",
                               archive.name = "IlluminaGA_DNASeq_curated.Level_2.1.4.0"))
    expect_equal(class(mutation),class(data.frame()))
    expect_null(TCGAquery_maf(tumor = "lgg", archive.name = "xxx"))
    capture.output(mutation <- TCGAquery_maf(tumor = "THYM",
                               archive.name = "genome.wustl.edu_THYM.IlluminaHiSeq_DNASeq_automated.Level_2.1.1.0"))
    expect_equal(class(mutation),class(data.frame()))

})
