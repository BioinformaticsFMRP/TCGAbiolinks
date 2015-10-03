context("Query")

test_that("Results of query by center return only this center", {
    res <- TCGAquery(center = "broad.mit.edu")
    expect_that(length(unique(res$Center)), equals(1))
    expect_that(unique(res$Center), equals("broad.mit.edu"))
})

test_that("Results of query by platform return only this platform", {
    res <- TCGAquery(platform = "pathology_reports")
    expect_that(length(unique(res$Platform)), equals(1))
    expect_that(unique(res$Platform), equals("pathology_reports"))
})

test_that("Results of query by tumor return only this tumor", {
    res <- TCGAquery( tumor = "gbm")
    expect_that(unique(res$Disease), equals("GBM"))
    expect_that(length(unique(res$Disease)), equals(1))
    res <- TCGAquery( tumor = "GBM")
    expect_that(unique(res$Disease), equals("GBM"))
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
    expect_that(unique(res$Disease), equals("GBM"))
    expect_that(unique(res$Platform), equals("bio"))
    expect_that(length(unique(res$Disease)), equals(1))
    expect_that(length(unique(res$Platform)), equals(1))
    expect_identical(unique(res$Center),"nationwidechildrens.org")
    expect_that(length(unique(res$Center)), equals(1))
})

test_that("Results of query by center/platform/disease does not exist
          return empty data frame", {
              res <- TCGAquery(tumor = "gbm", platform = "bio",
                               center = "jhu-usc.edu")
              expect_that(nrow(res), equals(0))
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

