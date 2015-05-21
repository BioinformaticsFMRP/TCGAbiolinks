context("Query")

test_that("Results of query by center return only this center", {
    res <- TCGAQuery(center = "broad.mit.edu")
    expect_that(length(unique(res$Center)), equals(1))
    expect_that(unique(res$Center), equals("broad.mit.edu"))
})

test_that("Results of query by platform return only this platform", {
    res <- TCGAQuery(platform = "pathology_reports")
    expect_that(length(unique(res$Platform)), equals(1))
    expect_that(unique(res$Platform), equals("pathology_reports"))
})

test_that("Results of query by tumor return only this tumor", {
    res <- TCGAQuery( tumor = "gbm")
    expect_that(unique(res$Disease), equals("GBM"))
    expect_that(length(unique(res$Disease)), equals(1))
    res <- TCGAQuery( tumor = "GBM")
    expect_that(unique(res$Disease), equals("GBM"))
})

test_that("Results of query by barcode return this barcode", {
    sample <- "TCGA-06-0125-01A-01D-A45W-05"
    res <- TCGAQuery( sample = sample )
    expect_that(all(grepl(sample,res$barcode)), equals(TRUE))
})

test_that("Results of query by barcode return this barcode", {
    samples <- c("TCGA-06-0125-01A-01D-A45W-05","TCGA-06-0152-01A-02D-A45W-05")
    res <- TCGAQuery( sample = samples )
    expect_true(all(apply(sapply(samples,
                                 function(x) grepl(x,res$barcode)),1,any)))
})

test_that("Results of query by center/platform/disease return only them", {
    res <- TCGAQuery(tumor = "gbm", platform = "bio",
                     center = "nationwidechildrens.org")
    expect_that(unique(res$Disease), equals("GBM"))
    expect_that(unique(res$Platform), equals("bio"))
    expect_that(length(unique(res$Disease)), equals(1))
    expect_that(length(unique(res$Platform)), equals(1))
    expect_that(unique(res$Center), equals("nationwidechildrens.org"))
    expect_that(length(unique(res$Center)), equals(1))
})

test_that("Results of query by center/platform/disease does not exist
          return empty data frame", {
              res <- TCGAQuery(tumor = "gbm", platform = "bio",
                               center = "jhu-usc.edu")
              expect_that(nrow(res), equals(0))
          })

test_that("If platform argument is wrong, result is NULL ", {
    capture.output(res <- TCGAQuery( platform = "omg"))
    expect_that(res, equals(NULL))
})

test_that("If disease argument is wrong, result is NULL ", {
    capture.output(res <- TCGAQuery( tumor = "omg"))
    expect_that(res, equals(NULL))
})
test_that("If level argument is wrong, result is NULL ", {
    suppressMessages(res <- TCGAQuery( level = 4))
    expect_that(res, equals(NULL))
})

test_that("If center argument is wrong, result is NULL ", {
    capture.output(res <- TCGAQuery( center = "omg"))
    expect_that(res, equals(NULL))
})

