context("query_clinical")

test_that("TCGAquery_SampleTypes returns the correct barcodes", {
    sink("/dev/null");
    barcode <- c("TCGA-B0-4698-01Z-00-DX1","TCGA-CZ-4863-02Z-00-DX1")
    expect_equal(TCGAquery_SampleTypes(barcode,c("TR")),barcode[2])
    expect_equal(TCGAquery_SampleTypes(barcode,c("TP","TR")),barcode)
    expect_equal(TCGAquery_SampleTypes(barcode,c("TP")),barcode[1])
    expect_equal(TCGAquery_SampleTypes(barcode,c("TN")),"Error message: one or more sample types do not exist")
})

test_that("TCGAquery_MatchedCoupledSampleTypes returns the correct barcodes", {
    barcode <-c("TCGA-B0-4698-01Z-00-DX1",
                "TCGA-B0-4698-02Z-00-DX1",
                "TCGA-CD-4698-02Z-00-DX1",
                "TCGA-XO-4698-01Z-00-DX1",
                "TCGA-XO-4698-11Z-00-DX1")

    expect_equal(TCGAquery_MatchedCoupledSampleTypes(barcode, c("TP","TR")),barcode[1:2])
    expect_equal(TCGAquery_MatchedCoupledSampleTypes(barcode, c("TP","NT")),barcode[4:5])
    expect_equal(TCGAquery_MatchedCoupledSampleTypes(barcode,c("TN","TP")),"Error message: one or more sample types do not exist")
    expect_equal(TCGAquery_MatchedCoupledSampleTypes(barcode,c("TN")),"Error message: exactly two types need to be provided")

})

test_that("TCGAquery_subtype returns the a data frame if exists data", {
    expect_equal(class(TCGAquery_subtype("lgg")),class(data.frame()))
    expect_equal(class(TCGAquery_subtype("all")),class(data.frame()))
    expect_error(TCGAquery_subtype("f"))
})
