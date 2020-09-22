context("query_clinical")

test_that("TCGAquery_SampleTypes returns the correct barcodes", {
    barcode <- c("TCGA-B0-4698-01Z-00-DX1","TCGA-CZ-4863-02Z-00-DX1")
    expect_equal(TCGAquery_SampleTypes(barcode,c("TR")),barcode[2])
    expect_equal(TCGAquery_SampleTypes(barcode,c("TP","TR")),barcode)
    expect_equal(TCGAquery_SampleTypes(barcode,c("TP")),barcode[1])
    expect_equal(TCGAquery_SampleTypes(barcode,c("TN")),"Error message: one or more sample types do not exist")
})


test_that("GDCquery_clinic populates correctly the data", {
    results <- GDCquery_clinic( "BEATAML1.0-COHORT")
    results.2028 <- results[results$submitter_id == "2028",]
    expect_equal(results.2028$vital_status,"Alive")
    expect_true(all(c("BA2486D","BA2144D") %in%
                        (str_split(results.2028$submitter_sample_ids,",") %>% unlist())))
    expect_equal(results.2028$age_at_diagnosis %% 365.25,134)
    expect_equal(as.integer(results.2028$age_at_diagnosis / 365.25),56)

    results <- GDCquery_clinic(project = "ORGANOID-PANCREATIC", type = "clinical")
    results.42 <- results[results$submitter_id == "42",]
    expect_equal(results.42$tumor_stage,"2A")
    expect_equal(results.42$ethnicity,"not hispanic or latino")
    expect_equal(as.integer(results.2028$age_at_diagnosis / 365.25),56)

    results <- GDCquery_clinic( "TCGA-LUAD")
    results.sample <- results[results$submitter_id == "TCGA-80-5608",]
    expect_equal(results.sample$vital_status,"Alive")
    expect_equal(results.sample$gender,"female")
    expect_equal(results.sample$ajcc_pathologic_stage,"Stage IA")

    # Sample without clinical data
    results.sample <- results[results$submitter_id == "TCGA-17-Z031",]
    expect_true(is.na(results.sample$vital_status))
    expect_true(is.na(results.sample$gender))
    expect_true(is.na(results.sample$ajcc_pathologic_stage))
})


test_that("TCGAquery_MatchedCoupledSampleTypes returns the correct barcodes", {
    barcode <- c("TCGA-B0-4698-01Z-00-DX1",
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


