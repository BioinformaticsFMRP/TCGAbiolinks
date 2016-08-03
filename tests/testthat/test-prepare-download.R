context("Download AND PREPARE")

test_that("GDCdownload works", {
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

    query <- GDCquery(project = "TARGET-AML",
                      data.category = "Transcriptome Profiling",
                      data.type = "miRNA Expression Quantification",
                      workflow.type = "BCGSC miRNA Profiling",
                      barcode = c("TARGET-20-PARUDL-03A-01R","TARGET-20-PASRRB-03A-01R"))
    # data will be saved in  GDCdata/TARGET-AML/legacy/Copy_number_variation/Copy_number_segmentation
    #GDCdownload(query, method = "client", directory = "example_data_dir")
    GDCdownload(query, method = "api", directory = "example_data_dir")
    files <- file.path("example_data_dir/TARGET-AML/harmonized/",
                       gsub(" ","_",query$results[[1]]$data_category),
                       gsub(" ","_",query$results[[1]]$data_type),
                       query$results[[1]]$file_id,
                       query$results[[1]]$file_name)
    expect_true(all(file.exists(files)))
    expect_true(file.exists("GDCdata/TCGA-ACC/legacy/Copy_number_variation/Copy_number_segmentation/3f7804da-3731-4943-928c-bc6290af63d2/AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_D10_1348300.nocnv_hg19.seg.txt"))
    #unlink("GDCdata",recursive = TRUE, force = TRUE)
    #unlink("gdc-client*",recursive = TRUE, force = TRUE)
})
