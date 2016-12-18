context("Query")


test_that("GDCquery can filter by data.category", {
    query <- GDCquery(project = "TCGA-ACC",data.category = "Copy Number Variation")
    expect_equal(length(unique(query$results[[1]]$data_type)),2)
    query <- GDCquery(project = "TCGA-ACC",data.category = "Copy Number Variation", data.type = "Copy Number Segment")
    expect_equal(unique(query$results[[1]]$data_type),"Copy Number Segment")
})

test_that("GDCquery can filter by sample.type", {
    sample.type <- "Primary solid Tumor"
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy Number Variation",
                      data.type = "Masked Copy Number Segment",
                      sample.type = sample.type)
    expect_equal(as.character(unique(query$results[[1]]$tissue.definition)),sample.type)

    sample.type <- "Solid Tissue Normal"
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy Number Variation",
                      data.type = "Masked Copy Number Segment",
                      sample.type = sample.type)
    expect_equal(as.character(unique(query$results[[1]]$tissue.definition)),sample.type)

    sample.type <- c("Solid Tissue Normal", "Primary solid Tumor")
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy Number Variation",
                      data.type = "Masked Copy Number Segment",
                      sample.type = sample.type)
    expect_true(all(sample.type %in% unique(query$results[[1]]$tissue.definition)))
})

test_that("GDCquery can filter by barcode", {
    barcode <- c("TARGET-20-PADZCG-04A-01R","TARGET-20-PARJCR-09A-01R")
    query <- GDCquery(project = "TARGET-AML",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      barcode = barcode)
    expect_true(all(sort(barcode) == sort(unique(query$results[[1]]$cases))))
    barcode <- c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01")
    query <- GDCquery(project = "TCGA-ACC",
                      data.category = "Copy Number Variation",
                      data.type = "Copy Number Segment",
                      barcode = barcode)
    expect_true(all(sort(barcode) == sort(unique(query$results[[1]]$cases))))
    barcode <- c("TCGA-OR-A5KU", "TCGA-OR-A5JK")
    query <- GDCquery(project = "TCGA-ACC",
                      data.category = "Clinical",
                      barcode = barcode)
    expect_true(all(sort(barcode) == sort(unique(query$results[[1]]$cases))))

    # Will work if barcode was not found
    query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical",
                           barcode = c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"))
    expect_true(!all(c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q") %in% query$results[[1]]$cases))
})

test_that("GDCquery can filter copy number from legacy data by file type. Case: nocnv_hg18", {
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy number variation",
                      legacy = TRUE,
                      file.type = "nocnv_hg18.seg",
                      barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
    expect_equal(query$results[[1]]$file_name,"AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_D10_1348300.nocnv_hg18.seg.txt")
})

test_that("GDCquery can filter copy number from legacy data by file type. Case: hg18", {
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy number variation",
                      legacy = TRUE,
                      file.type = "hg18.seg",
                      barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
    expect_equal(query$results[[1]]$file_name,"AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_D10_1348300.hg18.seg.txt")
})

test_that("GDCquery can filter copy number from legacy data by file type. Case: hg19", {
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy number variation",
                      legacy = TRUE,
                      file.type = "hg19.seg",
                      barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
    expect_equal(query$results[[1]]$file_name,"AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_D10_1348300.hg19.seg.txt")
})


test_that("GDCquery can filter copy number from legacy data by file type. Case: nocnv_hg19", {
    query <- GDCquery(project = "TCGA-ACC",
                      data.category =  "Copy number variation",
                      legacy = TRUE,
                      file.type = "nocnv_hg19.seg",
                      barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
    expect_equal(query$results[[1]]$file_name,"AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_D10_1348300.nocnv_hg19.seg.txt")

})


test_that("GDCquery can filter by access level", {
    query <- GDCquery(project = "TCGA-KIRP",
                      data.category = "Simple Nucleotide Variation",
                      access = "open")
    expect_equal(unique(query$results[[1]]$access),"open")
    query <- GDCquery(project = "TCGA-KIRP",
                      data.category = "Simple Nucleotide Variation",
                      access = "controlled")
    expect_equal(unique(query$results[[1]]$access),"controlled")
})

test_that("GDCquery_Maf works", {
    acc.maf <- GDCquery_Maf("ACC",pipelines = "muse")
    expect_true(nrow(acc.maf) > 0)
    acc.maf <- GDCquery_Maf("ACC", directory = "maf", pipelines = "muse")
    expect_true(nrow(acc.maf) > 0)
    unlink("GDCdata",recursive = TRUE, force = TRUE)
    unlink("maf",recursive = TRUE, force = TRUE)
})
