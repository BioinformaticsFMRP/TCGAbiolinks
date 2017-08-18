context("Visualize")

test_that("EAbarplot works", {
    Genelist <- c("FN1","COL1A1")
    ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
    TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                            GOBPTab = ansEA$ResBP,
                            GOCCTab = ansEA$ResCC,
                            GOMFTab = ansEA$ResMF,
                            PathTab = ansEA$ResPat,
                            nRGTab = Genelist,
                            nBar = 10,
                            filename="a.pdf")
    expect_true(file.exists("a.pdf"))
    unlink("a.pdf")
})

test_that("TCGAvisualize_meanMethylation works", {
    query <- GDCquery(project = "TCGA-GBM",
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 27",
                      legacy = TRUE,
                      barcode = c("TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
                                  "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
                                  "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
                                  "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"))
    tryCatch(GDCdownload(query, method = "api",chunks.per.download = 20),
             error = function(e) GDCdownload(query, method = "api"))
    data <- GDCprepare(query)

    # setting y limits: lower 0, upper 1
    TCGAvisualize_meanMethylation(data,groupCol = "shortLetterCode",filename = "tcga-gbm-1.pdf", y.limits = c(0,1))
    # setting y limits: lower 0
    TCGAvisualize_meanMethylation(data,groupCol = "shortLetterCode",filename = "tcga-gbm-2.pdf", y.limits = 0)
    # Changing shapes of jitters to show subgroups
    TCGAvisualize_meanMethylation(data,
                                  groupCol  = "subtype_Pan.Glioma.DNA.Methylation.Cluster",
                                  subgroupCol="vital_status",
                                  filename = "tcga-gbm-3.pdf")
    # Sorting bars by descending mean methylation
    TCGAvisualize_meanMethylation(data,
                                  groupCol  = "subtype_Pan.Glioma.DNA.Methylation.Cluster",
                                  sort="mean.desc",
                                  filename="tcga-gbm-4.pdf")
    # Sorting bars by asc mean methylation
    TCGAvisualize_meanMethylation(data,groupCol  = "subtype_Pan.Glioma.DNA.Methylation.Cluster",
                                  sort="mean.asc",
                                  filename=NULL)
    TCGAvisualize_meanMethylation(data,
                                  groupCol  = "vital_status",
                                  sort="mean.asc",
                                  filename="tcga-gbm-5.pdf")

    expect_true(file.exists("tcga-gbm-1.pdf"))
    expect_true(file.exists("tcga-gbm-2.pdf"))
    expect_true(file.exists("tcga-gbm-3.pdf"))
    expect_true(file.exists("tcga-gbm-4.pdf"))
    expect_true(file.exists("tcga-gbm-5.pdf"))
    unlink("tcga-gbm-*")
})

#test_that("TCGAvisualize_oncoprint works", {
#    mut <- GDCquery_Maf(tumor = "ACC",pipelines = "muse")
#    clin <- GDCquery_clinic("TCGA-ACC","clinical")
#    clin <- clin[,c("bcr_patient_barcode","disease","gender","tumor_stage","race","vital_status")]
#    TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:20],
#                            filename = "oncoprint.pdf",
#                            annotation = clin,
#                            color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
#                            rows.font.size=10,
#                            heatmap.legend.side = "right",
#                            dist.col = 0,
#                            label.font.size = 10)
#    unlink("oncoprint.pdf")
#
#})


