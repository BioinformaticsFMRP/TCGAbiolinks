context("Analyse")

test_that("TCGAanalyze_survival creates pdf", {
    clin <- GDCquery_clinic("TCGA-ACC", type = "clinical", save.csv = FALSE)
    TCGAanalyze_survival(clin,clusterCol="gender",filename = "test.pdf")
    expect_true(file.exists("test.pdf"))
    unlink("test.pdf")
})


test_that("TCGAanalyze_DMR ask for the missing parameters", {
    nrows <- 2; ncols <- 20
    counts <- matrix(c(rep(0.9,20),rep(0.1,20)), nrows)
    rowRanges <- GenomicRanges::GRanges((rep("chr1",2)),
                                        IRanges::IRanges(c(2000,2000), width=100),
                                        strand=c("+","-"),
                                        feature_id=sprintf("ID%03d", 1:2))
    colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
                                    row.names=LETTERS[1:20],
                                    group=rep(c("group1","group2","group3","group4"),c(5,5,5,5)))
    data <- SummarizedExperiment::SummarizedExperiment(
        assays=S4Vectors::SimpleList(counts=counts),
        rowRanges=rowRanges,
        colData=colData)
    expect_null(TCGAanalyze_DMR(data, p.cut = 0.85))
    expect_message(TCGAanalyze_DMR(data, p.cut = 0.85),"Please, set the groupCol parameter")
    expect_null(TCGAanalyze_DMR(data, p.cut = 0.85,"group"))
    expect_message(TCGAanalyze_DMR(data, p.cut = 0.85,"group"),"Please, set the group1 and group2 parameters")
})

test_that("Results of TCGAanalyze_DEA inverting groups changes signal and order of the signals are right", {

    dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
    dataFilt <- TCGAanalyze_Filtering(tabDF = dataBRCA, method = "quantile", qnt.cut =  0.25)

    # 5 samples
    samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))

    # 5 samples
    samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

    # Get one line for example
    A <- rowMeans(dataFilt["CLDN6|9074",samplesNT])
    B <- rowMeans(dataFilt["CLDN6|9074",samplesTP])

    # Should give the same signal as  dataDEGs["CLDN6|9074",]
    log2FC <- log2(B) - log2(A)

    # Should give the same signal as  dataDEGs.inv["CLDN6|9074",]
    log2FC.inv <- log2(A) - log2(B)

    suppressMessages({
        dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                    mat2 = dataFilt[,samplesTP],
                                    Cond1type = "Normal",
                                    Cond2type =  "Tumor")
    })
    expect_equal(dataDEGs["CLDN6|9074",]$logFC > 0,(log2FC > 0)[[1]])

    suppressMessages({
        dataDEGs.inv <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesTP],
                                        mat2 = dataFilt[,samplesNT],
                                        Cond1type =  "Tumor",
                                        Cond2type = "Normal")
    })
    expect_equal(dataDEGs$logFC,-1*dataDEGs.inv$logFC)
    expect_equal(dataDEGs.inv["CLDN6|9074",]$logFC > 0,(log2FC.inv > 0)[[1]])
    suppressMessages({
        dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                    mat2 = dataFilt[,samplesTP],
                                    Cond1type = "Normal",
                                    Cond2type =  "Tumor", method = "glmLRT")
    })
    expect_equal(dataDEGs["CLDN6|9074",]$logFC > 0,(log2FC > 0)[[1]])

    suppressMessages({
        dataDEGs.inv <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesTP],
                                        mat2 = dataFilt[,samplesNT],
                                        Cond1type =  "Tumor",
                                        Cond2type = "Normal",
                                        method = "glmLRT")
    })
    expect_equal(dataDEGs$logFC,-1*dataDEGs.inv$logFC)
    expect_equal(dataDEGs.inv["CLDN6|9074",]$logFC > 0,(log2FC.inv > 0)[[1]])

})

test_that("Results from TCGAanalyze_DMR are correct", {
    nrows <- 2; ncols <- 20
    counts <- matrix(c(rep(0.9,20),rep(0.1,20)), nrows)
    rowRanges <- GenomicRanges::GRanges((rep("chr1",2)),
                                        IRanges::IRanges(c(2000,2000), width=100),
                                        strand=c("+","-"),
                                        feature_id=sprintf("ID%03d", 1:2))
    colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
                                    row.names=LETTERS[1:20],
                                    group=rep(c("group1","group2"),c(10,10)))
    data <- SummarizedExperiment::SummarizedExperiment(
        assays=S4Vectors::SimpleList(counts=counts),
        rowRanges=rowRanges,
        colData=colData)
    SummarizedExperiment::colData(data)$group <- c(rep("group1",10),  rep("group2",10))
    hypo.hyper <- TCGAanalyze_DMR(data, p.cut = 0.85,"group","group1","group2")
    result <- values(hypo.hyper)[1,]
    expect_equal(result$mean.group1,0.9)
    expect_equal(result$mean.group2,0.1)
    expect_equal(result$diffmean.group1.group2 , -0.8)
    expect_equal(result$status.group1.group2 , "Hypomethylated")
    expect_equal(result$diffmean.group2.group1 , 0.8)
    expect_equal(result$status.group2.group1 , "Hypermethylated")

})


test_that("Results from TCGAanalyze_DEA and DMR in starburst plot are correct", {
    nrows <- 4; ncols <- 20
    counts <- matrix(nrow = nrows,ncol = ncols)
    counts[1,] <-c(rep(0.9,ncols/2),rep(0.1,ncols/2))
    counts[2,] <-c(rep(0.3,ncols/2),rep(0.1,ncols/2))
    counts[3,] <-c(rep(0.1,ncols/2),rep(0.9,ncols/2))
    counts[4,] <-c(rep(0.1,ncols/2),rep(0.3,ncols/2))

    rowRanges <- GenomicRanges::GRanges((rep("chr1",nrows)),
                                        IRanges::IRanges(rep(2000,nrows), width=100),
                                        strand=rep("+",nrows),
                                        feature_id=sprintf("ID%03d", 1:nrows),
                                        Gene_Symbol = sprintf("Gene%03d", 1:4))
    colData <- S4Vectors::DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
                                    row.names=LETTERS[1:20],
                                    group=rep(c("group1","group2"),c(10,10)))
    data <- SummarizedExperiment::SummarizedExperiment(
        assays=S4Vectors::SimpleList(counts=counts),
        rowRanges=rowRanges,
        colData=colData)
    SummarizedExperiment::colData(data)$group <- c(rep("group1",10),  rep("group2",10))
    met <- TCGAanalyze_DMR(data, p.cut = 0.85,"group","group1","group2")

    # Expression results
    #  logFC     FDR
    #    2     0.00001
    #    2       0.1
    #   -3       0.1
    #   -3    0.00001
    exp <- data.frame(logFC = c(2,2,-3,-3),
                      FDR = c(0.00001,0.1,0.00001,0.1))
    rownames(exp) <- sprintf("Gene%03d", 1:4)

    result.no.cut <- TCGAvisualize_starburst(met,exp,
                                             exp.p.cut = 1, met.p.cut = 1,
                                             group1="group1",group2="group2",
                                             diffmean.cut=0.0,logFC.cut = 0,
                                             names=TRUE, circle = FALSE,return.plot = TRUE)$starburst

    result.fdr.cut <- TCGAvisualize_starburst(met,exp,
                                              exp.p.cut = 0.05, met.p.cut = 0.05,
                                              group1="group1",group2="group2",
                                              diffmean.cut=0.0,
                                              names=TRUE, circle = FALSE,return.plot = TRUE)$starburst
    result.fc.cut <- TCGAvisualize_starburst(met,exp,
                                             exp.p.cut = 1, met.p.cut = 1,
                                             group1="group1",group2="group2",
                                             diffmean.cut=0.0,logFC.cut = 2.5,
                                             names=TRUE, circle = FALSE,return.plot = TRUE)$starburst

    result.met.cut <- TCGAvisualize_starburst(met,exp,
                                              exp.p.cut = 1, met.p.cut = 1,
                                              group1="group1",group2="group2",
                                              diffmean.cut=0.5,logFC.cut = 0,
                                              names=TRUE, circle = FALSE,return.plot = TRUE)$starburst
    result.met.exp.cut <- TCGAvisualize_starburst(met,exp,
                                                  exp.p.cut = 1, met.p.cut = 1,
                                                  group1="group1",group2="group2",
                                                  diffmean.cut=0.5,logFC.cut = 2.5,
                                                  names=TRUE, circle = FALSE,return.plot = TRUE)$starburst


    # Threshold are respected
    expect_equal(nrow(result.fdr.cut), 2)
    expect_equal(nrow(result.no.cut), 4)
    expect_equal(nrow(result.fc.cut), 2)
    expect_equal(nrow(result.met.cut), 2)
    expect_equal(nrow(result.met.exp.cut), 1)

    # group1 vs groups2 (logFC = log(group2) - log(group1), diffmean = group2 - group1 )
    expect_true(result.met.exp.cut$starburst.status == "Down regulated & Hyper methylated" &
                    result.met.exp.cut$logFC < 0 & result.met.exp.cut$diffmean.group1.group2 > 0)
    expect_true(result.met.cut[1,]$starburst.status == "Up regulated & Hypo methylated" &
                    result.met.cut[1,]$logFC > 0 & result.met.cut[1,]$diffmean.group1.group2 < 0)
    expect_true(result.met.cut[2,]$starburst.status == "Down regulated & Hyper methylated" &
                    result.met.cut[2,]$logFC < 0 & result.met.cut[2,]$diffmean.group1.group2 > 0)

    # --- Changing groups order
    result.met.cut.inv <- TCGAvisualize_starburst(met,exp,
                                                  exp.p.cut = 1, met.p.cut = 1,
                                                  group1="group2",group2="group1",
                                                  diffmean.cut=0.5,logFC.cut = 0,
                                                  names=TRUE, circle = FALSE,return.plot = TRUE)$starburst
    result.met.exp.cut.inv <- TCGAvisualize_starburst(met,exp,
                                                      exp.p.cut = 1, met.p.cut = 1,
                                                      group1="group2",group2="group1",
                                                      diffmean.cut=0.5,logFC.cut = 2.5,
                                                      names=TRUE, circle = FALSE,return.plot = TRUE)$starburst

    # group2 vs groups1 (logFC = log(group1) - log(group2), diffmean = group1 - group2 )
    expect_true(result.met.exp.cut.inv$starburst.status == "Down regulated & Hypo methylated" &
                    result.met.exp.cut.inv$logFC < 0 & result.met.exp.cut.inv$diffmean.group2.group1 < 0)
    expect_true(result.met.cut.inv[1,]$starburst.status == "Up regulated & Hyper methylated" &
                    result.met.cut.inv[1,]$logFC > 0 & result.met.cut.inv[1,]$diffmean.group2.group1 > 0)
    expect_true(result.met.cut.inv[2,]$starburst.status == "Down regulated & Hypo methylated" &
                    result.met.cut.inv[2,]$logFC < 0 & result.met.cut.inv[2,]$diffmean.group2.group1 < 0)
    unlink("DMR_results_group_group1_group2_pcut_0.85_meancut_0.2.csv")
    unlink("group_group1_group2_pcut_0.85_meancut_0.2.rda")
    unlink("histogram_diffmean.png")
    unlink("histogram_pvalues.png")
    unlink("histogram_pvalues_adj.png")
    unlink("methylation_volcano.pdf")
})
