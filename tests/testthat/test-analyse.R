context("Analyse")

test_that("TCGAanalyze_survival creates pdf", {
    clin <- data.frame(
           vital_status = c("alive","alive","alive","dead","alive","alive","dead","alive","dead","alive"),
           days_to_death = c(NA,NA,NA,172,NA,NA,3472,NA,786,NA),
           days_to_last_follow_up = c(3011,965,718,NA,1914,423,NA,5,656,1417),
           gender = c(rep("male",5),rep("female",5))
    )
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

test_that("TCGAanalyze_DMR is handling NAs correctly", {
    nrows <- 2; ncols <- 20
    counts <- matrix(c(rep(0.9,20),rep(0.1,20)), nrows)
    counts[1,1] <- NA
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

    counts[1,] <- NA
    data <- SummarizedExperiment::SummarizedExperiment(
        assays=S4Vectors::SimpleList(counts=counts),
        rowRanges=rowRanges,
        colData=colData)
    expect_error(TCGAanalyze_DMR(data, p.cut = 0.85,"group","group1","group2"),
                   "Sorry, but we found some probes with NA for all samples in your data, please either remove/or replace them")
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
